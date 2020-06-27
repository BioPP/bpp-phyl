//
// File: SubstitutionMappingToolsForASite.cpp
// Created by: Laurent Guéguen
// Created on: mercredi 26 septembre 2018, à 15h 11
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#include "SubstitutionMappingToolsForASite.h"
#include "DecompositionSubstitutionCount.h"
#include "DecompositionReward.h"
#include "ProbabilisticRewardMapping.h"
#include "ProbabilisticSubstitutionMapping.h"
#include "RewardMappingToolsForASite.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>

using namespace bpp;
using namespace std;

SubstitutionMappingToolsForASite::t_Sr_Sm_Sc SubstitutionMappingToolsForASite::m_Sr_Sm_Sc;

// From the STL:
#include <iomanip>

/******************************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingToolsForASite::computeCounts(
  size_t site,
  LikelihoodCalculationSingleProcess& rltc,
  const vector<uint>& edgeIds,
  const SubstitutionRegister& reg,
  std::shared_ptr<const AlphabetIndex2> weights,
  std::shared_ptr<const AlphabetIndex2> distances,
  double threshold,
  bool verbose)
{
  
  // Preamble:
  if (!rltc.isInitialized())
    throw Exception("SubstitutionMappingToolsForASite::computeSubstitutionVectors(). Likelihood object is not initialized.");

  const SubstitutionProcess& sp=rltc.getSubstitutionProcess();

  if (edgeIds.size()==0)
  {
    return new ProbabilisticSubstitutionMapping(sp.getParametrizablePhyloTree(), reg.getNumberOfSubstitutionTypes(), 0);
  }
  
  // Map from models to counts
  if (m_Sr_Sm_Sc.find(&reg)==m_Sr_Sm_Sc.end())
    m_Sr_Sm_Sc[&reg]=t_Sm_Sc();
  t_Sm_Sc& m_Sm_Sc=m_Sr_Sm_Sc[&reg];
  
  auto processTree = rltc.getTreeNode(0);

  for (auto speciesId :edgeIds)
  {
    const auto& dagIndexes = rltc.getEdgesIds(speciesId, 0);

    for (auto id:dagIndexes)
    {
      const auto& edge = processTree->getEdge(id);
      if (edge->getBrLen()) // if edge with model on it
      {
        auto model = edge->getModel();

        auto nMod = edge->getNMod();
        
        auto tm = dynamic_cast<const TransitionModel*>(model->getTargetValue());

        const SubstitutionModel* sm(0);
        
        if (nMod==0)
        {
          sm = dynamic_cast<const SubstitutionModel*>(tm);
          
          if (sm == NULL)
            throw Exception("SubstitutionMappingTools::computeCounts : SubstitutionVectors possible only for SubstitutionModels, not in branch " + TextTools::toString(speciesId) + ". Got model " + tm->getName());
        }
        else
        {
          size_t nmod = nMod->getTargetValue();

          auto ttm = dynamic_cast<const MixedTransitionModel*>(tm);
          if (ttm == NULL)
            throw Exception("SubstitutionMappingTools::computeCounts : Expecting Mixed model in branch " + TextTools::toString(speciesId) + ". Got model " + tm->getName());

          sm = dynamic_cast<const SubstitutionModel*>(ttm->getNModel(nmod));
          
          if (sm==NULL)
            throw Exception("SubstitutionMappingTools::computeCounts : Expecting Substitution model for submodel " + TextTools::toString(nmod) + " of mixed model " + tm->getName() + " in branch " + TextTools::toString(speciesId));
        }
            
        if (m_Sm_Sc.find(sm)==m_Sm_Sc.end())
        {
          m_Sm_Sc[sm] = std::make_shared<DecompositionSubstitutionCount>(sm, reg.clone(), weights, distances);
        }
      }
    }
  }

  // We create a Mapping objects

  if (m_Sm_Sc.size()==0)
    throw Exception("SubstitutionMappingToolsForASite::computeSubstitutionVectors not possible with null model.");

  return computeCounts(site, rltc, m_Sm_Sc, edgeIds, threshold, verbose);
}


ProbabilisticSubstitutionMapping* SubstitutionMappingToolsForASite::computeCounts(
  size_t site,
  LikelihoodCalculationSingleProcess& rltc,
  t_Sm_Sc& m_Sm_Sc,
  const vector<uint>& edgeIds,
  double threshold,
  bool verbose)
{
  return 0;
  /*
  if (!rltc.isInitialized())
    throw Exception("SubstitutionMappingToolsForASite::computeSubstitutionVectors(). Likelihood object is not initialized.");
  
  const SubstitutionProcess& sp=rltc.getSubstitutionProcess();

  if (edgeIds.size()==0)
  {
    return new ProbabilisticSubstitutionMapping(sp.getParametrizablePhyloTree();, m_Sm_Sc.begin()->second->getNumberOfSubstitutionTypes(), 0);
  }

  ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();
  
  // A few variables we'll need:

  size_t nbStates        = sp.getNumberOfStates();
  size_t nbClasses       = sp.getNumberOfClasses();

  size_t nbTypes         = m_Sm_Sc.begin()->second->getNumberOfSubstitutionTypes();
  size_t nbNodes         = edgeIds.size();
  
  // We create a Mapping objects
  
  unique_ptr<ProbabilisticSubstitutionMapping> substitutions(new ProbabilisticSubstitutionMapping(ppt, nbTypes, 1));

  // Find the corresponding index in the compressed array :

  size_t siteIndex=rltc.getRootArrayPosition(site);
  
  double Lr=rltc.getLikelihoodForASite(site);

  // Compute the number of substitutions for each class and each branch in the tree:
  if (verbose)
    ApplicationTools::displayTask("Compute counts", true);

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt=substitutions->allEdgesIterator();

  size_t nn=0;
  for (;!brIt->end();brIt->next())
  {
    if (verbose)
      ApplicationTools::displayGauge(nn++, nbNodes - 1);
    
    shared_ptr<PhyloBranchMapping> br=**brIt;
    
    // For each branch
    uint speciesId = substitutions->getEdgeIndex(br);

    if (edgeIds.size() > 0 && !VectorTools::contains(edgeIds, (int)speciesId))
      continue;

    vector<double> substitutionsForCurrentNode(nbTypes);

    for (size_t ncl=0; ncl<nbClasses; ncl++)
    {
      processTree = rltc.getTreeNode(ncl);
      double pr = sp.getProbabilityForModel(ncl);

      vector<double> substitutionsForCurrentClass(nbTypes);
      
      const auto& dagIndexes = rltc.getEdgesIds(speciesId, ncl);

      Eigen::MatrixXd npxy;

      // Sum on all dag edges for this speciesId
      for (auto id:dagIndexes)
      {
        auto edge = processTree->getEdge(id);

        auto tm = dynamic_cast<const TransitionModel*>(edge->getModel()->getTargetValue());

        auto nMod = edge->getNMod();
        
        const SubstitutionModel* sm(0);
        
        if (nMod==0)
          sm = dynamic_cast<const SubstitutionModel*>(tm);
        else
        {
          size_t nmod = nMod->getTargetValue();

          auto ttm = dynamic_cast<const MixedTransitionModel*>(tm);
          sm = dynamic_cast<const SubstitutionModel*>(ttm->getNModel(nmod));
        }

        auto subCount = mModCount[sm];

        const auto& likelihoodsFather = rltc.getBackwardLikelihoodsAtEdgeForClass(id, ncl)->getTargetValue();

        auto sonid = rltc.getForwardLikelihoodTree(ncl)->getSon(id);

        const auto& likelihoodsSon = rltc.getForwardLikelihoodsAtNodeForClass(sonid, ncl)->getTargetValue();

        const Eigen::MatrixXd& pxy = edge->getTransitionMatrix()->getTargetValue();

        for (size_t t = 0; t < nbTypes; ++t)
        {
          // compute all nxy * pxy first:

          subCount->storeAllNumbersOfSubstitutions(edge->getBrLen()->getValue(), t + 1, npxy);
          
          npxy.array() *= pxy.array();
          
          // Now loop over sites:

          auto counts = npxy * likelihoodsSon;

          auto bb = (likelihoodsFather.cwiseProduct(counts)).colwise().sum();

          substitutionsForCurrentClass[t] += bb;
        }
      }
      
      for (size_t t = 0; t < nbTypes; ++t)
        substitutionsForCurrentNode[t] += substitutionsForCurrentClass[t] * pr;
    }
    if (verbose)
      ApplicationTools::displayGauge(nn++, nbNodes - 1);
    
    shared_ptr<PhyloBranchMapping> br=**brIt;
    
    // For each branch
    uint speciesId = substitutions->getEdgeIndex(br);

    if (edgeIds.size() > 0 && !VectorTools::contains(edgeIds, (int)edid))
      continue;

    for (size_t ncl=0; ncl<nbClasses; ncl++)
    {
      processTree = rltc.getTreeNode(ncl);
      double pr = sp.getProbabilityForModel(ncl);


      // Then, we deal with the node of interest.
      // We first average upon 'y' to save computations, and then upon 'x'.
      // ('y' is the state at 'node' and 'x' the state at 'father'.)

      const Vdouble& likelihoodsFather_node = ici->getBelowLikelihoodArray(ComputingNode::D0)[siteIndex];

      const SubstitutionModel* sm=dynamic_cast<const SubstitutionModel*>(sp.getModel(icid,ncl));
      if (!sm)
        throw Exception("SubstitutionMappingToolsForASite:: non substitution model in node " + TextTools::toString(icid));

      if (m_Sm_Sc.find(sm)==m_Sm_Sc.end())
        throw Exception("SubstitutionMappingToolsForASite::computeCounts : Unknown model for substitution count. Should not happen.");
        
      SubstitutionCount& psc=*m_Sm_Sc.find(sm)->second;
      
      // compute all nxy * pxy first:
      
      const Matrix<double>& pxy = sp.getTransitionProbabilities(icid,ncl);
      
      for (size_t t = 0; t < nbTypes; ++t)
      {
        Matrix<double>* nijt = psc.getAllNumbersOfSubstitutions(d * rate, t + 1);
        MatrixTools::hadamardMult((*nijt),pxy,npxy[t]);        
        delete nijt;
      }
      
      // Now compute:
      
      for (size_t x = 0; x < nbStates; ++x)
      {
        double likelihoodsFatherConstantPart_x = likelihoodsFatherConstantPart[x];
        for (size_t y = 0; y < nbStates; ++y)
        {
          double likelihood_xy = usesLog
            ?likelihoodsFatherConstantPart_x + likelihoodsFather_node[y]
            :likelihoodsFatherConstantPart_x * likelihoodsFather_node[y];
          
          if (!usesLog)
          {
            if (likelihood_xy!=0) // to avoid multiplication per nan
              // (stop codons)
              for (size_t t = 0; t < nbTypes; ++t)
              {
                substitutionsForCurrentNode[t] += likelihood_xy * npxy[t](x,y);
                //                                <------------>  <----------->
                // Posterior probability              |               |
                // likelihood               ----------+               |
                //                                                    |
                // Substitution function            ------------------+
              }
          }
          else
          {
            if (likelihood_xy!=NumConstants::MINF())  // to avoid add per -inf
              for (size_t t = 0; t < nbTypes; ++t)
              {
                if (npxy[t](x,y)< -NumConstants::MILLI())
                {
                  ApplicationTools::displayWarning("These counts are negative, their logs could not be computed:" + TextTools::toString(npxy[t](x,y)));
                  throw Exception("Stop in SubstitutionMappingToolsForASite");
                }
                else
                {
                  if (npxy[t](x,y)>0)
                  {
                    double ll=likelihood_xy + log(npxy[t](x,y));
                    if (ll>substitutionsForCurrentNode[t])
                      substitutionsForCurrentNode[t] = ll + log(1 + exp(substitutionsForCurrentNode[t] - ll));
                    else
                      substitutionsForCurrentNode[t] += log(1 + exp(ll - substitutionsForCurrentNode[t]));
                  }
                } 
              }
          }
        }
      }
    }

    // Now we just have to copy the substitutions into the result vector:
    for (size_t t = 0; t < nbTypes; ++t)
    {
      double x = usesLog?exp(substitutionsForCurrentNode[t] - Lr):substitutionsForCurrentNode[t]/exp(Lr);
      
      if (std::isnan(x) || std::isinf(x))
      {
        if (verbose)
          ApplicationTools::displayWarning("On branch " + TextTools::toString(edid) + ", site " + TextTools::toString(site) + ", and type " + TextTools::toString(t) + ", counts could not be computed.");
        (*br)(0,t)=0;
      }
      else
      {
        if (threshold>=0 && x > threshold)
        {
          if (verbose)
            ApplicationTools::displayWarning("On branch " + TextTools::toString(edid) + ", site " + TextTools::toString(site) + ", and type " + TextTools::toString(t) + " count has been ignored because it is presumably saturated.");
          (*br)(0,t)=0;
        }
        else     
          (*br)(0,t)=x;;
      }
    }
  }
  
  if (verbose)
  {
    if (ApplicationTools::message)
      *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }
    
  return substitutions.release();
  */
}

/**************************************************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingToolsForASite::computeNormalizations(
  size_t site,
  LikelihoodCalculationSingleProcess& rltc,
  const vector<uint>& edgeIds,
  const BranchedModelSet* nullModels,
  const SubstitutionRegister& reg,
  std::shared_ptr<const AlphabetIndex2> distances,
  bool verbose)
{
  return 0;
  /*
  // Preamble:
  if (!rltc.isInitialized())
    throw Exception("SubstitutionMappingToolsForASite::computeNormalizations(). Likelihood object is not initialized.");

  const SubstitutionModel* sm(0);
  const SubstitutionProcess& sp=rltc.getSubstitutionProcess();

  if (edgeIds.size()==0)
  {
    const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();
    return new ProbabilisticSubstitutionMapping(ppt, reg.getNumberOfSubstitutionTypes(), 0);    
  }
  
  for (auto id :edgeIds)
  {
    if (dynamic_cast<const SubstitutionModel*>(sp.getModel(id,0))==NULL)
      throw Exception("SubstitutionMappingToolsForASite::computeNormalizations possible only for SubstitutionModels, not in branch " + TextTools::toString(id));
    else
      if (!sm)
        sm=dynamic_cast<const SubstitutionModel*>(sp.getModel(id,0));
  }
  
  size_t nbTypes = reg.getNumberOfSubstitutionTypes();
  size_t nbStates = sm->getAlphabet()->getSize();
  vector<int> supportedStates = sm->getAlphabetStates();
  
  const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();

  unique_ptr<ProbabilisticSubstitutionMapping> normalizations(new ProbabilisticSubstitutionMapping(ppt, nbTypes, 1));
  
  vector<size_t> vMod=nullModels->getModelNumbers();
  vector<UserAlphabetIndex1 >  usai(nbTypes, UserAlphabetIndex1(sm->getAlphabet()));

  for (auto& nbm : vMod)
  {
    const SubstitutionModel* modn = dynamic_cast<const SubstitutionModel*>(nullModels->getModel(nbm));
    if (!modn)
      throw Exception("SubstitutionMappingToolsForASite::computeNormalizations possible only for SubstitutionModels, not for model " + nullModels->getModel(nbm)->getName());
  }

  for (auto& nbm : vMod)
  {
    vector<uint> mids = VectorTools::vectorIntersection(edgeIds, nullModels->getBranchesWithModel(nbm));
    
    if (mids.size()>0)
    {
      const SubstitutionModel* modn = dynamic_cast<const SubstitutionModel*>(nullModels->getModel(nbm));

      for (size_t nbt = 0; nbt < nbTypes; nbt++)
        for (size_t i = 0; i < nbStates; i++)
          usai[nbt].setIndex(supportedStates[i], 0);

      for (size_t i = 0; i < nbStates; i++)
      {
        for (size_t j = 0; j < nbStates; j++)
        {
          if (i != j)
          {
            size_t nbt = reg.getType(i, j);
            if (nbt != 0)
              usai[nbt - 1].setIndex(supportedStates[i], usai[nbt - 1].getIndex(supportedStates[i]) + modn->Qij(i, j)*(distances?distances->getIndex(supportedStates[i],supportedStates[j]):1));
          }
        }
      }
      
      for (size_t nbt = 0; nbt < nbTypes; nbt++)
      {
        unique_ptr<Reward> reward(new DecompositionReward(modn, &usai[nbt]));
        
        unique_ptr<ProbabilisticRewardMapping> mapping(RewardMappingToolsForASite::computeRewardVectors(site, rltc, mids, *reward, verbose));

        for (size_t k = 0; k < mids.size(); k++)
        {
          shared_ptr<PhyloBranchMapping> brn=normalizations->getEdge(mids[k]);
          shared_ptr<PhyloBranchReward> brr=mapping->getEdge(mids[k]);

          (*brn)(0,nbt)=(*brr).getSiteReward(0);
        }
      }
    }
  }

  return normalizations.release();
  */
}

/************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingToolsForASite::computeNormalizedCounts(
  size_t site,
  LikelihoodCalculationSingleProcess& rltc,
  const vector<uint>& edgeIds,
  const BranchedModelSet* nullModels,
  const SubstitutionRegister& reg,
  std::shared_ptr<const AlphabetIndex2> weights,
  std::shared_ptr<const AlphabetIndex2> distances,
  bool perTimeUnit,
  uint siteSize,
  double threshold,
  bool verbose)
{
  unique_ptr<ProbabilisticSubstitutionMapping> counts(computeCounts(site,rltc,edgeIds,reg,weights,distances,threshold,verbose));
  
  unique_ptr<ProbabilisticSubstitutionMapping> factors(computeNormalizations(site,rltc,edgeIds,nullModels,reg,distances,verbose));
  
  return computeNormalizedCounts(counts.get(), factors.get(), edgeIds, perTimeUnit, siteSize);  
}
    
/************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingToolsForASite::computeNormalizedCounts(
  const ProbabilisticSubstitutionMapping* counts,
  const ProbabilisticSubstitutionMapping* factors,
  const vector<uint>& edgeIds,
  bool perTimeUnit,
  uint siteSize)
{
  unique_ptr<ProbabilisticSubstitutionMapping> normCounts(counts->clone());
  
  size_t nbTypes = counts->getNumberOfSubstitutionTypes();

  // Iterate on branches
  
  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt=normCounts->allEdgesIterator();

  for (;!brIt->end();brIt->next())
  {
    shared_ptr<PhyloBranchMapping> brNormCount=**brIt;

    Vdouble& ncou=brNormCount->getCounts()[0];

    // For each branch
    uint edid=normCounts->getEdgeIndex(brNormCount);

    if (edgeIds.size() > 0 && !VectorTools::contains(edgeIds, (int)edid))
    {
      VectorTools::fill(ncou,0.);
      continue;
    }
        
    shared_ptr<PhyloBranchMapping> brFactor=factors->getEdge(edid);
    shared_ptr<PhyloBranchMapping> brCount=counts->getEdge(edid);

    const Vdouble& cou=brCount->getCounts()[0];
    const Vdouble& fac=brFactor->getCounts()[0];


    // if not per time, multiply by the lengths of the branches of
    // the input tree
    
    double slg=(!perTimeUnit?brCount->getLength():1)/siteSize;
    
    for (size_t t = 0; t < nbTypes; ++t)
      ncou[t]=(fac[t]!=0? cou[t]/fac[t]*slg : 0);
  }

  return normCounts.release();
}

