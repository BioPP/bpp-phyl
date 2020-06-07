//
// File: SubstitutionMappingTools.cpp
// Created by: Julien Dutheil, Laurent Guéguen
// Created on: Wed Apr 5 13:04 2006
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

#include "SubstitutionMappingTools.h"
#include "DecompositionSubstitutionCount.h"
#include "DecompositionReward.h"
#include "ProbabilisticRewardMapping.h"
#include "ProbabilisticSubstitutionMapping.h"
#include "RewardMappingTools.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>

using namespace bpp;
using namespace std;

// From the STL:
#include <iomanip>

/******************************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeCounts(
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
    throw Exception("SubstitutionMappingTools::computeSubstitutionVectors(). Likelihood object is not initialized.");

  const SubstitutionModel* sm(0);
  const SubstitutionProcess& sp=rltc.getSubstitutionProcess();

  if (edgeIds.size()==0)
    return new ProbabilisticSubstitutionMapping(sp.getParametrizablePhyloTree(), reg.getNumberOfSubstitutionTypes(), rltc.getRootArrayPositions(), rltc.getNumberOfDistinctSites());
  
  for (auto id : edgeIds)
  {
    if (dynamic_cast<const SubstitutionModel*>(sp.getModelForNode(id))==NULL)
      throw Exception("SubstitutionMappingTools::computeSubstitutionVectors possible only for SubstitutionModels, not in branch " + TextTools::toString(id));
    else
      if (!sm)
        sm=dynamic_cast<const SubstitutionModel*>(sp.getModelForNode(id));
  }

  // We create a Mapping objects

  if (!sm)
    throw Exception("SubstitutionMappingTools::computeSubstitutionVectors not possible with null model.");

  unique_ptr<SubstitutionCount> substitutionCount(new DecompositionSubstitutionCount(sm, reg.clone(), weights, distances));

  return computeCounts(rltc, edgeIds, *substitutionCount, threshold, verbose);
}


ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeCounts(
  LikelihoodCalculationSingleProcess& rltc,
  const vector<uint>& edgeIds,
  SubstitutionCount& substitutionCount,
  double threshold,
  bool verbose)
{
  // Preamble:
  if (!rltc.isInitialized())
    throw Exception("SubstitutionMappingTools::computeCounts(). Likelihood object is not initialized.");

  const SubstitutionProcess& sp=rltc.getSubstitutionProcess();

  if (edgeIds.size()==0)
    return new ProbabilisticSubstitutionMapping(sp.getParametrizablePhyloTree(),
                                                substitutionCount.getNumberOfSubstitutionTypes(),
                                                rltc.getRootArrayPositions(),
                                                rltc.getNumberOfDistinctSites());
  
  auto processTree = rltc.getTreeNode(0);

  // Map from models to counts
  std::map<shared_ptr<ConfiguredModel> , std::shared_ptr<SubstitutionCount> > mModCount;
  
  for (auto speciesId :edgeIds)
  {
    const auto& dagIndexes = rltc.getEdgesIds(speciesId, 0);

    for (auto id:dagIndexes)
    {
      const auto& edge = processTree->getEdge(id);
      if (edge->getBrLen()) // if edge with model on it
      {
        auto model = edge->getModel();

        auto sm = dynamic_cast<const SubstitutionModel*>(model->getTargetValue());
        
        if (sm == NULL)
          throw Exception("SubstitutionMappingTools::computeCounts : SubstitutionVectors possible only for SubstitutionModels, not in branch " + TextTools::toString(speciesId));
        if (mModCount.find(model)==mModCount.end())
        {
          mModCount[model] = std::shared_ptr<SubstitutionCount>(substitutionCount.clone());
          mModCount[model]->setSubstitutionModel(sm);
        }
      }
    }
  }

  const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();
  
  // A few variables we'll need:

  size_t nbDistinctSites = rltc.getNumberOfDistinctSites();
  size_t nbClasses       = sp.getNumberOfClasses();

  size_t nbTypes         = substitutionCount.getNumberOfSubstitutionTypes();
  size_t nbNodes         = edgeIds.size();
  
  const auto& rootPatternLinks = rltc.getRootArrayPositions();

  // We create a Mapping objects
  
  unique_ptr<ProbabilisticSubstitutionMapping> substitutions(new ProbabilisticSubstitutionMapping(ppt, nbTypes, rootPatternLinks, nbDistinctSites));

  // Store likelihood for each compressed site :

  auto Lr = rltc.getSiteLikelihoods(true)->getTargetValue();

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

    vector<Eigen::RowVectorXd> substitutionsForCurrentNode(nbTypes);
     for (auto& sub:substitutionsForCurrentNode)
      sub.setZero(nbDistinctSites);

    for (size_t ncl=0; ncl<nbClasses; ncl++)
    {
      processTree = rltc.getTreeNode(ncl);

      double pr = sp.getProbabilityForModel(ncl);

      vector<Eigen::RowVectorXd> substitutionsForCurrentClass(nbTypes);
      for (auto& sub:substitutionsForCurrentClass)
        sub.setZero(nbDistinctSites);
      
      const auto& dagIndexes = rltc.getEdgesIds(speciesId, ncl);

      Eigen::MatrixXd npxy;

      // Sum on all dag edges for this speciesId
      for (auto id:dagIndexes)
      {
        auto edge = processTree->getEdge(id);

        auto subCount = mModCount[edge->getModel()];

        const auto& likelihoodsFather = rltc.getBackwardLikelihoodsAtEdgeForClass(id, ncl)->getTargetValue();

        auto sonid = processTree->getSon(id);
          
        const auto& likelihoodsSon = rltc.getForwardLikelihoodsAtNodeForClass(sonid, ncl)->getTargetValue();

        const Eigen::MatrixXd& pxy = edge->getTransitionMatrix()->accessValueConst();

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

    
    // Now we just have to copy the substitutions into the result vector:
    for (size_t i = 0; i < nbDistinctSites; ++i)
    {
      for (size_t t = 0; t < nbTypes; ++t)
      {
        double x = substitutionsForCurrentNode[t](i)/Lr(i);
        
        if (std::isnan(x) || std::isinf(x))
        {
          if (verbose)
            ApplicationTools::displayWarning("On branch " + TextTools::toString(speciesId) + ", site index " + TextTools::toString(i) + ", and type " + TextTools::toString(t) + ", counts could not be computed.");
          (*br)(i,t)=0;
        }
        else
        {
          if (threshold>=0 && x > threshold)
          {
            if (verbose)
              ApplicationTools::displayWarning("On branch " + TextTools::toString(speciesId) + ", site index" + TextTools::toString(i) + ", and type " + TextTools::toString(t) + " count has been ignored because it is presumably saturated.");
            (*br)(i,t)=0;
          }
          else     
            (*br)(i,t)= x;;
        }
      }
    }
  } // end of loop on branches
  
  if (verbose)
  {
    if (ApplicationTools::message)
      *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }

  return substitutions.release();
}

/**************************************************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeNormalizations(
  LikelihoodCalculationSingleProcess& rltc,
  const vector<uint>& edgeIds,
  const BranchedModelSet* nullModels,
  const SubstitutionRegister& reg,
  std::shared_ptr<const AlphabetIndex2> distances,
  bool verbose)
{
  // Preamble:
  if (!rltc.isInitialized())
    throw Exception("SubstitutionMappingTools::computeNormalizations(). Likelihood object is not initialized.");

  const SubstitutionProcess& sp=rltc.getSubstitutionProcess();
  const auto& statemap = sp.getStateMap();
  const auto* alphabet = statemap.getAlphabet();
  
  if (edgeIds.size()==0)
    return new ProbabilisticSubstitutionMapping(sp.getParametrizablePhyloTree(),
                                                reg.getNumberOfSubstitutionTypes(),
                                                rltc.getRootArrayPositions(),
                                                rltc.getNumberOfDistinctSites());

  // A few variables we'll need:

  size_t nbTypes = reg.getNumberOfSubstitutionTypes();

  size_t nbStates = statemap.getNumberOfModelStates();
  vector<int> supportedStates = statemap.getAlphabetStates();
  
  size_t nbDistinctSites = rltc.getNumberOfDistinctSites();
  
  unique_ptr<ProbabilisticSubstitutionMapping> normalizations(new ProbabilisticSubstitutionMapping(sp.getParametrizablePhyloTree(), nbTypes, rltc.getRootArrayPositions(), nbDistinctSites));

  
  vector<size_t> vMod=nullModels->getModelNumbers();

  for (auto& nbm : vMod)
  {
    const SubstitutionModel* modn = dynamic_cast<const SubstitutionModel*>(nullModels->getModel(nbm));
    if (!modn)
      throw Exception("SubstitutionMappingTools::computeNormalizations possible only for SubstitutionModels, not for model " + nullModels->getModel(nbm)->getName());
  }


  vector<UserAlphabetIndex1 >  usai(nbTypes, UserAlphabetIndex1(alphabet));

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
        
        unique_ptr<ProbabilisticRewardMapping> mapping(RewardMappingTools::computeRewardVectors(rltc, mids, *reward, verbose));

        for (size_t k = 0; k < mids.size(); k++)
        {
          shared_ptr<PhyloBranchMapping> brn=normalizations->getEdge(mids[k]);
          shared_ptr<PhyloBranchReward> brr=mapping->getEdge(mids[k]);

          for (size_t i = 0; i < nbDistinctSites; ++i)
            (*brn)(i,nbt)=(*brr).getSiteReward(i);
        }
      }
    }
  }

  return normalizations.release();
}

/************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeNormalizedCounts(
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
  unique_ptr<ProbabilisticSubstitutionMapping> counts(computeCounts(rltc,edgeIds,reg,weights,distances,threshold,verbose));
  
  unique_ptr<ProbabilisticSubstitutionMapping> factors(computeNormalizations(rltc,edgeIds,nullModels,reg,distances,verbose));
  
  return computeNormalizedCounts(counts.get(), factors.get(), edgeIds, perTimeUnit, siteSize);  
}
    
/************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeNormalizedCounts(
  const ProbabilisticSubstitutionMapping* counts,
  const ProbabilisticSubstitutionMapping* factors,
  const vector<uint>& edgeIds,
  bool perTimeUnit,
  uint siteSize)
{
  unique_ptr<ProbabilisticSubstitutionMapping> normCounts(counts->clone());
  
  size_t nbTypes = counts->getNumberOfSubstitutionTypes();
  size_t nbDistinctSites = counts->getNumberOfDistinctSites();

  // Iterate on branches
  
  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt=normCounts->allEdgesIterator();

  for (;!brIt->end();brIt->next())
  {
    shared_ptr<PhyloBranchMapping> brNormCount=**brIt;

    VVdouble& brnCou=brNormCount->getCounts();

    // For each branch
    uint edid=normCounts->getEdgeIndex(brNormCount);

    if (edgeIds.size() > 0 && !VectorTools::contains(edgeIds, (int)edid))
    {
      for (auto& brk : brnCou)
        VectorTools::fill(brk,0.);
      continue;
    }
        
    shared_ptr<PhyloBranchMapping> brFactor=factors->getEdge(edid);
    shared_ptr<PhyloBranchMapping> brCount=counts->getEdge(edid);

    const VVdouble& cou=brCount->getCounts();
    const VVdouble& fac=brFactor->getCounts();


    // if not per time, multiply by the lengths of the branches of
    // the input tree
    
    double slg=(!perTimeUnit?brCount->getLength():1)/siteSize;
    
    for (size_t k = 0; k < nbDistinctSites; k++)
    {
      Vdouble& ncou_k = brnCou[k];
      
      const Vdouble& cou_k=cou[k];
      const Vdouble& fac_k=fac[k];
      
      for (size_t t = 0; t < nbTypes; ++t)
        ncou_k[t]=(fac_k[t]!=0? cou_k[t]/fac_k[t]*slg : 0);
    }
  }

  return normCounts.release();
}

/*******************************************************************************/
/* Get trees of counts */
/*******************************************************************************/

PhyloTree* SubstitutionMappingTools::getTreeForType(const ProbabilisticSubstitutionMapping& counts,
                                                    size_t type)
{
  unique_ptr<PhyloTree> pt(new PhyloTree(counts));
  size_t nbSites = counts.getNumberOfSites();

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt=counts.allEdgesIterator();

  for (;!brIt->end();brIt->next())
  {
    shared_ptr<PhyloBranchMapping> brm=(**brIt);
    double x=0;
    
    for (size_t i = 0; i < nbSites; ++i)
      x += brm->getSiteTypeCount(counts.getSiteIndex(i), type);

    pt->getEdge(counts.getEdgeIndex(brm))->setLength(x);
  }
 
  return pt.release();
}

PhyloTree* SubstitutionMappingTools::getTreeForType(const ProbabilisticSubstitutionMapping& counts,
                                                    const ProbabilisticSubstitutionMapping& factors,
                                                    size_t type)
{
  unique_ptr<PhyloTree> pt(new PhyloTree(counts));
  size_t nbSites = counts.getNumberOfSites();

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt=counts.allEdgesIterator();

  for (;!brIt->end();brIt->next())
  {
    shared_ptr<PhyloBranchMapping> brm=(**brIt);
    shared_ptr<PhyloBranchMapping> brf=factors.getEdge(counts.getEdgeIndex(brm));
    
    double x=0,f=0;
    
    for (size_t i = 0; i < nbSites; ++i)
    {
      x += brm->getSiteTypeCount(counts.getSiteIndex(i), type);
      f += brf->getSiteTypeCount(counts.getSiteIndex(i), type);
    }
    
    pt->getEdge(counts.getEdgeIndex(brm))->setLength(x/f);
  }
 
  return pt.release();
}
  

/********************************************************************************/
/*  Get vectors of counts  */
/********************************************************************************/

VVVdouble SubstitutionMappingTools::getCountsPerSitePerBranchPerType(
  const ProbabilisticSubstitutionMapping& counts,
  const vector<uint>& ids)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);
  size_t nbBr = idc.size();
  
  size_t nbSites = counts.getNumberOfSites();

  size_t nbTypes= counts.getNumberOfSubstitutionTypes();

  VVVdouble result;
  VectorTools::resize3(result, nbSites, nbBr, nbTypes);

  for (size_t j = 0; j < nbSites; ++j)
  {
    VVdouble& resS=result[j];
    size_t siteIndex=counts.getSiteIndex(j);

    for (size_t k = 0; k < nbBr; ++k)
    {
      Vdouble& resSB=resS[k];

      for (size_t i = 0; i < nbTypes; ++i)
        resSB[i] = counts(idc[k], siteIndex, i);
    }
  }

  return result;
}



/**************************************************************************************************/

Vdouble SubstitutionMappingTools::getCountsForSitePerBranch(
  const ProbabilisticSubstitutionMapping& counts,
  size_t site)
{
  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt=counts.allEdgesIterator();
  size_t siteIndex = counts.getSiteIndex(site);
  
  Vdouble v(counts.getNumberOfBranches(),0);
  for (;!brIt->end();brIt->next())
    v[counts.getEdgeIndex(**brIt)] = VectorTools::sum((***brIt).getSiteCount(siteIndex));
  
  return v;
}

Vdouble SubstitutionMappingTools::getCountsForSitePerBranch(
  const ProbabilisticSubstitutionMapping& counts,
  const ProbabilisticSubstitutionMapping& factors,
  size_t site)
{
  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt=counts.allEdgesIterator();
  size_t siteIndex = counts.getSiteIndex(site);
  
  Vdouble v(counts.getNumberOfBranches(),0);
  for (;!brIt->end();brIt->next())
  {
    shared_ptr<PhyloBranchMapping> brm=(**brIt);
    
    uint edid=counts.getEdgeIndex(brm);
    shared_ptr<PhyloBranchMapping> brf=factors.getEdge(edid);

    
    v[edid] = VectorTools::sum(brm->getSiteCount(siteIndex))/VectorTools::sum(brf->getSiteCount(siteIndex));
  }
  
  return v;
}

/**************************************************************************************************/

VVdouble SubstitutionMappingTools::getCountsPerSitePerBranch(
  const ProbabilisticSubstitutionMapping& counts,
  const vector<uint>& ids)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);
  size_t nbBr = idc.size();

  size_t nbSites = counts.getNumberOfSites();
  
  VVdouble result;
  VectorTools::resize2(result,nbSites, nbBr);

  for (size_t k = 0; k < nbSites; ++k)
  {
    vector<double> countsf(SubstitutionMappingTools::getCountsForSitePerBranch(counts, k));
    Vdouble* resS=&result[k];
    
    for (size_t i = 0; i < nbBr; ++i)
      (*resS)[i]= countsf[idc[i]];
  }
  return result;
}

VVdouble SubstitutionMappingTools::getCountsPerSitePerBranch(
  const ProbabilisticSubstitutionMapping& counts,
  const ProbabilisticSubstitutionMapping& factors,
  const vector<uint>& ids)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);
  size_t nbBr = idc.size();

  size_t nbSites = counts.getNumberOfSites();
  
  VVdouble result;
  VectorTools::resize2(result,nbSites, nbBr);

  for (size_t k = 0; k < nbSites; ++k)
  {
    vector<double> countsf(SubstitutionMappingTools::getCountsForSitePerBranch(counts, factors, k));
    Vdouble* resS=&result[k];
    
    for (size_t i = 0; i < nbBr; ++i)
      (*resS)[i]= countsf[idc[i]];
  }
  return result;
}

/**************************************************************************************************/

Vdouble SubstitutionMappingTools::getCountsForBranchPerType(
  const ProbabilisticSubstitutionMapping& counts,
  uint branchId)
{
  size_t nbSites = counts.getNumberOfSites();
  size_t nbTypes = counts.getNumberOfSubstitutionTypes();
  Vdouble v(nbTypes, 0);
  shared_ptr<PhyloBranchMapping> br=counts.getEdge(branchId);

  for (size_t i = 0; i < nbSites; ++i)
    for (size_t t = 0; t < nbTypes; ++t)
      v[t] += br->getSiteTypeCount(counts.getSiteIndex(i), t);
 
  return v;
}

/**************************************************************************************************/

VVdouble SubstitutionMappingTools::getCountsPerBranchPerType(
  const ProbabilisticSubstitutionMapping& counts,
  const vector<uint>& ids)
{
  VVdouble result;
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);

  for (auto id : idc)
    result.push_back(SubstitutionMappingTools::getCountsForBranchPerType(counts, id));
    
  return result;
}

/**************************************************************************************************/

VVdouble SubstitutionMappingTools::getCountsPerTypePerBranch(
  const ProbabilisticSubstitutionMapping& counts,
  const vector<uint>& ids)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);
  size_t nbBr = idc.size();

  size_t nbTypes = counts.getNumberOfSubstitutionTypes();

  VVdouble result;
  VectorTools::resize2(result, nbTypes, nbBr);

  for (size_t i=0;i<idc.size();i++)
  {
    Vdouble cou=SubstitutionMappingTools::getCountsForBranchPerType(counts, idc[i]);
    for (size_t nbt=0;nbt<nbTypes;nbt++)
      result[nbt][i]=cou[nbt];
  }
  
  return result;
}


/************************************************************/


VVdouble SubstitutionMappingTools::computeCountsPerTypePerBranch(
  LikelihoodCalculationSingleProcess& rltc,
  const vector<uint>& ids,
  const SubstitutionRegister& reg,
  std::shared_ptr<const AlphabetIndex2> weights,
  std::shared_ptr<const AlphabetIndex2> distances,
  double threshold,
  bool verbose)
{
  ProbabilisticSubstitutionMapping psm(computeCounts(rltc,ids,reg,weights,distances,threshold,verbose));
  
  VVdouble result = getCountsPerTypePerBranch(psm, ids);

  const CategorySubstitutionRegister* creg = dynamic_cast<const CategorySubstitutionRegister*>(&reg);

  if ((creg!=NULL) && !creg->isStationary())
  {
    size_t nbTypes = result[0].size();

    for (size_t k = 0; k < ids.size(); ++k)
    {
      vector<double> freqs(0);// = rltc.getPosteriorStateFrequencies(ids[k]);
      // Compute frequencies for types:
      
      vector<double> freqsTypes(creg->getNumberOfCategories());
      for (size_t i = 0; i < freqs.size(); ++i)
      {
        size_t c = creg->getCategory(i);
        freqsTypes[c - 1] += freqs[i];
      }

      // We devide the counts by the frequencies and rescale:
      
      double s = VectorTools::sum(result[k]);
      for (size_t t = 0; t < nbTypes; ++t)
      {
        result[k][t] /= freqsTypes[creg->getCategoryFrom(t + 1) - 1];
      }

      double s2 = VectorTools::sum(result[k]);
      // Scale:
      
      result[k] = (result[k] / s2) * s;
    }
  }
  return result;
}  


/**************************************************************************************************/

Vdouble SubstitutionMappingTools::getCountsForSitePerType(const ProbabilisticSubstitutionMapping& counts, size_t site, const vector<uint>& ids)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);

  size_t siteIndex = counts.getSiteIndex(site);

  size_t nbTypes    = counts.getNumberOfSubstitutionTypes();
  Vdouble v(nbTypes,0);

  for (auto id : idc)
  {
    shared_ptr<PhyloBranchMapping> br=counts.getEdge(id);
    for (size_t t = 0; t < nbTypes; ++t)
      v[t] += (*br)(siteIndex, t);
  }
  
  return v;
}

/**************************************************************************************************/

Vdouble SubstitutionMappingTools::getCountsForSitePerType(
  const ProbabilisticSubstitutionMapping& counts,
  const ProbabilisticSubstitutionMapping& factors,
  size_t site,
  bool perTimeUnit,
  uint siteSize)
 {
  size_t siteIndex = counts.getSiteIndex(site);
  size_t nbTypes   = counts.getNumberOfSubstitutionTypes();
  
  Vdouble v(nbTypes,0);
  Vdouble n(nbTypes,0);

  double lg=(!perTimeUnit?0:1);
  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt=counts.allEdgesIterator();
  for (;!brIt->end();brIt->next())
  {
    for (size_t t = 0; t < nbTypes; ++t)
      v[t] += (***brIt)(siteIndex, t);
    if (!perTimeUnit)
      lg+=(**brIt)->getLength();
  }

  double slg=lg/siteSize;

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> nrIt=factors.allEdgesIterator();
  for (;!nrIt->end();nrIt->next())
    for (size_t t = 0; t < nbTypes; ++t)
      n[t] += (***nrIt)(siteIndex, t);

  for (size_t t = 0; t < nbTypes; ++t)
    v[t] = v[t]/n[t]*slg;

  return v;
}

/**************************************************************************************************/

Vdouble SubstitutionMappingTools::getCountsForSitePerType(
  const ProbabilisticSubstitutionMapping& counts,
  const ProbabilisticSubstitutionMapping& factors,
  size_t site,
  const vector<uint>& ids,
  bool perTimeUnit,
  uint siteSize)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);

  size_t siteIndex = counts.getSiteIndex(site);
  size_t nbTypes   = counts.getNumberOfSubstitutionTypes();
  
  Vdouble v(nbTypes,0);
  Vdouble n(nbTypes,0);

  double lg=(!perTimeUnit?0:1);
  for (auto id : idc)
  {
    shared_ptr<PhyloBranchMapping> br=counts.getEdge(id);
    for (size_t t = 0; t < nbTypes; ++t)
      v[t] += (*br)(siteIndex, t);
    if (!perTimeUnit)
      lg+=br->getLength();
  }

  double slg=lg/siteSize;
  
  for (auto id : idc)
  {
    shared_ptr<PhyloBranchMapping> br=factors.getEdge(id);
    for (size_t t = 0; t < nbTypes; ++t)
      n[t] += (*br)(siteIndex, t);
  }
  
  for (size_t t = 0; t < nbTypes; ++t)
    v[t] = v[t]/n[t]*slg;

  return v;
}

/**************************************************************************************************/

VVdouble SubstitutionMappingTools::getCountsPerSitePerType(
  const ProbabilisticSubstitutionMapping& counts,
  const ProbabilisticSubstitutionMapping& factors,
  bool perTimeUnit,
  uint siteSize)
{
  size_t nbSites = counts.getNumberOfSites();
  size_t nbTypes = counts.getNumberOfSubstitutionTypes();
  VVdouble result;
  VectorTools::resize2(result, nbSites, nbTypes);

  for (size_t k = 0; k < nbSites; ++k)
    result[k]= SubstitutionMappingTools::getCountsForSitePerType(counts, factors, k, perTimeUnit, siteSize);

  return result;
}

/**************************************************************************************************/

VVdouble SubstitutionMappingTools::getCountsPerSitePerType(const ProbabilisticSubstitutionMapping& counts, const vector<uint>& ids)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);
  
  size_t nbSites = counts.getNumberOfSites();
  size_t nbTypes = counts.getNumberOfSubstitutionTypes();
  VVdouble result;
  VectorTools::resize2(result, nbSites, nbTypes);

  for (size_t k = 0; k < nbSites; ++k)
    result[k]= SubstitutionMappingTools::getCountsForSitePerType(counts,k,idc);

  return result;
}

/**************************************************************************************************/

VVdouble SubstitutionMappingTools::getCountsPerSitePerType(
  const ProbabilisticSubstitutionMapping& counts,
  const ProbabilisticSubstitutionMapping& factors,
  const vector<uint>& ids,
  bool perTimeUnit,
  uint siteSize)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);

  size_t nbSites = counts.getNumberOfSites();
  size_t nbTypes = counts.getNumberOfSubstitutionTypes();
  VVdouble result;
  VectorTools::resize2(result, nbSites, nbTypes);

  for (size_t k = 0; k < nbSites; ++k)
    result[k]= SubstitutionMappingTools::getCountsForSitePerType(counts, factors, k, idc, perTimeUnit, siteSize);

  return result;
}

/**************************************************************************************************/

double SubstitutionMappingTools::getNormForSite(const ProbabilisticSubstitutionMapping& counts, size_t site)
{
  size_t siteIndex = counts.getSiteIndex(site);
  double sumSquare = 0;

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt=counts.allEdgesIterator();

  size_t nbTypes    = counts.getNumberOfSubstitutionTypes();
  for (;!brIt->end();brIt->next())
  {
    double sum = 0;
    for (size_t t = 0; t < nbTypes; ++t)
      sum += (***brIt)(siteIndex, t);

    sumSquare += sum * sum;
  }
  return sqrt(sumSquare);
}





/******************************************************/
                 /* OUTPUT */
/******************************************************/

void SubstitutionMappingTools::outputPerSitePerBranch(
  const string& filename,
  const vector<uint>& ids,
  const VVdouble& counts)
{
  size_t nbSites=counts.size();
  if (nbSites==0)
    return;
  size_t nbBr=counts[0].size();
  
  ofstream file;
  file.open(filename.c_str());

  file << "sites";
  for (size_t i = 0; i < nbBr; ++i)
  {
    file << "\t" << ids[i];
  }
  file << endl;

  for (size_t k = 0; k < nbSites; ++k)
  {
    const Vdouble& countS=counts[k];
    file << k;
    for (size_t i = 0; i < nbBr; ++i)
    {
      file << "\t" << countS[i];
    }
    file << endl;
  }
  file.close();
}


/**************************************************************************************************/

void SubstitutionMappingTools::outputPerSitePerType(
  const string& filename,
  const SubstitutionRegister& reg,
  const VVdouble& counts)
{
  
  size_t nbSites=counts.size();
  if (nbSites==0)
    return;
  size_t nbTypes=counts[0].size();
  
  ofstream file;
  file.open(filename.c_str());
  
  file << "sites";
  for (size_t i = 0; i < nbTypes; ++i)
  {
    file << "\t" << reg.getTypeName(i + 1);
  }
  file << endl;

  for (size_t k = 0; k < nbSites; ++k)
  {
    file << k;
    const Vdouble& resS=counts[k];
    for (size_t i = 0; i < nbTypes; ++i)
    {
      file << "\t" << resS[i];
    }
    file << endl;
  }
  file.close();
}


/**************************************************************************************************/

void SubstitutionMappingTools::outputPerSitePerBranchPerType(
  const string& filenamePrefix,
  const vector<uint>& ids,
  const SubstitutionRegister& reg,
  const VVVdouble& counts)
{
  size_t nbSites=counts.size();
  if (nbSites==0)
    return;
  size_t nbBr = counts[0].size();
  if (nbBr==0)
    return;
  size_t nbTypes = counts[0][0].size();

  ofstream file;

  for (size_t i = 0; i < nbTypes; ++i)
  {
    string name=reg.getTypeName(i+1);
    if (name=="")
      name=TextTools::toString(i + 1);

    string path = filenamePrefix + name + string(".count");
    
    ApplicationTools::displayResult(string("Output counts of type ") + TextTools::toString(i + 1) + string(" to file"), path);
    file.open(path.c_str());
  
    file << "sites";
    for (size_t k = 0; k < nbBr; ++k)
    {
      file << "\t" << ids[k];
    }
    file << endl;

    for (size_t j = 0; j < nbSites; ++j)
    {
      const VVdouble& resS=counts[j];
      
      file << j;
      for (size_t k = 0; k < nbBr; ++k)
      {
        file << "\t" << resS[k][i];
      }
      file << endl;
    }
    file.close();
  }
}


/**************************************************************************************************/

void SubstitutionMappingTools::writeToStream(
  const ProbabilisticSubstitutionMapping& substitutions,
  const AlignedValuesContainer& sites,
  size_t type,
  ostream& out)
{
  if (!out)
    throw IOException("SubstitutionMappingTools::writeToFile. Can't write to stream.");
  out << "Branches";
  out << "\tMean";
  for (size_t i = 0; i < substitutions.getNumberOfSites(); i++)
  {
    out << "\tSite" << sites.getSymbolListSite(i).getPosition();
  }
  out << endl;

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt=substitutions.allEdgesIterator();

  for (;!brIt->end();brIt->next())
  {
    const shared_ptr<PhyloBranchMapping> br=**brIt;
    
    out << substitutions.getEdgeIndex(br) << "\t" << br->getLength();

    for (size_t i = 0; i < substitutions.getNumberOfSites(); i++)
      out << "\t" << (*br)(substitutions.getSiteIndex(i), type);

    out << endl;
  }
}

/**************************************************************************************************/

void SubstitutionMappingTools::readFromStream(istream& in, ProbabilisticSubstitutionMapping& substitutions, size_t type)
{
  try
  {
    DataTable* data = DataTable::read(in, "\t", true, -1);
    vector<string> ids = data->getColumn(0);
    data->deleteColumn(0); // Remove ids
    data->deleteColumn(0); // Remove means
    // Now parse the table:
    size_t nbSites = data->getNumberOfColumns();
    substitutions.setNumberOfSites(nbSites);

    unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt=substitutions.allEdgesIterator();

    for (;!brIt->end();brIt->next())
    {
      const shared_ptr<PhyloBranchMapping> br=**brIt;
      uint brid = substitutions.getEdgeIndex(br);
      for (size_t j = 0; j < nbSites; j++)
        (*br)(j,type) = TextTools::toDouble((*data)(brid, j));
    }
    
    // Parse the header:
    for (size_t i = 0; i < nbSites; i++)
    {
      string siteTxt = data->getColumnName(i);
      int site = 0;
      if (siteTxt.substr(0, 4) == "Site")
        site = TextTools::to<int>(siteTxt.substr(4));
      else
        site = TextTools::to<int>(siteTxt);
      substitutions.setSitePosition(i, site);
    }

    delete data;
  }
  catch (Exception& e)
  {
    throw IOException(string("Bad input file. ") + e.what());
  }
}


