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
  RecursiveLikelihoodTreeCalculation& rltc,
  const vector<uint>& nodeIds,
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
  const SubstitutionProcess& sp=*rltc.getSubstitutionProcess();

  if (nodeIds.size()==0)
  {
    const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();
    const vector<size_t>& rootPatternLinks = rltc.getLikelihoodData().getRootArrayPositions();
    size_t nbDistinctSites = rltc.getLikelihoodData().getNumberOfDistinctSites();
    
    return new ProbabilisticSubstitutionMapping(ppt, reg.getNumberOfSubstitutionTypes(), rootPatternLinks, nbDistinctSites);    
  }
  
  for (auto id : nodeIds)
  {
    if (dynamic_cast<const SubstitutionModel*>(sp.getModel(id,0))==NULL)
      throw Exception("SubstitutionMappingTools::computeSubstitutionVectors possible only for SubstitutionModels, not in branch " + TextTools::toString(id));
    else
      if (!sm)
        sm=dynamic_cast<const SubstitutionModel*>(sp.getModel(id,0));
  }

  // We create a Mapping objects

  if (!sm)
    throw Exception("SubstitutionMappingTools::computeSubstitutionVectors not possible with null model.");

  unique_ptr<SubstitutionCount> substitutionCount(new DecompositionSubstitutionCount(sm, reg.clone(), weights, distances));

  return computeCounts(rltc, nodeIds, *substitutionCount, threshold, verbose);
}


ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeCounts(
  RecursiveLikelihoodTreeCalculation& rltc,
  const vector<uint>& nodeIds,
  SubstitutionCount& substitutionCount,
  double threshold,
  bool verbose)
{
  // Preamble:
  if (!rltc.isInitialized())
    throw Exception("SubstitutionMappingTools::computeSubstitutionVectors(). Likelihood object is not initialized.");
  rltc.computeTreeLikelihood();
  
  const SubstitutionProcess& sp=*rltc.getSubstitutionProcess();

  if (nodeIds.size()==0)
  {
    const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();
    const vector<size_t>& rootPatternLinks = rltc.getLikelihoodData().getRootArrayPositions();
    size_t nbDistinctSites = rltc.getLikelihoodData().getNumberOfDistinctSites();
    
    return new ProbabilisticSubstitutionMapping(ppt, substitutionCount.getNumberOfSubstitutionTypes(), rootPatternLinks, nbDistinctSites);
  }

  for (auto id :nodeIds)
  {
    if (dynamic_cast<const SubstitutionModel*>(sp.getModel(id,0))==NULL)
      throw Exception("SubstitutionMappingTools::computeSubstitutionVectors possible only for SubstitutionModels, not in branch " + TextTools::toString(id));
  }

  const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();
  
  // A few variables we'll need:

  const RecursiveLikelihoodTree& rlt=dynamic_cast<const RecursiveLikelihoodTree&>(rltc.getLikelihoodData());

  size_t nbDistinctSites = rlt.getNumberOfDistinctSites();
  size_t nbStates        = sp.getNumberOfStates();
  size_t nbClasses       = sp.getNumberOfClasses();

  size_t nbTypes         = substitutionCount.getNumberOfSubstitutionTypes();
  size_t nbNodes         = nodeIds.size();
  
  const vector<size_t>& rootPatternLinks = rlt.getRootArrayPositions();

  // We create a Mapping objects
  
  unique_ptr<ProbabilisticSubstitutionMapping> substitutions(new ProbabilisticSubstitutionMapping(ppt, nbTypes, rootPatternLinks, nbDistinctSites));

  // Store likelihood for each compressed site :

  Vdouble Lr(nbDistinctSites, 0);

  for (size_t i = 0; i < nbDistinctSites; i++)
    Lr[i]=rltc.getLogLikelihoodForASiteIndex(i);

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
    uint edid=substitutions->getEdgeIndex(br);

    if (nodeIds.size() > 0 && !VectorTools::contains(nodeIds, (int)edid))
      continue;

    uint fathid=substitutions->getFather(edid);
    uint icid=substitutions->getSon(edid);

    double d=br->getLength();
    
    VVdouble likelihoodsFatherConstantPart;    
    VectorTools::resize2(likelihoodsFatherConstantPart, nbDistinctSites, nbStates);

    // now the counts
      
    VVdouble substitutionsForCurrentNode;
    VectorTools::resize2(substitutionsForCurrentNode,nbDistinctSites,nbTypes);

    bool usesLog=false;
    
    for (size_t ncl=0; ncl<nbClasses; ncl++)
    {
      const RecursiveLikelihoodTree::LikTree& rlt_c=rlt[ncl];
      
      shared_ptr<RecursiveLikelihoodNode> ici = rlt_c.getNode(icid);

      // reinit substitutionsForCurrentNode for log 
      if (!usesLog && ici->usesLog())
      {
        for (auto& vi : substitutionsForCurrentNode)
          std::fill(vi.begin(), vi.end(), NumConstants::MINF());
      }
        
      usesLog=ici->usesLog();

      double pr=sp.getProbabilityForModel(ncl);
      double rate=sp.getRateForModel(ncl);
      
      for (size_t i=0; i<nbDistinctSites; i++)
        VectorTools::fill(likelihoodsFatherConstantPart[i],usesLog?log(pr):pr);

      rltc.computeLikelihoodsAtNode(fathid);

      // up father
      shared_ptr<RecursiveLikelihoodNode> father = rlt_c.getNode(fathid);
      
      // down the brothers

      unique_ptr<RecursiveLikelihoodTree::LikTree::EdgeIterator> brothIt=rlt_c.branchesIterator(father);

      for (;!brothIt->end();brothIt->next())
      {
        if (rlt_c.getEdgeIndex(**brothIt)!=edid)
        {
          bool slog=rlt_c.getSon(**brothIt)->usesLog();
          
          if (!usesLog)
          {
            if (!slog)
              likelihoodsFatherConstantPart *= rlt_c.getSon(**brothIt)->getToFatherBelowLikelihoodArray(ComputingNode::D0);
            else
              likelihoodsFatherConstantPart *= VectorTools::exp(rlt_c.getSon(**brothIt)->getToFatherBelowLikelihoodArray(ComputingNode::D0));
          }
          else {
            if (slog)
              likelihoodsFatherConstantPart += rlt_c.getSon(**brothIt)->getToFatherBelowLikelihoodArray(ComputingNode::D0);
            else
              likelihoodsFatherConstantPart += VectorTools::log(rlt_c.getSon(**brothIt)->getToFatherBelowLikelihoodArray(ComputingNode::D0));
          }
        }
      }

      bool flog=father->usesLog();

      if (!usesLog)
      {
        if (!flog)
          likelihoodsFatherConstantPart *= father->getAboveLikelihoodArray();
        else
          likelihoodsFatherConstantPart *= VectorTools::exp(father->getAboveLikelihoodArray());
      }
      else {
        if (flog)
          likelihoodsFatherConstantPart += father->getAboveLikelihoodArray();
        else
          likelihoodsFatherConstantPart += VectorTools::log(father->getAboveLikelihoodArray());
      }
      
      // Then, we deal with the node of interest.
      // We first average upon 'y' to save computations, and then upon 'x'.
      // ('y' is the state at 'node' and 'x' the state at 'father'.)

      const VVdouble& likelihoodsFather_node = ici->getBelowLikelihoodArray(ComputingNode::D0);

      const SubstitutionModel* sm=dynamic_cast<const SubstitutionModel*>(sp.getModel(icid,ncl));
      if (!sm)
        throw Exception("SubstitutionMappingTools:: non substitution model in node " + TextTools::toString(icid));
      
      substitutionCount.setSubstitutionModel(sm);
      
      // compute all nxy * pxy first:
      VVVdouble npxy;

      VectorTools::resize3(npxy,nbTypes,nbStates,nbStates);

      const Matrix<double>& pxy = sp.getTransitionProbabilities(icid,ncl);
      
      for (size_t t = 0; t < nbTypes; ++t)
      {
        VVdouble& npxy_t=npxy[t];
        Matrix<double>* nijt = substitutionCount.getAllNumbersOfSubstitutions(d * rate, t + 1);
        MatrixTools::hadamardMult((*nijt),pxy,(*nijt));

        for (size_t x=0; x<nbStates; x++)
          for (size_t y=0; y<nbStates; y++)
            npxy_t[x][y]=(*nijt)(x,y);
        
        delete nijt;
      }
      
      // Now loop over sites:

      for (size_t i=0; i< nbDistinctSites; i++)
      {
        const Vdouble* likelihoodsFather_node_i = &(likelihoodsFather_node[i]);
        const Vdouble* likelihoodsFatherConstantPart_i = &(likelihoodsFatherConstantPart[i]);
        
        for (size_t x = 0; x < nbStates; ++x)
        {
          double likelihoodsFatherConstantPart_i_x = (*likelihoodsFatherConstantPart_i)[x];
          for (size_t y = 0; y < nbStates; ++y)
          {
            double likelihood_xy = usesLog
              ?likelihoodsFatherConstantPart_i_x + (*likelihoodsFather_node_i)[y]
              :likelihoodsFatherConstantPart_i_x * (*likelihoodsFather_node_i)[y];
            
            if (!usesLog)
            {
              if (likelihood_xy!=0) // to avoid multiplication per nan
                                    // (stop codons)
                for (size_t t = 0; t < nbTypes; ++t)
                {
                  substitutionsForCurrentNode[i][t] += likelihood_xy * npxy[t][x][y];
                 //                                <------------>  <----------->
                 // Posterior probability              |               |
                 // for site i *                       |               |
                 // likelihood for this site ----------+               |
                 //                                                    |
                 // Substitution function for site i ------------------+
                }
            }
            else
            {
              if (likelihood_xy!=NumConstants::MINF())  // to avoid add per -inf
                for (size_t t = 0; t < nbTypes; ++t)
                {
                  if (npxy[t][x][y]< -NumConstants::MILLI())
                  {
                    ApplicationTools::displayWarning("These counts are negative, their logs could not be computed:" + TextTools::toString(npxy[t][x][y]));
                    throw Exception("Stop in SubstitutionMappingTools");
                  }
                  else
                  {
                    if (npxy[t][x][y]>0)
                    {
                      double ll=likelihood_xy + log(npxy[t][x][y]);
                      if (ll>substitutionsForCurrentNode[i][t])
                        substitutionsForCurrentNode[i][t] = ll + log(1 + exp(substitutionsForCurrentNode[i][t] - ll));
                      else
                        substitutionsForCurrentNode[i][t] += log(1 + exp(ll - substitutionsForCurrentNode[i][t]));
                    }
                  } 
                }
            }
          }
        }
      }
    }

    // Now we just have to copy the substitutions into the result vector:
    for (size_t i = 0; i < nbDistinctSites; ++i)
    {
      for (size_t t = 0; t < nbTypes; ++t)
      {
        double x = usesLog?exp(substitutionsForCurrentNode[i][t] - Lr[i]):substitutionsForCurrentNode[i][t]/exp(Lr[i]);
        
        if (std::isnan(x) || std::isinf(x))
        {
          if (verbose)
            ApplicationTools::displayWarning("On branch " + TextTools::toString(edid) + ", site index " + TextTools::toString(i) + ", and type " + TextTools::toString(t) + ", counts could not be computed.");
          (*br)(i,t)=0;
        }
        else
        {
          if (threshold>=0 && x > threshold)
          {
            if (verbose)
              ApplicationTools::displayWarning("On branch " + TextTools::toString(edid) + ", site index" + TextTools::toString(i) + ", and type " + TextTools::toString(t) + " count has been ignored because it is presumably saturated.");
            (*br)(i,t)=0;
          }
          else     
            (*br)(i,t)= x;;
        }
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
}

/**************************************************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeNormalizations(
  RecursiveLikelihoodTreeCalculation& rltc,
  const vector<uint>& nodeIds,
  const BranchedModelSet* nullModels,
  const SubstitutionRegister& reg,
  std::shared_ptr<const AlphabetIndex2> distances,
  bool verbose)
{
  // Preamble:
  if (!rltc.isInitialized())
    throw Exception("SubstitutionMappingTools::computeNormalizations(). Likelihood object is not initialized.");
  rltc.computeTreeLikelihood();

  const SubstitutionModel* sm(0);
  const SubstitutionProcess& sp=*rltc.getSubstitutionProcess();

  if (nodeIds.size()==0)
  {
    const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();
    const vector<size_t>& rootPatternLinks = rltc.getLikelihoodData().getRootArrayPositions();
    size_t nbDistinctSites = rltc.getLikelihoodData().getNumberOfDistinctSites();
    
    return new ProbabilisticSubstitutionMapping(ppt, reg.getNumberOfSubstitutionTypes(), rootPatternLinks, nbDistinctSites);    
  }
  
  for (auto id :nodeIds)
  {
    if (dynamic_cast<const SubstitutionModel*>(sp.getModel(id,0))==NULL)
      throw Exception("SubstitutionMappingTools::computeNormalizations possible only for SubstitutionModels, not in branch " + TextTools::toString(id));
    else
      if (!sm)
        sm=dynamic_cast<const SubstitutionModel*>(sp.getModel(id,0));
  }

  const RecursiveLikelihoodTree& rlt=dynamic_cast<const RecursiveLikelihoodTree&>(rltc.getLikelihoodData());
  const vector<size_t>& rootPatternLinks = rlt.getRootArrayPositions();

  size_t nbTypes = reg.getNumberOfSubstitutionTypes();
  size_t nbStates = sm->getAlphabet()->getSize();
  vector<int> supportedStates = sm->getAlphabetStates();
  size_t nbDistinctSites = rlt.getNumberOfDistinctSites();
  
  const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();

  unique_ptr<ProbabilisticSubstitutionMapping> normalizations(new ProbabilisticSubstitutionMapping(ppt, nbTypes, rootPatternLinks, nbDistinctSites));
  
  vector<size_t> vMod=nullModels->getModelNumbers();
  vector<UserAlphabetIndex1 >  usai(nbTypes, UserAlphabetIndex1(sm->getAlphabet()));

  for (auto& nbm : vMod)
  {
    const SubstitutionModel* modn = dynamic_cast<const SubstitutionModel*>(nullModels->getModel(nbm));
    if (!modn)
      throw Exception("SubstitutionMappingTools::computeNormalizations possible only for SubstitutionModels, not for model " + nullModels->getModel(nbm)->getName());
  }

  for (auto& nbm : vMod)
  {
    vector<uint> mids = VectorTools::vectorIntersection(nodeIds, nullModels->getBranchesWithModel(nbm));
    
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
  RecursiveLikelihoodTreeCalculation& rltc,
  const vector<uint>& nodeIds,
  const BranchedModelSet* nullModels,
  const SubstitutionRegister& reg,
  std::shared_ptr<const AlphabetIndex2> weights,
  std::shared_ptr<const AlphabetIndex2> distances,
  bool perTimeUnit,
  uint siteSize,
  double threshold,
  bool verbose)
{
  unique_ptr<ProbabilisticSubstitutionMapping> counts(computeCounts(rltc,nodeIds,reg,weights,distances,threshold,verbose));
  
  unique_ptr<ProbabilisticSubstitutionMapping> factors(computeNormalizations(rltc,nodeIds,nullModels,reg,distances,verbose));
  
  return computeNormalizedCounts(counts.get(), factors.get(), nodeIds, perTimeUnit, siteSize);  
}
    
/************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeNormalizedCounts(
  const ProbabilisticSubstitutionMapping* counts,
  const ProbabilisticSubstitutionMapping* factors,
  const vector<uint>& nodeIds,
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

    if (nodeIds.size() > 0 && !VectorTools::contains(nodeIds, (int)edid))
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
  RecursiveLikelihoodTreeCalculation& rltc,
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
      vector<double> freqs = rltc.getLikelihoodData().getPosteriorStateFrequencies(ids[k]);
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


/**************************************************************************************************/

/*ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(
  const DRTreeLikelihood& drtl,
  SubstitutionCount& substitutionCount,
  bool verbose) throw (Exception)
  {
  // Preamble:
  if (!drtl.isInitialized())
  throw Exception("SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(). Likelihood object is not initialized.");

  // A few variables we'll need:
  const TreeTemplate<Node> tree(drtl.getTree());
  const AlignedValuesContainer*    sequences = drtl.getData();
  const DiscreteDistribution* rDist = drtl.getRateDistribution();

  size_t nbSites         = sequences->getNumberOfSites();
  size_t nbDistinctSites = drtl.getLikelihoodData()->getNumberOfDistinctSites();
  size_t nbStates        = sequences->getAlphabet()->getSize();
  size_t nbClasses       = rDist->getNumberOfCategories();
  size_t nbTypes         = substitutionCounts.getNumberOfSubstitutionTypes();
  vector<const Node*> nodes   = tree.getNodes();
  const vector<size_t>* rootPatternLinks
  = &drtl.getLikelihoodData()->getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  size_t nbNodes = nodes.size();

  // We create a new ProbabilisticSubstitutionMapping object:
  ProbabilisticSubstitutionMapping* substitutions = new ProbabilisticSubstitutionMapping(tree, &substitutionCount, nbSites);

Vdouble rcRates = rDist->getCategories();

// Compute the number of substitutions for each class and each branch in the tree:
if (verbose)
ApplicationTools::displayTask("Compute joint node-pairs likelihood", true);

for (size_t l = 0; l < nbNodes; ++l)
{
// For each node,
const Node* currentNode = nodes[l];

const Node* father = currentNode->getFather();

double d = currentNode->getDistanceToFather();

if (verbose)
ApplicationTools::displayGauge(l, nbNodes - 1);
VVdouble substitutionsForCurrentNode(nbDistinctSites);
for (size_t i = 0; i < nbDistinctSites; ++i)
{
substitutionsForCurrentNode[i].resize(nbTypes);
}

// Now we've got to compute likelihoods in a smart manner... ;)
VVVdouble likelihoodsFatherConstantPart(nbDistinctSites);
for (size_t i = 0; i < nbDistinctSites; ++i)
{
VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
likelihoodsFatherConstantPart_i->resize(nbClasses);
for (size_t c = 0; c < nbClasses; ++c)
{
Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
likelihoodsFatherConstantPart_i_c->resize(nbStates);
double rc = rDist->getProbability(c);
for (size_t s = 0; s < nbStates; ++s)
{
// (* likelihoodsFatherConstantPart_i_c)[s] = rc * model->freq(s);
// freq is already accounted in the array
(*likelihoodsFatherConstantPart_i_c)[s] = rc;
}
}
}

// First, what will remain constant:
size_t nbSons =  father->getNumberOfSons();
for (size_t n = 0; n < nbSons; ++n)
{
const Node* currentSon = father->getSon(n);
if (currentSon->getId() != currentNode->getId())
{
const VVVdouble* likelihoodsFather_son = &drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentSon->getId());

// Now iterate over all site partitions:
unique_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(currentSon->getId()));
VVVdouble pxy;
bool first;
while (mit->hasNext())
{
TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
unique_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
first = true;
while (sit->hasNext())
{
size_t i = sit->next();
// We retrieve the transition probabilities for this site partition:
if (first)
{
pxy = drtl.getTransitionProbabilitiesPerRateClass(currentSon->getId(), i);
first = false;
}
const VVdouble* likelihoodsFather_son_i = &(*likelihoodsFather_son)[i];
VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
for (size_t c = 0; c < nbClasses; ++c)
{
const Vdouble* likelihoodsFather_son_i_c = &(*likelihoodsFather_son_i)[c];
Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
VVdouble* pxy_c = &pxy[c];
for (size_t x = 0; x < nbStates; ++x)
{
Vdouble* pxy_c_x = &(*pxy_c)[x];
double likelihood = 0.;
for (size_t y = 0; y < nbStates; ++y)
{
likelihood += (*pxy_c_x)[y] * (*likelihoodsFather_son_i_c)[y];
}
(*likelihoodsFatherConstantPart_i_c)[x] *= likelihood;
}
}
}
}
}
}
if (father->hasFather())
{
const Node* currentSon = father->getFather();
const VVVdouble* likelihoodsFather_son = &drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentSon->getId());
// Now iterate over all site partitions:
unique_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(father->getId()));
VVVdouble pxy;
bool first;
while (mit->hasNext())
{
TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
unique_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
first = true;
while (sit->hasNext())
{
size_t i = sit->next();
// We retrieve the transition probabilities for this site partition:
if (first)
{
pxy = drtl.getTransitionProbabilitiesPerRateClass(father->getId(), i);
first = false;
}
const VVdouble* likelihoodsFather_son_i = &(*likelihoodsFather_son)[i];
VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
for (size_t c = 0; c < nbClasses; ++c)
{
const Vdouble* likelihoodsFather_son_i_c = &(*likelihoodsFather_son_i)[c];
Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
VVdouble* pxy_c = &pxy[c];
for (size_t x = 0; x < nbStates; ++x)
{
double likelihood = 0.;
for (size_t y = 0; y < nbStates; ++y)
{
Vdouble* pxy_c_x = &(*pxy_c)[y];
likelihood += (*pxy_c_x)[x] * (*likelihoodsFather_son_i_c)[y];
}
(*likelihoodsFatherConstantPart_i_c)[x] *= likelihood;
}
}
}
}
}
else
{
// Account for root frequencies:
for (size_t i = 0; i < nbDistinctSites; ++i)
{
vector<double> freqs = drtl.getRootFrequencies(i);
VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
for (size_t c = 0; c < nbClasses; ++c)
{
Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
for (size_t x = 0; x < nbStates; ++x)
{
(*likelihoodsFatherConstantPart_i_c)[x] *= freqs[x];
}
}
}
}

// Then, we deal with the node of interest.
// We first average uppon 'y' to save computations, and then uppon 'x'.
// ('y' is the state at 'node' and 'x' the state at 'father'.)

// Iterate over all site partitions:
const VVVdouble* likelihoodsFather_node = &drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentNode->getId());
unique_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(currentNode->getId()));
VVVdouble pxy;
bool first;
while (mit->hasNext())
{
TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
substitutionCounts.setSubstitutionModel(bmd->getSubstitutionModel());
// compute all nxy first:
VVVVdouble nxy(nbClasses);
for (size_t c = 0; c < nbClasses; ++c)
{
double rc = rcRates[c];
VVVdouble* nxy_c = &nxy[c];
nxy_c->resize(nbTypes);
for (size_t t = 0; t < nbTypes; ++t)
{
VVdouble* nxy_c_t = &(*nxy_c)[t];
nxy_c_t->resize(nbStates);
Matrix<double>* nijt = substitutionCounts.getAllNumbersOfSubstitutions(d * rc, t + 1);
for (size_t x = 0; x < nbStates; ++x)
{
Vdouble* nxy_c_t_x = &(*nxy_c_t)[x];
nxy_c_t_x->resize(nbStates);
for (size_t y = 0; y < nbStates; ++y)
{
(*nxy_c_t_x)[y] = (*nijt)(x, y);
}
}
delete nijt;
}
}

// Now loop over sites:
unique_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
first = true;
while (sit->hasNext())
{
size_t i = sit->next();
// We retrieve the transition probabilities and substitution counts for this site partition:
if (first)
{
pxy = drtl.getTransitionProbabilitiesPerRateClass(currentNode->getId(), i);
first = false;
}
const VVdouble* likelihoodsFather_node_i = &(*likelihoodsFather_node)[i];
VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
RowMatrix<double> pairProbabilities(nbStates, nbStates);
MatrixTools::fill(pairProbabilities, 0.);
VVVdouble subsCounts(nbStates);
for (size_t j = 0; j < nbStates; ++j)
{
subsCounts[j].resize(nbStates);
for (size_t k = 0; k < nbStates; ++k)
{
subsCounts[j][k].resize(nbTypes);
}
}
for (size_t c = 0; c < nbClasses; ++c)
{
const Vdouble* likelihoodsFather_node_i_c = &(*likelihoodsFather_node_i)[c];
Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
const VVdouble* pxy_c = &pxy[c];
VVVdouble* nxy_c = &nxy[c];
for (size_t x = 0; x < nbStates; ++x)
{
double* likelihoodsFatherConstantPart_i_c_x = &(*likelihoodsFatherConstantPart_i_c)[x];
const Vdouble* pxy_c_x = &(*pxy_c)[x];
for (size_t y = 0; y < nbStates; ++y)
{
double likelihood_cxy = (*likelihoodsFatherConstantPart_i_c_x)
* (*pxy_c_x)[y]
* (*likelihoodsFather_node_i_c)[y];
pairProbabilities(x, y) += likelihood_cxy; // Sum over all rate classes.
for (size_t t = 0; t < nbTypes; ++t)
{
subsCounts[x][y][t] += likelihood_cxy * (*nxy_c)[t][x][y];
}
}
}
}
// Now the vector computation:
// Here we do not average over all possible pair of ancestral states,
// We only consider the one with max likelihood:
vector<size_t> xy = MatrixTools::whichMax(pairProbabilities);
for (size_t t = 0; t < nbTypes; ++t)
{
substitutionsForCurrentNode[i][t] += subsCounts[xy[0]][xy[1]][t] / pairProbabilities(xy[0], xy[1]);
}
}
}
// Now we just have to copy the substitutions into the result vector:
for (size_t i = 0; i < nbSites; i++)
{
for (size_t t = 0; t < nbTypes; t++)
{
(*substitutions)(l, i, t) = substitutionsForCurrentNode[(*rootPatternLinks)[i]][t];
}
}
}
if (verbose)
{
if (ApplicationTools::message)
*ApplicationTools::message << " ";
ApplicationTools::displayTaskDone();
}
return substitutions;
}
*/

/**************************************************************************************************/

/*ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(
  const DRTreeLikelihood& drtl,
  SubstitutionCount& substitutionCount,
  bool verbose) throw (Exception)
  {
  // Preamble:
  if (!drtl.isInitialized())
  throw Exception("SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(). Likelihood object is not initialized.");

  // A few variables we'll need:

  const TreeTemplate<Node> tree(drtl.getTree());
  const AlignedValuesContainer*    sequences = drtl.getData();
  const DiscreteDistribution* rDist = drtl.getRateDistribution();
  const Alphabet*             alpha = sequences->getAlphabet();

  size_t nbSites         = sequences->getNumberOfSites();
  size_t nbDistinctSites = drtl.getLikelihoodData()->getNumberOfDistinctSites();
  size_t nbStates        = alpha->getSize();
  size_t nbTypes         = substitutionCounts.getNumberOfSubstitutionTypes();
  vector<const Node*> nodes    = tree.getNodes();
  const vector<size_t>* rootPatternLinks
  = &drtl.getLikelihoodData()->getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  size_t nbNodes = nodes.size();

  // We create a new ProbabilisticSubstitutionMapping object:
  ProbabilisticSubstitutionMapping* substitutions = new ProbabilisticSubstitutionMapping(tree, &substitutionCount, nbSites);

  // Compute the whole likelihood of the tree according to the specified model:

  Vdouble rcRates = rDist->getCategories();

  // Compute the number of substitutions for each class and each branch in the tree:
  if (verbose)
  ApplicationTools::displayTask("Compute marginal ancestral states");
  MarginalAncestralStateReconstruction masr(&drtl);
  map<int, vector<size_t> > ancestors = masr.getAllAncestralStates();
  if (verbose)
  ApplicationTools::displayTaskDone();

  // Now we just have to compute the substitution vectors:
  if (verbose)
  ApplicationTools::displayTask("Compute substitution vectors", true);

  for (size_t l = 0; l < nbNodes; l++)
  {
  const Node* currentNode = nodes[l];

  const Node* father = currentNode->getFather();

  double d = currentNode->getDistanceToFather();

  vector<size_t> nodeStates = ancestors[currentNode->getId()]; // These are not 'true' ancestors ;)
  vector<size_t> fatherStates = ancestors[father->getId()];

  // For each node,
  if (verbose)
  ApplicationTools::displayGauge(l, nbNodes - 1);
  VVdouble substitutionsForCurrentNode(nbDistinctSites);
  for (size_t i = 0; i < nbDistinctSites; ++i)
  {
  substitutionsForCurrentNode[i].resize(nbTypes);
  }

  // Here, we have no likelihood computation to do!

  // Then, we deal with the node of interest.
  // ('y' is the state at 'node' and 'x' the state at 'father'.)
  // Iterate over all site partitions:
  unique_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(currentNode->getId()));
  while (mit->hasNext())
  {
  TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
  substitutionCounts.setSubstitutionModel(bmd->getSubstitutionModel());
  // compute all nxy first:
  VVVdouble nxyt(nbTypes);
  for (size_t t = 0; t < nbTypes; ++t)
  {
  nxyt[t].resize(nbStates);
  Matrix<double>* nxy = substitutionCounts.getAllNumbersOfSubstitutions(d, t + 1);
  for (size_t x = 0; x < nbStates; ++x)
  {
  nxyt[t][x].resize(nbStates);
  for (size_t y = 0; y < nbStates; ++y)
  {
  nxyt[t][x][y] = (*nxy)(x, y);
  }
  }
  delete nxy;
  }
  // Now loop over sites:
  unique_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
  while (sit->hasNext())
  {
  size_t i = sit->next();
  size_t fatherState = fatherStates[i];
  size_t nodeState   = nodeStates[i];
  if (fatherState >= nbStates || nodeState >= nbStates)
  for (size_t t = 0; t < nbTypes; ++t)
  {
  substitutionsForCurrentNode[i][t] = 0;
  }                                                    // To be conservative! Only in case there are generic characters.
  else
  for (size_t t = 0; t < nbTypes; ++t)
  {
  substitutionsForCurrentNode[i][t] = nxyt[t][fatherState][nodeState];
  }
  }
  }

  // Now we just have to copy the substitutions into the result vector:
  for (size_t i = 0; i < nbSites; i++)
  {
  for (size_t t = 0; t < nbTypes; t++)
  {
  (*substitutions)(l, i, t) = substitutionsForCurrentNode[(*rootPatternLinks)[i]][t];
  }
  }
  }
  if (verbose)
  {
  if (ApplicationTools::message)
  *ApplicationTools::message << " ";
  ApplicationTools::displayTaskDone();
  }
  return substitutions;
  }
*/
/**************************************************************************************************/

/*ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeSubstitutionVectorsMarginal(
  const DRTreeLikelihood& drtl,
  SubstitutionCount& substitutionCount,
  bool verbose) throw (Exception)
  {
  // Preamble:
  if (!drtl.isInitialized())
  throw Exception("SubstitutionMappingTools::computeSubstitutionVectorsMarginal(). Likelihood object is not initialized.");

  // A few variables we'll need:

  const TreeTemplate<Node> tree(drtl.getTree());
  const AlignedValuesContainer*    sequences = drtl.getData();
  const DiscreteDistribution* rDist = drtl.getRateDistribution();

  size_t nbSites         = sequences->getNumberOfSites();
  size_t nbDistinctSites = drtl.getLikelihoodData()->getNumberOfDistinctSites();
  size_t nbStates        = sequences->getAlphabet()->getSize();
  size_t nbClasses       = rDist->getNumberOfCategories();
  size_t nbTypes         = substitutionCounts.getNumberOfSubstitutionTypes();
  vector<const Node*> nodes    = tree.getNodes();
  const vector<size_t>* rootPatternLinks
  = &drtl.getLikelihoodData()->getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  size_t nbNodes = nodes.size();

  // We create a new ProbabilisticSubstitutionMapping object:
  ProbabilisticSubstitutionMapping* substitutions = new ProbabilisticSubstitutionMapping(tree, &substitutionCount, nbSites);

  // Compute the whole likelihood of the tree according to the specified model:

  Vdouble rcProbs = rDist->getProbabilities();
  Vdouble rcRates = rDist->getCategories();

  // II) Compute the number of substitutions for each class and each branch in the tree:
  if (verbose)
  ApplicationTools::displayTask("Compute marginal node-pairs likelihoods", true);

  for (size_t l = 0; l < nbNodes; l++)
  {
  const Node* currentNode = nodes[l];

  const Node* father = currentNode->getFather();

  double d = currentNode->getDistanceToFather();

  // For each node,
  if (verbose)
  ApplicationTools::displayGauge(l, nbNodes - 1);
  VVdouble substitutionsForCurrentNode(nbDistinctSites);
  for (size_t i = 0; i < nbDistinctSites; ++i)
  {
  substitutionsForCurrentNode[i].resize(nbTypes);
  }

  // Then, we deal with the node of interest.
  // ('y' is the state at 'node' and 'x' the state at 'father'.)
  VVVdouble probsNode   = DRTreeLikelihoodTools::getPosteriorProbabilitiesPerStatePerRate(drtl, currentNode->getId());
  VVVdouble probsFather = DRTreeLikelihoodTools::getPosteriorProbabilitiesPerStatePerRate(drtl, father->getId());

  // Iterate over all site partitions:
  unique_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(currentNode->getId()));
  while (mit->hasNext())
  {
  TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
  substitutionCounts.setSubstitutionModel(bmd->getSubstitutionModel());
  // compute all nxy first:
  VVVVdouble nxy(nbClasses);
  for (size_t c = 0; c < nbClasses; ++c)
  {
  VVVdouble* nxy_c = &nxy[c];
  double rc = rcRates[c];
  nxy_c->resize(nbTypes);
  for (size_t t = 0; t < nbTypes; ++t)
  {
  VVdouble* nxy_c_t = &(*nxy_c)[t];
  Matrix<double>* nijt = substitutionCounts.getAllNumbersOfSubstitutions(d * rc, t + 1);
  nxy_c_t->resize(nbStates);
  for (size_t x = 0; x < nbStates; ++x)
  {
  Vdouble* nxy_c_t_x = &(*nxy_c_t)[x];
  nxy_c_t_x->resize(nbStates);
  for (size_t y = 0; y < nbStates; ++y)
  {
  (*nxy_c_t_x)[y] = (*nijt)(x, y);
  }
  }
  delete nijt;
  }
  }

  // Now loop over sites:
  unique_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
  while (sit->hasNext())
  {
  size_t i = sit->next();
  VVdouble* probsNode_i   = &probsNode[i];
  VVdouble* probsFather_i = &probsFather[i];
  for (size_t c = 0; c < nbClasses; ++c)
  {
  Vdouble* probsNode_i_c   = &(*probsNode_i)[c];
  Vdouble* probsFather_i_c = &(*probsFather_i)[c];
  VVVdouble* nxy_c = &nxy[c];
  for (size_t x = 0; x < nbStates; ++x)
  {
  for (size_t y = 0; y < nbStates; ++y)
  {
  double prob_cxy = (*probsFather_i_c)[x] * (*probsNode_i_c)[y];
  // Now the vector computation:
  for (size_t t = 0; t < nbTypes; ++t)
  {
  substitutionsForCurrentNode[i][t] += prob_cxy * (*nxy_c)[t][x][y];
  //                                   <------>   <--------------->
  // Posterior probability                 |                |
  // for site i and rate class c *         |                |
  // likelihood for this site--------------+                |
  //                                                        |
  // Substitution function for site i and rate class c-------+
  }
  }
  }
  }
  }
  }

  // Now we just have to copy the substitutions into the result vector:
  for (size_t i = 0; i < nbSites; ++i)
  {
  for (size_t t = 0; t < nbTypes; ++t)
  {
  (*substitutions)(l, i, t) = substitutionsForCurrentNode[(*rootPatternLinks)[i]][t];
  }
  }
  }
  if (verbose)
  {
  if (ApplicationTools::message)
  *ApplicationTools::message << " ";
  ApplicationTools::displayTaskDone();
  }
  return substitutions;
  }
*/
