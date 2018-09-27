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
#include "ProbabilisticRewardMappingForASite.h"
#include "ProbabilisticSubstitutionMappingForASite.h"
#include "RewardMappingToolsForASite.h"

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

ProbabilisticSubstitutionMappingForASite* SubstitutionMappingToolsForASite::computeCounts(
  size_t site,
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
    throw Exception("SubstitutionMappingToolsForASite::computeSubstitutionVectors(). Likelihood object is not initialized.");

  const SubstitutionModel* sm(0);
  const SubstitutionProcess& sp=*rltc.getSubstitutionProcess();

  if (nodeIds.size()==0)
  {
    const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();
    return new ProbabilisticSubstitutionMappingForASite(ppt, reg.getNumberOfSubstitutionTypes());    
  }
  
  for (auto id : nodeIds)
  {
    if (dynamic_cast<const SubstitutionModel*>(sp.getModel(id,0))==NULL)
      throw Exception("SubstitutionMappingToolsForASite::computeSubstitutionVectors possible only for SubstitutionModels, not in branch " + TextTools::toString(id));
    else
      if (!sm)
        sm=dynamic_cast<const SubstitutionModel*>(sp.getModel(id,0));
  }

  // We create a Mapping objects

  if (!sm)
    throw Exception("SubstitutionMappingToolsForASite::computeSubstitutionVectors not possible with null model.");

  unique_ptr<SubstitutionCount> substitutionCount(new DecompositionSubstitutionCount(sm, reg.clone(), weights, distances));

  return computeCounts(site, rltc, nodeIds, *substitutionCount, threshold, verbose);
}


ProbabilisticSubstitutionMappingForASite* SubstitutionMappingToolsForASite::computeCounts(
  size_t site,
  RecursiveLikelihoodTreeCalculation& rltc,
  const vector<uint>& nodeIds,
  SubstitutionCount& substitutionCount,
  double threshold,
  bool verbose)
{
  // Preamble:
  if (!rltc.isInitialized())
    throw Exception("SubstitutionMappingToolsForASite::computeSubstitutionVectors(). Likelihood object is not initialized.");
  rltc.computeTreeLikelihood();
  
  const SubstitutionProcess& sp=*rltc.getSubstitutionProcess();

  if (nodeIds.size()==0)
  {
    const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();
    return new ProbabilisticSubstitutionMappingForASite(ppt, substitutionCount.getNumberOfSubstitutionTypes());
  }

  for (auto id :nodeIds)
  {
    if (dynamic_cast<const SubstitutionModel*>(sp.getModel(id,0))==NULL)
      throw Exception("SubstitutionMappingToolsForASite::computeSubstitutionVectors possible only for SubstitutionModels, not in branch " + TextTools::toString(id));
  }

  const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();
  
  // A few variables we'll need:

  const RecursiveLikelihoodTree& rlt=dynamic_cast<const RecursiveLikelihoodTree&>(rltc.getLikelihoodData());

  size_t nbStates        = sp.getNumberOfStates();
  size_t nbClasses       = sp.getNumberOfClasses();

  size_t nbTypes         = substitutionCount.getNumberOfSubstitutionTypes();
  size_t nbNodes         = nodeIds.size();
  
  // We create a Mapping objects
  
  unique_ptr<ProbabilisticSubstitutionMappingForASite> substitutions(new ProbabilisticSubstitutionMappingForASite(ppt, nbTypes));

  // Find the corresponding indec in the compressed array :

  size_t siteIndex=rltc.getLikelihoodData().getRootArrayPosition(site);
  
  double Lr=rltc.getLogLikelihoodForASiteIndex(siteIndex);

  // Compute the number of substitutions for each class and each branch in the tree:
  if (verbose)
    ApplicationTools::displayTask("Compute counts", true);

  unique_ptr<ProbabilisticSubstitutionMappingForASite::mapTree::EdgeIterator> brIt=substitutions->allEdgesIterator();

  size_t nn=0;
  
  for (;!brIt->end();brIt->next())
  {
    if (verbose)
      ApplicationTools::displayGauge(nn++, nbNodes - 1);
    
    shared_ptr<PhyloBranchMappingForASite> br=**brIt;
    
    // For each branch
    uint edid=substitutions->getEdgeIndex(br);

    if (nodeIds.size() > 0 && !VectorTools::contains(nodeIds, (int)edid))
      continue;

    uint fathid=substitutions->getFather(edid);
    uint icid=substitutions->getSon(edid);

    double d=br->getLength();
    
    Vdouble likelihoodsFatherConstantPart(nbStates);

    // now the counts
      
    Vdouble substitutionsForCurrentNode(nbTypes);

    bool usesLog=false;
    
    for (size_t ncl=0; ncl<nbClasses; ncl++)
    {
      const RecursiveLikelihoodTree::LikTree& rlt_c=rlt[ncl];
      
      shared_ptr<RecursiveLikelihoodNode> ici = rlt_c.getNode(icid);

      // reinit substitutionsForCurrentNode for log 
      if (!usesLog && ici->usesLog())
      {
        std::fill(substitutionsForCurrentNode.begin(), substitutionsForCurrentNode.end(), NumConstants::MINF());
      }
        
      usesLog=ici->usesLog();

      double pr=sp.getProbabilityForModel(ncl);
      double rate=sp.getRateForModel(ncl);
      
      VectorTools::fill(likelihoodsFatherConstantPart,usesLog?log(pr):pr);

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
              likelihoodsFatherConstantPart *= rlt_c.getSon(**brothIt)->getToFatherBelowLikelihoodArray(ComputingNode::D0)[siteIndex];
            else
              likelihoodsFatherConstantPart *= VectorTools::exp(rlt_c.getSon(**brothIt)->getToFatherBelowLikelihoodArray(ComputingNode::D0)[siteIndex]);
          }
          else {
            if (slog)
              likelihoodsFatherConstantPart += rlt_c.getSon(**brothIt)->getToFatherBelowLikelihoodArray(ComputingNode::D0)[siteIndex];
            else
              likelihoodsFatherConstantPart += VectorTools::log(rlt_c.getSon(**brothIt)->getToFatherBelowLikelihoodArray(ComputingNode::D0)[siteIndex]);
          }
        }
      }

      bool flog=father->usesLog();

      if (!usesLog)
      {
        if (!flog)
          likelihoodsFatherConstantPart *= father->getAboveLikelihoodArray()[siteIndex];
        else
          likelihoodsFatherConstantPart *= VectorTools::exp(father->getAboveLikelihoodArray()[siteIndex]);
      }
      else {
        if (flog)
          likelihoodsFatherConstantPart += father->getAboveLikelihoodArray()[siteIndex];
        else
          likelihoodsFatherConstantPart += VectorTools::log(father->getAboveLikelihoodArray()[siteIndex]);
      }
      
      // Then, we deal with the node of interest.
      // We first average upon 'y' to save computations, and then upon 'x'.
      // ('y' is the state at 'node' and 'x' the state at 'father'.)

      const Vdouble& likelihoodsFather_node = ici->getBelowLikelihoodArray(ComputingNode::D0)[siteIndex];

      const SubstitutionModel* sm=dynamic_cast<const SubstitutionModel*>(sp.getModel(icid,ncl));
      if (!sm)
        throw Exception("SubstitutionMappingToolsForASite:: non substitution model in node " + TextTools::toString(icid));
      
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
                substitutionsForCurrentNode[t] += likelihood_xy * npxy[t][x][y];
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
                if (npxy[t][x][y]< -NumConstants::MILLI())
                {
                  ApplicationTools::displayWarning("These counts are negative, their logs could not be computed:" + TextTools::toString(npxy[t][x][y]));
                  throw Exception("Stop in SubstitutionMappingToolsForASite");
                }
                else
                {
                  if (npxy[t][x][y]>0)
                  {
                    double ll=likelihood_xy + log(npxy[t][x][y]);
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
        (*br)(t)=0;
      }
      else
      {
        if (threshold>=0 && x > threshold)
        {
          if (verbose)
            ApplicationTools::displayWarning("On branch " + TextTools::toString(edid) + ", site " + TextTools::toString(site) + ", and type " + TextTools::toString(t) + " count has been ignored because it is presumably saturated.");
          (*br)(t)=0;
        }
        else     
          (*br)(t)=x;;
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

ProbabilisticSubstitutionMappingForASite* SubstitutionMappingToolsForASite::computeNormalizations(
  size_t site,
  RecursiveLikelihoodTreeCalculation& rltc,
  const vector<uint>& nodeIds,
  const BranchedModelSet* nullModels,
  const SubstitutionRegister& reg,
  std::shared_ptr<const AlphabetIndex2> distances,
  bool verbose)
{
  // Preamble:
  if (!rltc.isInitialized())
    throw Exception("SubstitutionMappingToolsForASite::computeNormalizations(). Likelihood object is not initialized.");
  rltc.computeTreeLikelihood();

  const SubstitutionModel* sm(0);
  const SubstitutionProcess& sp=*rltc.getSubstitutionProcess();

  if (nodeIds.size()==0)
  {
    const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();
    return new ProbabilisticSubstitutionMappingForASite(ppt, reg.getNumberOfSubstitutionTypes());    
  }
  
  for (auto id :nodeIds)
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

  unique_ptr<ProbabilisticSubstitutionMappingForASite> normalizations(new ProbabilisticSubstitutionMappingForASite(ppt, nbTypes));
  
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
        
        unique_ptr<ProbabilisticRewardMappingForASite> mapping(RewardMappingToolsForASite::computeRewardVectors(site, rltc, mids, *reward, verbose));

        for (size_t k = 0; k < mids.size(); k++)
        {
          shared_ptr<PhyloBranchMappingForASite> brn=normalizations->getEdge(mids[k]);
          shared_ptr<PhyloBranchRewardForASite> brr=mapping->getEdge(mids[k]);

          (*brn)(nbt)=(*brr).getReward();
        }
      }
    }
  }

  return normalizations.release();
}

/************************************************************/

ProbabilisticSubstitutionMappingForASite* SubstitutionMappingToolsForASite::computeNormalizedCounts(
  size_t site,
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
  unique_ptr<ProbabilisticSubstitutionMappingForASite> counts(computeCounts(site,rltc,nodeIds,reg,weights,distances,threshold,verbose));
  
  unique_ptr<ProbabilisticSubstitutionMappingForASite> factors(computeNormalizations(site,rltc,nodeIds,nullModels,reg,distances,verbose));
  
  return computeNormalizedCounts(counts.get(), factors.get(), nodeIds, perTimeUnit, siteSize);  
}
    
/************************************************************/

ProbabilisticSubstitutionMappingForASite* SubstitutionMappingToolsForASite::computeNormalizedCounts(
  const ProbabilisticSubstitutionMappingForASite* counts,
  const ProbabilisticSubstitutionMappingForASite* factors,
  const vector<uint>& nodeIds,
  bool perTimeUnit,
  uint siteSize)
{
  unique_ptr<ProbabilisticSubstitutionMappingForASite> normCounts(counts->clone());
  
  size_t nbTypes = counts->getNumberOfTypes();

  // Iterate on branches
  
  unique_ptr<ProbabilisticSubstitutionMappingForASite::mapTree::EdgeIterator> brIt=normCounts->allEdgesIterator();

  for (;!brIt->end();brIt->next())
  {
    shared_ptr<PhyloBranchMappingForASite> brNormCount=**brIt;

    Vdouble& brnCou=brNormCount->getCounts();

    // For each branch
    uint edid=normCounts->getEdgeIndex(brNormCount);

    if (nodeIds.size() > 0 && !VectorTools::contains(nodeIds, (int)edid))
    {
      for (auto& brk : brnCou)
        brk=0;
      continue;
    }
        
    shared_ptr<PhyloBranchMappingForASite> brFactor=factors->getEdge(edid);
    shared_ptr<PhyloBranchMappingForASite> brCount=counts->getEdge(edid);

    const Vdouble& cou=brCount->getCounts();
    const Vdouble& fac=brFactor->getCounts();


    // if not per time, multiply by the lengths of the branches of
    // the input tree
    
    double slg=(!perTimeUnit?brCount->getLength():1)/siteSize;
    
    for (size_t t = 0; t < nbTypes; ++t)
      brnCou[t]=(fac[t]!=0? cou[t]/fac[t]*slg : 0);
  }

  return normCounts.release();
}

/*******************************************************************************/
/* Get trees of counts */
/*******************************************************************************/

PhyloTree* SubstitutionMappingToolsForASite::getTreeForType(const ProbabilisticSubstitutionMappingForASite& counts,
                                                            size_t type)
{
  unique_ptr<PhyloTree> pt(new PhyloTree(counts));

  unique_ptr<ProbabilisticSubstitutionMappingForASite::mapTree::EdgeIterator> brIt=counts.allEdgesIterator();

  for (;!brIt->end();brIt->next())
  {
    shared_ptr<PhyloBranchMappingForASite> brm=(**brIt);
    pt->getEdge(counts.getEdgeIndex(brm))->setLength(brm->getTypeCount(type));
  }
 
  return pt.release();
}

PhyloTree* SubstitutionMappingToolsForASite::getTreeForType(const ProbabilisticSubstitutionMappingForASite& counts,
                                                            const ProbabilisticSubstitutionMappingForASite& factors,
                                                            size_t type)
{
  unique_ptr<PhyloTree> pt(new PhyloTree(counts));

  unique_ptr<ProbabilisticSubstitutionMappingForASite::mapTree::EdgeIterator> brIt=counts.allEdgesIterator();

  for (;!brIt->end();brIt->next())
  {
    shared_ptr<PhyloBranchMappingForASite> brm=(**brIt);
    shared_ptr<PhyloBranchMappingForASite> brf=factors.getEdge(counts.getEdgeIndex(brm));
    
    pt->getEdge(counts.getEdgeIndex(brm))->setLength(brm->getTypeCount(type)/brf->getTypeCount(type));
  }
 
  return pt.release();
}
  

/********************************************************************************/
/*  Get vectors of counts  */
/********************************************************************************/

VVdouble SubstitutionMappingToolsForASite::getCountsPerBranchPerType(
  const ProbabilisticSubstitutionMappingForASite& counts,
  const vector<uint>& ids)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);
  size_t nbBr = idc.size();
  
  size_t nbTypes= counts.getNumberOfTypes();

  VVdouble result;
  VectorTools::resize2(result, nbBr, nbTypes);

  for (size_t k = 0; k < nbBr; ++k)
  {
    Vdouble& resB=result[k];

    for (size_t i = 0; i < nbTypes; ++i)
      resB[i] = counts(idc[k],i);
  }

  return result;
}

/**************************************************************************************************/

Vdouble SubstitutionMappingToolsForASite::getCountsPerBranch(
  const ProbabilisticSubstitutionMappingForASite& counts,
  const vector<uint>& ids)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);
  size_t nbBr = idc.size();

  Vdouble result(nbBr);

  vector<double> countsf(SubstitutionMappingToolsForASite::getCountsPerBranch(counts));

  for (size_t i = 0; i < nbBr; ++i)
    result[i]= countsf[idc[i]];

  return result;
}

Vdouble SubstitutionMappingToolsForASite::getCountsPerBranch(
  const ProbabilisticSubstitutionMappingForASite& counts,
  const ProbabilisticSubstitutionMappingForASite& factors,
  const vector<uint>& ids)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);
  size_t nbBr = idc.size();

  Vdouble result(nbBr);

  vector<double> countsf(SubstitutionMappingToolsForASite::getCountsPerBranch(counts, factors));

  for (size_t i = 0; i < nbBr; ++i)
    result[i]= countsf[idc[i]];

  return result;
}

/**************************************************************************************************/

Vdouble SubstitutionMappingToolsForASite::getCountsForBranchPerType(
  const ProbabilisticSubstitutionMappingForASite& counts,
  uint branchId)
{
  size_t nbTypes = counts.getNumberOfTypes();
  Vdouble v(nbTypes, 0);
  shared_ptr<PhyloBranchMappingForASite> br=counts.getEdge(branchId);

  for (size_t t = 0; t < nbTypes; ++t)
    v[t] = br->getTypeCount(t);
  
  return v;
}

/**************************************************************************************************/

VVdouble SubstitutionMappingToolsForASite::getCountsPerTypePerBranch(
  const ProbabilisticSubstitutionMappingForASite& counts,
  const vector<uint>& ids)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);
  size_t nbBr = idc.size();

  size_t nbTypes = counts.getNumberOfTypes();

  VVdouble result;
  VectorTools::resize2(result, nbTypes, nbBr);

  for (size_t i=0;i<idc.size();i++)
  {
    Vdouble cou=SubstitutionMappingToolsForASite::getCountsForBranchPerType(counts, idc[i]);
    for (size_t nbt=0;nbt<nbTypes;nbt++)
      result[nbt][i]=cou[nbt];
  }
  
  return result;
}

/**************************************************************************************************/

Vdouble SubstitutionMappingToolsForASite::getCountsPerType(const ProbabilisticSubstitutionMappingForASite& counts, const vector<uint>& ids)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);

  size_t nbTypes    = counts.getNumberOfTypes();
  Vdouble v(nbTypes,0);

  for (auto id : idc)
  {
    shared_ptr<PhyloBranchMappingForASite> br=counts.getEdge(id);
    for (size_t t = 0; t < nbTypes; ++t)
      v[t] += (*br)(t);
  }
  
  return v;
}

/**************************************************************************************************/

Vdouble SubstitutionMappingToolsForASite::getCountsPerType(
  const ProbabilisticSubstitutionMappingForASite& counts,
  const ProbabilisticSubstitutionMappingForASite& factors,
  const vector<uint>& ids,
  bool perTimeUnit,
  uint siteSize)
{
  const Vuint idc(ids.size()==0?counts.getAllEdgesIndexes():ids);

  size_t nbTypes   = counts.getNumberOfTypes();
  
  Vdouble v(nbTypes,0);
  Vdouble n(nbTypes,0);

  double lg=(!perTimeUnit?0:1);
  for (auto id : idc)
  {
    shared_ptr<PhyloBranchMappingForASite> br=counts.getEdge(id);
    for (size_t t = 0; t < nbTypes; ++t)
      v[t] += (*br)(t);
    if (!perTimeUnit)
      lg+=br->getLength();
  }

  double slg=lg/siteSize;
  
  for (auto id : idc)
  {
    shared_ptr<PhyloBranchMappingForASite> br=factors.getEdge(id);
    for (size_t t = 0; t < nbTypes; ++t)
      n[t] += (*br)(t);
  }
  
  for (size_t t = 0; t < nbTypes; ++t)
    v[t] = v[t]/n[t]*slg;

  return v;
}



/**************************************************************************************************/

void SubstitutionMappingToolsForASite::writeToStream(
  const ProbabilisticSubstitutionMappingForASite& substitutions,
  size_t type,
  ostream& out)
{
  if (!out)
    throw IOException("SubstitutionMappingToolsForASite::writeToFile. Can't write to stream.");
  out << "Branches";
  out << "\tMean";
  out << "\tSite";
  out << endl;

  unique_ptr<ProbabilisticSubstitutionMappingForASite::mapTree::EdgeIterator> brIt=substitutions.allEdgesIterator();

  for (;!brIt->end();brIt->next())
  {
    const shared_ptr<PhyloBranchMappingForASite> br=**brIt;
    
    out << substitutions.getEdgeIndex(br) << "\t" << br->getLength();

    out << "\t" << (*br)(type);
    
    out << endl;
  }
}

