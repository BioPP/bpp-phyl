//
// File: RewardMappingToolsForASite.cpp
// Created by: Laurent Guéguen
// Created on: mercredi 26 septembre 2018, à 16h 01
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

#include "RewardMappingToolsForASite.h"
#include "../Likelihood/DRTreeLikelihoodTools.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/DataTable.h>

using namespace bpp;

// From the STL:
#include <iomanip>

using namespace std;


/******************************************************************************/

ProbabilisticRewardMapping* RewardMappingToolsForASite::computeRewardVectors(
  size_t site,
  LikelihoodCalculationSingleProcess& rltc,
  const vector<uint>& nodeIds,
  Reward& reward,
  bool verbose)
{
  // Preamble:
  /*
    if (!rltc.isInitialized())
    throw Exception("RewardMappingToolsForASite::computeSubstitutionVectors(). Likelihood object is not initialized.");
  rltc.computeTreeLikelihood();

  const SubstitutionProcess& sp=*rltc.getSubstitutionProcess();
  const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();
  
  // A few variables we'll need:

  const RecursiveLikelihoodTree& rlt=dynamic_cast<const RecursiveLikelihoodTree&>(rltc.getLikelihoodData());

  size_t nbStates        = sp.getNumberOfStates();
  size_t nbClasses       = sp.getNumberOfClasses();
  size_t nbNodes         = nodeIds.size();
  
  // We create a new ProbabilisticRewardMapping object:
  unique_ptr<ProbabilisticRewardMapping> rewards(new ProbabilisticRewardMapping(ppt, 1));
  
  // Store likelihood for each site (here rootPatterns are managed):

  size_t siteIndex=rltc.getLikelihoodData().getRootArrayPosition(site);

  double Lr=rltc.getLogLikelihoodForASiteIndex(siteIndex);

  // Compute the reward for each class and each branch in the tree:
  if (verbose)
    ApplicationTools::displayTask("Compute rewards", true);

  unique_ptr<ProbabilisticRewardMapping::mapTree::EdgeIterator> brIt=rewards->allEdgesIterator();

  size_t nn=0;
  
  Vdouble likelihoodsFatherConstantPart(nbStates);

  for (;!brIt->end();brIt->next())
  {
    if (verbose)
      ApplicationTools::displayGauge(nn++, nbNodes - 1);
    
    shared_ptr<PhyloBranchReward> br=**brIt;
    
    // For each branch
    uint edid=rewards->getEdgeIndex(br);

    if (nodeIds.size() > 0 && !VectorTools::contains(nodeIds, (int)edid))
      continue;

    uint fathid=rewards->getFatherOfEdge(edid);
    uint icid=rewards->getSon(edid);

    double d=br->getLength();

    double rewardsForCurrentNode(0);

    bool usesLog=false;

    for (size_t ncl=0; ncl<nbClasses; ncl++)
    {
      const RecursiveLikelihoodTree::LikTree& rlt_c=rlt[ncl];
      
      shared_ptr<RecursiveLikelihoodNode> ici = rlt_c.getNode(icid);

      // reinit substitutionsForCurrentNode for log 
      if (!usesLog && ici->usesLog())
        rewardsForCurrentNode=NumConstants::MINF();
      
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
        throw Exception("RewardMappingToolsForASite:: non substitution model in node " + TextTools::toString(icid));
      
      reward.setSubstitutionModel(sm);
      
      // compute all nxy * pxy first:
      
      const Matrix<double>& pxy = sp.getTransitionProbabilities(icid,ncl);
      Matrix<double>* nij = reward.getAllRewards(d * rate);
      MatrixTools::hadamardMult((*nij),pxy,(*nij));

      
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
              rewardsForCurrentNode += likelihood_xy * (*nij)(x,y);
            //                       <------------>   <--------------->
            // Posterior probability         |                 |
            // for  rate class c *           |                 |
            // likelihood              ------+                 |
            //                                                 |
            // Reward function and rate class c ---------------+
          }
          else
          {
            if (likelihood_xy!=NumConstants::MINF())  // to avoid add per -inf
            {
              if ((*nij)(x,y)< -NumConstants::MILLI())
              {
                ApplicationTools::displayWarning("These rewards are negative, their logs could not be computed:" + TextTools::toString((*nij)(x,y)));
                throw Exception("Stop in RewardMappingToolsForASite");
              }
              else
              {
                if ((*nij)(x,y)>0)
                {
                  double ll=likelihood_xy + log((*nij)(x,y));
                  if (ll>rewardsForCurrentNode)
                    rewardsForCurrentNode = ll + log(1 + exp(rewardsForCurrentNode - ll));
                  else
                    rewardsForCurrentNode += log(1 + exp(ll - rewardsForCurrentNode));
                }
              }
            }
          }
        }
      }
    
    delete nij;
    }
    
    // Now we just have to copy the substitutions into the result vector:
    (*br)(0) = usesLog?exp(rewardsForCurrentNode - Lr):
      rewardsForCurrentNode / exp(Lr);

  }
  if (verbose)
  {
    if (ApplicationTools::message)
      *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }
  return rewards.release();
  */
  return 0;
}

