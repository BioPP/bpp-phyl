//
// File: RewardMappingTools.cpp
// Created by: Laurent Guéguen
// Created on: vendredi 29 mars 2013, à 15h 01
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

#include "RewardMappingTools.h"
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

ProbabilisticRewardMapping* RewardMappingTools::computeRewardVectors(
  RecursiveLikelihoodTreeCalculation& rltc,
  const vector<uint>& nodeIds,
  Reward& reward,
  bool verbose)
{
  // Preamble:
  if (!rltc.isInitialized())
    throw Exception("RewardMappingTools::computeSubstitutionVectors(). Likelihood object is not initialized.");
  rltc.computeTreeLikelihood();

  const SubstitutionProcess& sp=*rltc.getSubstitutionProcess();
  const ParametrizablePhyloTree& ppt=sp.getParametrizablePhyloTree();
  
  // A few variables we'll need:

  const RecursiveLikelihoodTree& rlt=dynamic_cast<const RecursiveLikelihoodTree&>(rltc.getLikelihoodData());

  size_t nbDistinctSites = rlt.getNumberOfDistinctSites();
  size_t nbStates        = sp.getNumberOfStates();
  size_t nbClasses       = sp.getNumberOfClasses();
  size_t nbNodes         = nodeIds.size();
  
  const vector<size_t>& rootPatternLinks = rlt.getRootArrayPositions();

  // We create a new ProbabilisticRewardMapping object:
  unique_ptr<ProbabilisticRewardMapping> rewards(new ProbabilisticRewardMapping(ppt, rootPatternLinks, nbDistinctSites));
  
  // Store likelihood for each site (here rootPatterns are managed):
  Vdouble Lr(nbDistinctSites, 0);
  
  for (size_t i = 0; i < nbDistinctSites; i++)
    Lr[i]=rltc.getLogLikelihoodForASiteIndex(i);

  // Compute the reward for each class and each branch in the tree:
  if (verbose)
    ApplicationTools::displayTask("Compute rewards", true);

  unique_ptr<ProbabilisticRewardMapping::mapTree::EdgeIterator> brIt=rewards->allEdgesIterator();

  size_t nn=0;
  
  for (;!brIt->end();brIt->next())
  {
    if (verbose)
      ApplicationTools::displayGauge(nn++, nbNodes - 1);
    
    shared_ptr<PhyloBranchReward> br=**brIt;
    
    // For each branch
    uint edid=rewards->getEdgeIndex(br);

    if (nodeIds.size() > 0 && !VectorTools::contains(nodeIds, (int)edid))
      continue;

    uint fathid=rewards->getFather(edid);
    uint icid=rewards->getSon(edid);

    double d=br->getLength();

    Vdouble rewardsForCurrentNode(nbDistinctSites);

    VVdouble likelihoodsFatherConstantPart;    
    VectorTools::resize2(likelihoodsFatherConstantPart, nbDistinctSites, nbStates);

    bool usesLog=false;

    for (size_t ncl=0; ncl<nbClasses; ncl++)
    {
      const RecursiveLikelihoodTree::LikTree& rlt_c=rlt[ncl];
      
      shared_ptr<RecursiveLikelihoodNode> ici = rlt_c.getNode(icid);

      // reinit substitutionsForCurrentNode for log 
      if (!usesLog && ici->usesLog())
      {
        std::fill(rewardsForCurrentNode.begin(), rewardsForCurrentNode.end(), NumConstants::MINF());
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
        throw Exception("RewardMappingTools:: non substitution model in node " + TextTools::toString(icid));
      
      reward.setSubstitutionModel(sm);
      
      // compute all nxy * pxy first:
      
      const Matrix<double>& pxy = sp.getTransitionProbabilities(icid,ncl);
      Matrix<double>* nij = reward.getAllRewards(d * rate);
      MatrixTools::hadamardMult((*nij),pxy,(*nij));

      
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
                rewardsForCurrentNode[i] += likelihood_xy * (*nij)(x,y);
            //                       <------------>   <--------------->
            // Posterior probability         |                 |
            // for site i and rate class c * |                 |
            // likelihood for this site------+                 |
            //                                                 |
            // Reward function for site i and rate class c------+
            }
            else
            {
              if (likelihood_xy!=NumConstants::MINF())  // to avoid add per -inf
              {
                if ((*nij)(x,y)< -NumConstants::MILLI())
                {
                  ApplicationTools::displayWarning("These rewards are negative, their logs could not be computed:" + TextTools::toString((*nij)(x,y)));
                  throw Exception("Stop in RewardMappingTools");
                }
                else
                {
                  if ((*nij)(x,y)>0)
                  {
                    double ll=likelihood_xy + log((*nij)(x,y));
                    if (ll>rewardsForCurrentNode[i])
                      rewardsForCurrentNode[i] = ll + log(1 + exp(rewardsForCurrentNode[i] - ll));
                    else
                      rewardsForCurrentNode[i] += log(1 + exp(ll - rewardsForCurrentNode[i]));
                  }
                }
              }
            }
          }
        }
      }

      delete nij;
    }
    
    // Now we just have to copy the substitutions into the result vector:
    for (size_t i = 0; i < nbDistinctSites; ++i)
      (*br)(i) = usesLog?exp(rewardsForCurrentNode[i] - Lr[i]):
        rewardsForCurrentNode[i] / exp(Lr[i]);

  }
  if (verbose)
  {
    if (ApplicationTools::message)
      *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }
  return rewards.release();
}

/**************************************************************************************************/

void RewardMappingTools::writeToStream(
  const ProbabilisticRewardMapping& rewards,
  const AlignedValuesContainer& sites,
  ostream& out)
{
  if (!out)
    throw IOException("RewardMappingTools::writeToFile. Can't write to stream.");
  out << "Branches";
  out << "\tMean";
  for (size_t i = 0; i < rewards.getNumberOfSites(); i++)
  {
    out << "\tSite" << sites.getSymbolListSite(i).getPosition();
  }
  out << endl;

  unique_ptr<ProbabilisticRewardMapping::mapTree::EdgeIterator> brIt=rewards.allEdgesIterator();

  for (;!brIt->end();brIt->next())
  {
    const shared_ptr<PhyloBranchReward> br=**brIt;
    
    out << rewards.getEdgeIndex(br) << "\t" << br->getLength();
    
    for (size_t i = 0; i < rewards.getNumberOfSites(); i++)
      out << "\t" << br->getSiteReward(rewards.getSiteIndex(i));
    
    out << endl;
  }
}

/**************************************************************************************************/

void RewardMappingTools::readFromStream(istream& in, ProbabilisticRewardMapping& rewards)
{
  try
  {
    DataTable* data = DataTable::read(in, "\t", true, -1);
    vector<string> ids = data->getColumn(0);
    data->deleteColumn(0); // Remove ids
    data->deleteColumn(0); // Remove means
    // Now parse the table:
    size_t nbSites = data->getNumberOfColumns();
    rewards.setNumberOfSites(nbSites);

    unique_ptr<ProbabilisticRewardMapping::mapTree::EdgeIterator> brIt=rewards.allEdgesIterator();

    for (;!brIt->end();brIt->next())
    {
      const shared_ptr<PhyloBranchReward> br=**brIt;

      uint brid = rewards.getEdgeIndex(br);
      for (size_t j = 0; j < nbSites; j++)
        (*br)(j) = TextTools::toDouble((*data)(brid, j));
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
      rewards.setSitePosition(i, site);
    }

    delete data;
  }
  catch (Exception& e)
  {
    throw IOException(string("Bad input file. ") + e.what());
  }
}

/**************************************************************************************************/

double RewardMappingTools::computeSumForBranch(const ProbabilisticRewardMapping& smap, size_t branchIndex)
{
  size_t nbSites = smap.getNumberOfSites();
  double v = 0;
  shared_ptr<PhyloBranchReward> br=smap.getEdge((uint)branchIndex);
  
  for (size_t i = 0; i < nbSites; ++i)
  {
    v += br->getSiteReward(smap.getSiteIndex(i));
  }
  return v;
}

/**************************************************************************************************/

double RewardMappingTools::computeSumForSite(const ProbabilisticRewardMapping& smap, size_t site)
{
  double v = 0;
  unique_ptr<ProbabilisticRewardMapping::mapTree::EdgeIterator> brIt=smap.allEdgesIterator();

  size_t siteIndex=smap.getSiteIndex(site);
  
  for (;!brIt->end();brIt->next())
  {
    v += (**brIt)->getSiteReward(siteIndex);
  }
  return v;
}

/**************************************************************************************************/
