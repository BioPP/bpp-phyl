// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Text/TextTools.h>

#include "../Likelihood/DataFlow/ForwardLikelihoodTree.h"
#include "../Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h"
#include "RewardMappingTools.h"

using namespace bpp;
using namespace numeric;

// From the STL:
#include <iomanip>

using namespace std;


/******************************************************************************/

unique_ptr<ProbabilisticRewardMapping> RewardMappingTools::computeRewardVectors(
  LikelihoodCalculationSingleProcess& rltc,
  const vector<uint>& edgeIds,
  Reward& reward,
  short unresolvedOption,
  bool verbose)
{
  // Preamble:
  if (!rltc.isInitialized())
    throw Exception("RewardMappingTools::computeSubstitutionVectors(). Likelihood object is not initialized.");

  const SubstitutionProcessInterface& sp = rltc.substitutionProcess();

  if (edgeIds.size() == 0)
    return make_unique<ProbabilisticRewardMapping>(
        sp.parametrizablePhyloTree(),
	rltc.getRootArrayPositions(),
	rltc.getNumberOfDistinctSites());

  auto processTree = rltc.getTreeNode(0);

  /* First, set substitution rewards */

  map<const SubstitutionModelInterface*, shared_ptr<Reward>> mModReward;

  for (auto speciesId : edgeIds)
  {
    const auto& dagIndexes = rltc.getEdgesIds(speciesId, 0);

    for (auto id : dagIndexes)
    {
      const auto& edge = processTree->getEdge(id);
      if (edge->getBrLen()) // if edge with model on it
      {
        auto model = edge->getModel();

        auto nMod = edge->getNMod();

        auto tm = dynamic_pointer_cast<const TransitionModelInterface>(model->targetValue());

        shared_ptr<const SubstitutionModelInterface> sm = nullptr;

        if (nMod == 0)
        {
          sm = dynamic_pointer_cast<const SubstitutionModelInterface>(tm);

          if (!sm)
            throw Exception("SubstitutionMappingTools::computeCounts : SubstitutionVectors possible only for SubstitutionModels, not in branch " + TextTools::toString(speciesId) + ". Got model " + tm->getName());
        }
        else
        {
          size_t nmod = nMod->targetValue();

          auto ttm = dynamic_pointer_cast<const MixedTransitionModelInterface>(tm);
          if (!ttm)
            throw Exception("SubstitutionMappingTools::computeCounts : Expecting Mixed model in branch " + TextTools::toString(speciesId) + ". Got model " + tm->getName());

          sm = dynamic_pointer_cast<const SubstitutionModelInterface>(ttm->getNModel(nmod));

          if (!sm)
            throw Exception("SubstitutionMappingTools::computeCounts : Expecting Substitution model for submodel " + TextTools::toString(nmod) + " of mixed model " + tm->getName() + " in branch " + TextTools::toString(speciesId));
        }

        if (mModReward.find(sm.get()) == mModReward.end())
        {
          mModReward[sm.get()] = shared_ptr<Reward>(reward.clone());
          mModReward[sm.get()]->setSubstitutionModel(sm);
        }
      }
    }
  }

  auto ppt = sp.getParametrizablePhyloTree();

  // A few variables we'll need:

  size_t nbDistinctSites = rltc.getNumberOfDistinctSites();
  size_t nbClasses       = sp.getNumberOfClasses();
  size_t nbNodes         = edgeIds.size();

  const auto& rootPatternLinks = rltc.getRootArrayPositions();

  // We create a new ProbabilisticRewardMapping object:
  unique_ptr<ProbabilisticRewardMapping> rewards(new ProbabilisticRewardMapping(*ppt, rootPatternLinks, nbDistinctSites));

  // Compute the reward for each class and each branch in the tree:

  // Get the DAG of probabilities of the edges
  ProbabilityDAG probaDAG(rltc.getForwardLikelihoodTree(0));

  if (verbose)
    ApplicationTools::displayTask("Compute rewards", true);

  unique_ptr<ProbabilisticRewardMapping::mapTree::EdgeIterator> brIt = rewards->allEdgesIterator();

  Eigen::MatrixXd rpxy;

  size_t nn = 0;
  for ( ; !brIt->end(); brIt->next())
  {
    if (verbose)
      ApplicationTools::displayGauge(nn++, nbNodes - 1);

    shared_ptr<PhyloBranchReward> br = **brIt;

    // For each branch
    uint speciesId = rewards->getEdgeIndex(br);

    if (edgeIds.size() > 0 && !VectorTools::contains(edgeIds, (int)speciesId))
      continue;

    RowLik rewardsForCurrentNode(RowLik::Zero(int(nbDistinctSites)));

    for (size_t ncl = 0; ncl < nbClasses; ncl++)
    {
      processTree = rltc.getTreeNode(ncl);

      double pr = sp.getProbabilityForModel(ncl);

      RowLik rewardsForCurrentClass(RowLik::Zero(Eigen::Index(nbDistinctSites)));

      const auto& dagIndexes = rltc.getEdgesIds(speciesId, ncl);

      // Sum on all dag edges for this speciesId
      for (auto id : dagIndexes)
      {
        auto edge = processTree->getEdge(id);

        auto tm = dynamic_pointer_cast<const TransitionModelInterface>(edge->model().targetValue());

        auto nMod = edge->getNMod();

        shared_ptr<const SubstitutionModelInterface> sm = nullptr;

        if (nMod == 0)
          sm = dynamic_pointer_cast<const SubstitutionModelInterface>(tm);
        else
        {
          size_t nmod = nMod->targetValue();

          auto ttm = dynamic_pointer_cast<const MixedTransitionModelInterface>(tm);
          sm = dynamic_pointer_cast<const SubstitutionModelInterface>(ttm->getNModel(nmod));
        }

        auto subReward = mModReward[sm.get()];

        const auto& likelihoodsTopEdge = rltc.getBackwardLikelihoodsAtEdgeForClass(id, ncl)->targetValue();

        auto sonid = rltc.getForwardLikelihoodTree(ncl)->getSon(id);
        auto fatid = rltc.getForwardLikelihoodTree(ncl)->getFatherOfEdge(id);

        const auto& likelihoodsBotEdge = rltc.getForwardLikelihoodsAtNodeForClass(sonid, ncl)->targetValue();

        // compute all nxy * pxy first:

        const Eigen::MatrixXd& pxy = edge->getTransitionMatrix()->accessValueConst();
        const auto& likelihoodsFather = rltc.getLikelihoodsAtNodeForClass(fatid, ncl)->targetValue();


        subReward->storeAllRewards(edge->getBrLen()->getValue(), rpxy);

        rpxy.array() *= pxy.array();

        // Now loop over sites:

        auto rew = rpxy * likelihoodsBotEdge;

        auto bb = (cwise(likelihoodsTopEdge) * cwise(rew)).colwise().sum();

        
        Eigen::VectorXd ff(likelihoodsBotEdge.cols());
        switch(unresolvedOption){
        case SubstitutionMappingTools::UNRESOLVED_ZERO:
        case SubstitutionMappingTools::UNRESOLVED_AVERAGE:
            
          // Nullify counts where sum likelihoods > 1 : ie unknown
          for (auto i=0;i<ff.size();i++)
          {
            const auto& s=likelihoodsBotEdge.col(i).sum();
            if (s>=2.)
              ff[i]=(unresolvedOption==SubstitutionMappingTools::UNRESOLVED_ZERO)?0.:1./convert(s);
            else
              ff[i]=1;
          }
          
          // Normalizes by likelihood on this node

          bb *= ff.array();
        default:
          ;
        }

        // Normalizes by likelihood on this node
        auto cc = bb / cwise(likelihoodsFather);

        // adds, with branch ponderation  ( * edge / edge * father) probs
        cwise(rewardsForCurrentClass) += cc * probaDAG.getProbaAtNode(fatid);
      }

      rewardsForCurrentNode += rewardsForCurrentClass * pr;
    }

    // Now we just have to copy the substitutions into the result vector:
    for (size_t i = 0; i < nbDistinctSites; ++i)
    {
      (*br)(i) = convert(rewardsForCurrentNode(Eigen::Index(i)));
    }
  }
  if (verbose)
  {
    if (ApplicationTools::message)
      *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }
  return rewards;
}

/**************************************************************************************************/

void RewardMappingTools::writeToStream(
  const ProbabilisticRewardMapping& rewards,
  const AlignmentDataInterface& sites,
  ostream& out)
{
  if (!out)
    throw IOException("RewardMappingTools::writeToFile. Can't write to stream.");
  out << "Branches";
  out << "\tMean";
  for (size_t i = 0; i < rewards.getNumberOfSites(); ++i)
  {
    out << "\tSite" << sites.site(i).getCoordinate();
  }
  out << endl;

  unique_ptr<ProbabilisticRewardMapping::mapTree::EdgeIterator> brIt = rewards.allEdgesIterator();

  for ( ; !brIt->end(); brIt->next())
  {
    const shared_ptr<PhyloBranchReward> br = **brIt;

    out << rewards.getEdgeIndex(br) << "\t" << br->getLength();

    for (size_t i = 0; i < rewards.getNumberOfSites(); ++i)
    {
      out << "\t" << br->getSiteReward(rewards.getSiteIndex(i));
    }

    out << endl;
  }
}

/**************************************************************************************************/

void RewardMappingTools::readFromStream(istream& in, ProbabilisticRewardMapping& rewards)
{
  try
  {
    auto data = DataTable::read(in, "\t", true, -1);
    vector<string> ids = data->getColumn(0);
    data->deleteColumn(0); // Remove ids
    data->deleteColumn(0); // Remove means
    // Now parse the table:
    size_t nbSites = data->getNumberOfColumns();
    rewards.setNumberOfSites(nbSites);

    unique_ptr<ProbabilisticRewardMapping::mapTree::EdgeIterator> brIt = rewards.allEdgesIterator();

    for ( ; !brIt->end(); brIt->next())
    {
      const shared_ptr<PhyloBranchReward> br = **brIt;

      uint brid = rewards.getEdgeIndex(br);
      for (size_t j = 0; j < nbSites; j++)
      {
        (*br)(j) = TextTools::toDouble((*data)(brid, j));
      }
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
  shared_ptr<PhyloBranchReward> br = smap.getEdge((uint)branchIndex);

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
  unique_ptr<ProbabilisticRewardMapping::mapTree::EdgeIterator> brIt = smap.allEdgesIterator();

  size_t siteIndex = smap.getSiteIndex(site);

  for ( ; !brIt->end(); brIt->next())
  {
    v += (**brIt)->getSiteReward(siteIndex);
  }
  return v;
}

/**************************************************************************************************/
