//
// File: SubstitutionMappingTools.cpp
// Authors:
//   Julien Dutheil, Laurent GuÃÂ©guen
// Created: 2006-04-05 13:04:00
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>
#include <Bpp/Text/TextTools.h>

#include "../Likelihood/DataFlow/ForwardLikelihoodTree.h"
#include "DecompositionReward.h"
#include "DecompositionSubstitutionCount.h"
#include "ProbabilisticRewardMapping.h"
#include "ProbabilisticSubstitutionMapping.h"
#include "RewardMappingTools.h"
#include "SubstitutionMappingTools.h"

using namespace bpp;
using namespace numeric;
using namespace std;

// From the STL:
#include <iomanip>

/******************************************************************************/

unique_ptr<ProbabilisticSubstitutionMapping> SubstitutionMappingTools::computeCounts(
  LikelihoodCalculationSingleProcess& rltc,
  const vector<uint>& edgeIds,
  shared_ptr<const SubstitutionRegisterInterface> reg,
  shared_ptr<const AlphabetIndex2> weights,
  shared_ptr<const AlphabetIndex2> distances,
  short unresolvedOption,
  double threshold,
  bool verbose)
{
  // Preamble:
  if (!rltc.isInitialized())
    throw Exception("SubstitutionMappingTools::computeSubstitutionVectors(). Likelihood object is not initialized.");

  auto& sp = rltc.substitutionProcess();

  if (edgeIds.size() == 0)
    return make_unique<ProbabilisticSubstitutionMapping>(
	sp.parametrizablePhyloTree(),
       	reg->getNumberOfSubstitutionTypes(), 
	rltc.getRootArrayPositions(), 
	rltc.getNumberOfDistinctSites());

  DecompositionSubstitutionCount substitutionCount(reg, weights, distances);
  return computeCounts(rltc, edgeIds, substitutionCount, unresolvedOption, threshold, verbose);
}

/******************************************************************************/

unique_ptr<ProbabilisticSubstitutionMapping> SubstitutionMappingTools::computeCounts(
  LikelihoodCalculationSingleProcess& rltc,
  const vector<uint>& edgeIds,
  SubstitutionCountInterface& substitutionCount,
  short unresolvedOption,
  double threshold,
  bool verbose)
{
  // Preamble:
  if (!rltc.isInitialized())
    throw Exception("SubstitutionMappingTools::computeCounts(). Likelihood object is not initialized.");

  auto& sp = rltc.substitutionProcess();

  if (edgeIds.size() == 0)
    return make_unique<ProbabilisticSubstitutionMapping>(
        sp.parametrizablePhyloTree(),
	substitutionCount.getNumberOfSubstitutionTypes(),
	rltc.getRootArrayPositions(),
	rltc.getNumberOfDistinctSites());

  auto processTree = rltc.getTreeNode(0);

  /* First, set substitution counts */

  // Map from models to counts
  std::map<const SubstitutionModelInterface*, std::shared_ptr<SubstitutionCountInterface> > mModCount;

  for (auto speciesId : edgeIds)
  {
    const auto& dagIndexes = rltc.getEdgesIds(speciesId, 0);

    for (auto id:dagIndexes)
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

        if (mModCount.find(sm.get()) == mModCount.end())
        {
          mModCount[sm.get()] = std::shared_ptr<SubstitutionCountInterface>(substitutionCount.clone());
          mModCount[sm.get()]->setSubstitutionModel(sm);
        }
      }
    }
  }

  auto ppt = sp.getParametrizablePhyloTree();

  // A few variables we'll need:

  size_t nbDistinctSites = rltc.getNumberOfDistinctSites();
  size_t nbClasses       = sp.getNumberOfClasses();

  size_t nbTypes         = substitutionCount.getNumberOfSubstitutionTypes();
  size_t nbNodes         = edgeIds.size();

  const auto& rootPatternLinks = rltc.getRootArrayPositions();

  // We create a Mapping objects

  unique_ptr<ProbabilisticSubstitutionMapping> substitutions(new ProbabilisticSubstitutionMapping(*ppt, nbTypes, rootPatternLinks, nbDistinctSites));

  // Compute the number of substitutions for each class and each branch in the tree

  // Get the DAG of probabilities of the edges
  ProbabilityDAG probaDAG(rltc.getForwardLikelihoodTree(0));

  if (verbose)
    ApplicationTools::displayTask("Compute counts", true);

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt = substitutions->allEdgesIterator();

  size_t nn = 0;
  for ( ; !brIt->end(); brIt->next())
  {
    if (verbose)
      ApplicationTools::displayGauge(nn++, nbNodes - 1);

    shared_ptr<PhyloBranchMapping> br = **brIt;

    // For each branch
    uint speciesId = substitutions->getEdgeIndex(br);

    if (edgeIds.size() > 0 && !VectorTools::contains(edgeIds, (int)speciesId))
      continue;

//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wdelete-non-abstract-non-virtual-dtor" //Disable warning coming from Eigen because of non-virtual destructor.
    vector<RowLik> substitutionsForCurrentNode(nbTypes);
//#pragma GCC diagnostic pop
    for (auto& sub : substitutionsForCurrentNode)
    {
      sub = RowLik::Zero((int)nbDistinctSites);
    }

    for (size_t ncl = 0; ncl < nbClasses; ncl++)
    {
      processTree = rltc.getTreeNode(ncl);
      double pr = sp.getProbabilityForModel(ncl);

      vector<RowLik> substitutionsForCurrentClass(nbTypes);
      for (auto& sub:substitutionsForCurrentClass)
      {
        sub = RowLik::Zero((int)nbDistinctSites);
      }

      const auto& dagIndexes = rltc.getEdgesIds(speciesId, ncl);

      Eigen::MatrixXd npxy;

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

        auto subCount = mModCount[sm.get()];

        const auto& likelihoodsTopEdge = rltc.getBackwardLikelihoodsAtEdgeForClass(id, ncl)->targetValue();

        auto sonid = rltc.getForwardLikelihoodTree(ncl)->getSon(id);
        auto fatid = rltc.getForwardLikelihoodTree(ncl)->getFatherOfEdge(id);

        const auto& likelihoodsBotEdge = rltc.getForwardLikelihoodsAtNodeForClass(sonid, ncl)->targetValue();

        const Eigen::MatrixXd& pxy = edge->getTransitionMatrix()->targetValue();

        const auto& likelihoodsFather = rltc.getLikelihoodsAtNodeForClass(fatid, ncl)->targetValue();

        for (size_t t = 0; t < nbTypes; ++t)
        {
          // compute all nxy * pxy first:

          subCount->storeAllNumbersOfSubstitutions(edge->getBrLen()->getValue(), t + 1, npxy);

          npxy.array() *= pxy.array();

          // Now loop over sites:

          auto counts = npxy * likelihoodsBotEdge;

          auto bb = (cwise(likelihoodsTopEdge) * cwise(counts)).colwise().sum();

          Eigen::VectorXd ff(likelihoodsBotEdge.cols());
          switch(unresolvedOption){
          case SubstitutionMappingTools::UNRESOLVED_ZERO:
          case SubstitutionMappingTools::UNRESOLVED_AVERAGE:
            
            // Nullify counts where sum likelihoods > 1 : ie unknown
            for (auto i=0;i<ff.size();i++)
            {
              const DataLik s=likelihoodsBotEdge.col(i).sum();
              if (s>=2.)
                ff[i]=(unresolvedOption==SubstitutionMappingTools::UNRESOLVED_ZERO)?0.:1./convert(s);
              else
                ff[i]=1;
            }

            bb *= ff.array();
          default:
            ;
          }

          // Normalizes by likelihood on this node
          auto cc = bb / cwise(likelihoodsFather);

          // adds, with branch ponderation  ( * edge / edge * father) probs
          cwise(substitutionsForCurrentClass[t]) += cc * probaDAG.getProbaAtNode(fatid);
        }
      }

      // sum for all rate classes, with class ponderation
      for (size_t t = 0; t < nbTypes; ++t)
      {
        substitutionsForCurrentNode[t] += substitutionsForCurrentClass[t] * pr;
      }
    }

    // Now we just have to copy the substitutions into the result vector:

    for (size_t i = 0; i < nbDistinctSites; ++i)
    {
      for (size_t t = 0; t < nbTypes; ++t)
      {
        double x = convert(substitutionsForCurrentNode[t](Eigen::Index(i)));
        if (std::isnan(x) || std::isinf(x))
        {
          if (verbose)
            ApplicationTools::displayWarning("On branch " + TextTools::toString(speciesId) + ", site index " + TextTools::toString(i) + ", and type " + TextTools::toString(t) + ", counts could not be computed.");
          (*br)(i, t) = 0;
        }
        else
        {
          if (threshold >= 0 && x > threshold)
          {
            if (verbose)
              ApplicationTools::displayWarning("On branch " + TextTools::toString(speciesId) + ", site index" + TextTools::toString(i) + ", and type " + TextTools::toString(t) + " count has been ignored because it is presumably saturated.");
            (*br)(i, t) = 0;
          }
          else
            (*br)(i, t) = x;
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

  return substitutions;
}

/**************************************************************************************************/

unique_ptr<ProbabilisticSubstitutionMapping> SubstitutionMappingTools::computeNormalizations(
  LikelihoodCalculationSingleProcess& rltc,
  const vector<uint>& edgeIds,
  std::shared_ptr<const BranchedModelSet> nullModels,
  std::shared_ptr<const SubstitutionRegisterInterface> reg,
  std::shared_ptr<const AlphabetIndex2> distances,
  short unresolvedOption,
  bool verbose)
{
  // Preamble:
  if (!rltc.isInitialized())
    throw Exception("SubstitutionMappingTools::computeNormalizations(). Likelihood object is not initialized.");

  auto& sp = rltc.substitutionProcess();

  if (edgeIds.size() == 0)
    return make_unique<ProbabilisticSubstitutionMapping>(
	sp.parametrizablePhyloTree(),
	reg->getNumberOfSubstitutionTypes(),
	rltc.getRootArrayPositions(),
	rltc.getNumberOfDistinctSites());

  auto processTree = rltc.getTreeNode(0);


  /* First, set substitution rewards */

  // Map from models to type-vector of rewards
  std::map<const SubstitutionModelInterface*, std::vector<std::shared_ptr<DecompositionReward> > > mModRewards;

  const auto& stateMap = sp.stateMap();

  size_t nbTypes = reg->getNumberOfSubstitutionTypes();

  size_t nbStates = stateMap.getNumberOfModelStates();
  vector<int> supportedStates = stateMap.getAlphabetStates();

  vector<shared_ptr<UserAlphabetIndex1>> vusai(nbTypes);
  for (auto& usai : vusai) {  
    usai.reset(new UserAlphabetIndex1(stateMap.getAlphabet()));
  }

  for (auto speciesId : edgeIds)
  {
    const auto& dagIndexes = rltc.getEdgesIds(speciesId, 0);

    // look for matching null model
    auto nullmodel = nullModels->getModelForBranch(speciesId);

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

        // Sets vector of reward counts (depend on nullmodel)
        if (mModRewards.find(sm.get()) == mModRewards.end())
        {
          // Look for matching substitution nullmodel
          shared_ptr<const SubstitutionModelInterface> nullsm = nullptr;

          if (nMod == 0)
          {
            nullsm = dynamic_pointer_cast<const SubstitutionModelInterface>(nullmodel);

            if (!nullsm)
              throw Exception("SubstitutionMappingTools::computeNormalizations : SubstitutionVectors possible only for SubstitutionModels, not in branch " + TextTools::toString(speciesId) + "for null model " + nullmodel->getName());
          }
          else
          {
            size_t nmod = nMod->targetValue();

            auto nullttm = dynamic_pointer_cast<const MixedTransitionModelInterface>(nullmodel);
            if (!nullttm)
              throw Exception("SubstitutionMappingTools::computeNormalizations : Expecting Mixed model in branch " + TextTools::toString(speciesId) + " for null model " + nullttm->getName());

            nullsm = dynamic_pointer_cast<const SubstitutionModelInterface>(nullttm->getNModel(nmod));

            if (!nullsm)
              throw Exception("SubstitutionMappingTools::computeNormalizations : Expecting Substitution model for submodel " + TextTools::toString(nmod) + " of null mixed model " + nullmodel->getName() + " in branch " + TextTools::toString(speciesId));
          }

          for (auto& usai : vusai)
          {
            for (size_t i = 0; i < nbStates; ++i)
            {
              usai->setIndex(supportedStates[i], 0);
            }
          }

          for (size_t i = 0; i < nbStates; ++i)
          {
            for (size_t j = 0; j < nbStates; ++j)
            {
              if (i != j)
              {
                size_t nbt = reg->getType(i, j);
                if (nbt != 0)
                  vusai[nbt - 1]->setIndex(supportedStates[i], vusai[nbt - 1]->getIndex(supportedStates[i]) + nullsm->Qij(i, j) * (distances ? distances->getIndex(supportedStates[i], supportedStates[j]) : 1));
              }
            }
          }

          mModRewards[sm.get()] = vector<std::shared_ptr<DecompositionReward>>(nbTypes);
          auto& mMdodsm = mModRewards[sm.get()];
          for (size_t t = 0; t < nbTypes; ++t)
          {
            mMdodsm[t] = make_shared<DecompositionReward>(nullsm, vusai[t]);
          }
        }
      }
    }
  }

  //////////////////////////////////////////////////////
  //// Now the computation

  // A few variables we'll need:

  size_t nbDistinctSites = rltc.getNumberOfDistinctSites();
  size_t nbNodes         = edgeIds.size();
  size_t nbClasses       = sp.getNumberOfClasses();

  const auto& rootPatternLinks = rltc.getRootArrayPositions();

  // We create a Mapping objects

  unique_ptr<ProbabilisticSubstitutionMapping> normalizations(new ProbabilisticSubstitutionMapping(*sp.getParametrizablePhyloTree(), nbTypes, rootPatternLinks, nbDistinctSites));

  // Compute the reward for each class and each branch in the tree:

  // Get the DAG of probabilities of the edges
  ProbabilityDAG probaDAG(rltc.getForwardLikelihoodTree(0));

  if (verbose)
    ApplicationTools::displayTask("Compute rewards", true);

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt = normalizations->allEdgesIterator();

  Eigen::MatrixXd rpxy;

  size_t nn = 0;
  for ( ; !brIt->end(); brIt->next())
  {
    if (verbose)
      ApplicationTools::displayGauge(nn++, nbNodes - 1);

    shared_ptr<PhyloBranchMapping> br = **brIt;

    // For each branch
    uint speciesId = normalizations->getEdgeIndex(br);

    if (edgeIds.size() > 0 && !VectorTools::contains(edgeIds, (int)speciesId))
      continue;

    vector<RowLik> rewardsForCurrentNode(nbTypes);
    for (auto& sub:rewardsForCurrentNode)
    {
      sub = RowLik::Zero((int)nbDistinctSites);
    }

    for (size_t ncl = 0; ncl < nbClasses; ncl++)
    {
      processTree = rltc.getTreeNode(ncl);
      double pr = sp.getProbabilityForModel(ncl);

      vector<RowLik> rewardsForCurrentClass(nbTypes);
      for (auto& sub:rewardsForCurrentClass)
      {
        sub = RowLik::Zero((int)nbDistinctSites);
      }

      const auto& dagIndexes = rltc.getEdgesIds(speciesId, ncl);

      // Sum on all dag edges for this speciesId
      for (auto id:dagIndexes)
      {
        auto edge = processTree->getEdge(id);

        auto tm = dynamic_pointer_cast<const TransitionModelInterface>(edge->getModel()->targetValue());

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

        // Rewards for this model
        auto subReward = mModRewards[sm.get()];

        const auto& likelihoodsTopEdge = rltc.getBackwardLikelihoodsAtEdgeForClass(id, ncl)->targetValue();

        auto sonid = rltc.getForwardLikelihoodTree(ncl)->getSon(id);
        auto fatid = rltc.getForwardLikelihoodTree(ncl)->getFatherOfEdge(id);

        const auto& likelihoodsBotEdge = rltc.getForwardLikelihoodsAtNodeForClass(sonid, ncl)->targetValue();

        const Eigen::MatrixXd& pxy = edge->getTransitionMatrix()->accessValueConst();

        const auto& likelihoodsFather = rltc.getLikelihoodsAtNodeForClass(fatid, ncl)->targetValue();

        for (size_t t = 0; t < nbTypes; ++t)
        {
          // compute all rxy * pxy first:
          subReward[t]->storeAllRewards(edge->getBrLen()->getValue(), rpxy);

          rpxy.array() *= pxy.array();

          // Now loop over sites:

          auto rew = rpxy * likelihoodsBotEdge;

          auto bb = (cwise(likelihoodsTopEdge) * cwise(rew)).colwise().sum();

          // Nullify counts where sum likelihoods > 1 : ie unknown
          Eigen::VectorXd ff(likelihoodsBotEdge.cols());

          switch(unresolvedOption){
          case SubstitutionMappingTools::UNRESOLVED_ZERO:
          case SubstitutionMappingTools::UNRESOLVED_AVERAGE:
            
            // Nullify counts where sum likelihoods > 1 : ie unknown
            for (auto i=0;i<ff.size();i++)
            {
              const DataLik s=likelihoodsBotEdge.col(i).sum();
              if (s>=2.)
                ff[i]=(unresolvedOption==SubstitutionMappingTools::UNRESOLVED_ZERO)?0.:1./convert(s);
              else
                ff[i]=1;
            }
                    
            bb *= ff.array();
          default:
            ;
          }

          // Normalizes by likelihood on this node
          auto cc = bb / cwise(likelihoodsFather);

          // adds, with branch ponderation  ( * edge / edge * father) probs
          cwise(rewardsForCurrentClass[t]) += cc  * probaDAG.getProbaAtNode(fatid);
        }
      }

      for (size_t t = 0; t < nbTypes; ++t)
      {
        rewardsForCurrentNode[t] += rewardsForCurrentClass[t] * pr;
      }
    }

    // Now we just have to copy the substitutions into the result vector:
    for (size_t i = 0; i < nbDistinctSites; ++i)
    {
      for (size_t t = 0; t < nbTypes; ++t)
      {
        (*br)(i, t) = convert(rewardsForCurrentNode[t](Eigen::Index(i)));
      }
    }
  } // end of loop on branches

  if (verbose)
  {
    if (ApplicationTools::message)
      *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }

  return normalizations;
}


/************************************************************/

unique_ptr<ProbabilisticSubstitutionMapping> SubstitutionMappingTools::computeNormalizedCounts(
  LikelihoodCalculationSingleProcess& rltc,
  const vector<uint>& edgeIds,
  shared_ptr<const BranchedModelSet> nullModels,
  shared_ptr<const SubstitutionRegisterInterface> reg,
  shared_ptr<const AlphabetIndex2> weights,
  shared_ptr<const AlphabetIndex2> distances,
  bool perTimeUnit,
  uint siteSize,
  short unresolvedOption,
  double threshold,
  bool verbose)
{
    unique_ptr<ProbabilisticSubstitutionMapping> counts(computeCounts(rltc, edgeIds, reg, weights, distances, unresolvedOption, threshold, verbose));
    unique_ptr<ProbabilisticSubstitutionMapping> factors(computeNormalizations(rltc, edgeIds, nullModels, reg, distances, unresolvedOption, verbose));
  return computeNormalizedCounts(*counts, *factors, edgeIds, perTimeUnit, siteSize);
}

/************************************************************/

unique_ptr<ProbabilisticSubstitutionMapping> SubstitutionMappingTools::computeNormalizedCounts(
  const ProbabilisticSubstitutionMapping& counts,
  const ProbabilisticSubstitutionMapping& factors,
  const vector<uint>& edgeIds,
  bool perTimeUnit,
  uint siteSize)
{
  unique_ptr<ProbabilisticSubstitutionMapping> normCounts(counts.clone());

  size_t nbTypes = counts.getNumberOfSubstitutionTypes();
  size_t nbDistinctSites = counts.getNumberOfDistinctSites();

  // Iterate on branches

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt = normCounts->allEdgesIterator();

  for ( ; !brIt->end(); brIt->next())
  {
    shared_ptr<PhyloBranchMapping> brNormCount = **brIt;

    VVdouble& brnCou = brNormCount->getCounts();

    // For each branch
    uint edid = normCounts->getEdgeIndex(brNormCount);

    if (edgeIds.size() > 0 && !VectorTools::contains(edgeIds, (int)edid))
    {
      for (auto& brk : brnCou)
      {
        VectorTools::fill(brk, 0.);
      }
      continue;
    }

    shared_ptr<PhyloBranchMapping> brFactor = factors.getEdge(edid);
    shared_ptr<PhyloBranchMapping> brCount = counts.getEdge(edid);

    const VVdouble& cou = brCount->getCounts();
    const VVdouble& fac = brFactor->getCounts();


    // if not per time, multiply by the lengths of the branches of
    // the input tree

    double slg = (!perTimeUnit ? brCount->getLength() : 1) / siteSize;

    for (size_t k = 0; k < nbDistinctSites; ++k)
    {
      Vdouble& ncou_k = brnCou[k];

      const Vdouble& cou_k = cou[k];
      const Vdouble& fac_k = fac[k];

      for (size_t t = 0; t < nbTypes; ++t)
      {
        ncou_k[t] = (fac_k[t] != 0 ? cou_k[t] / fac_k[t] * slg : 0);
      }
    }
  }

  return normCounts;
}

/*******************************************************************************/
/* Get trees of counts */
/*******************************************************************************/

unique_ptr<PhyloTree> SubstitutionMappingTools::getTreeForType(
    const ProbabilisticSubstitutionMapping& counts,
    size_t type)
{
  auto pt = make_unique<PhyloTree>(counts);
  size_t nbSites = counts.getNumberOfSites();

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt = counts.allEdgesIterator();

  for ( ; !brIt->end(); brIt->next())
  {
    shared_ptr<PhyloBranchMapping> brm = (**brIt);
    double x = 0;

    for (size_t i = 0; i < nbSites; ++i)
    {
      x += brm->getSiteTypeCount(counts.getSiteIndex(i), type);
    }

    pt->getEdge(counts.getEdgeIndex(brm))->setLength(x);
  }

  return pt;
}

/*******************************************************************************/

unique_ptr<PhyloTree> SubstitutionMappingTools::getTreeForType(
    const ProbabilisticSubstitutionMapping& counts,
    const ProbabilisticSubstitutionMapping& factors,
    size_t type)
{
  auto pt = make_unique<PhyloTree>(counts);
  size_t nbSites = counts.getNumberOfSites();

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt = counts.allEdgesIterator();

  for ( ; !brIt->end(); brIt->next())
  {
    shared_ptr<PhyloBranchMapping> brm = (**brIt);
    shared_ptr<PhyloBranchMapping> brf = factors.getEdge(counts.getEdgeIndex(brm));

    double x = 0, f = 0;

    for (size_t i = 0; i < nbSites; ++i)
    {
      x += brm->getSiteTypeCount(counts.getSiteIndex(i), type);
      f += brf->getSiteTypeCount(counts.getSiteIndex(i), type);
    }

    pt->getEdge(counts.getEdgeIndex(brm))->setLength(x / f);
  }

  return pt;
}


/********************************************************************************/
/*  Get vectors of counts  */
/********************************************************************************/

VVVdouble SubstitutionMappingTools::getCountsPerSitePerBranchPerType(
  const ProbabilisticSubstitutionMapping& counts,
  const vector<uint>& ids)
{
  const Vuint idc(ids.size() == 0 ? counts.getAllEdgesIndexes() : ids);
  size_t nbBr = idc.size();

  size_t nbSites = counts.getNumberOfSites();

  size_t nbTypes = counts.getNumberOfSubstitutionTypes();

  VVVdouble result;
  VectorTools::resize3(result, nbSites, nbBr, nbTypes);

  for (size_t j = 0; j < nbSites; ++j)
  {
    VVdouble& resS = result[j];
    size_t siteIndex = counts.getSiteIndex(j);

    for (size_t k = 0; k < nbBr; ++k)
    {
      Vdouble& resSB = resS[k];

      for (size_t i = 0; i < nbTypes; ++i)
      {
        resSB[i] = counts(idc[k], siteIndex, i);
      }
    }
  }

  return result;
}


/**************************************************************************************************/

Vdouble SubstitutionMappingTools::getCountsForSitePerBranch(
  const ProbabilisticSubstitutionMapping& counts,
  size_t site)
{
  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt = counts.allEdgesIterator();
  size_t siteIndex = counts.getSiteIndex(site);

  Vdouble v(counts.getNumberOfBranches(), 0);
  for ( ; !brIt->end(); brIt->next())
  {
    v[counts.getEdgeIndex(**brIt)] = VectorTools::sum((***brIt).getSiteCount(siteIndex));
  }

  return v;
}

Vdouble SubstitutionMappingTools::getCountsForSitePerBranch(
  const ProbabilisticSubstitutionMapping& counts,
  const ProbabilisticSubstitutionMapping& factors,
  size_t site)
{
  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt = counts.allEdgesIterator();
  size_t siteIndex = counts.getSiteIndex(site);

  Vdouble v(counts.getNumberOfBranches(), 0);
  for ( ; !brIt->end(); brIt->next())
  {
    shared_ptr<PhyloBranchMapping> brm = (**brIt);

    uint edid = counts.getEdgeIndex(brm);
    shared_ptr<PhyloBranchMapping> brf = factors.getEdge(edid);


    v[edid] = VectorTools::sum(brm->getSiteCount(siteIndex)) / VectorTools::sum(brf->getSiteCount(siteIndex));
  }

  return v;
}

/**************************************************************************************************/

VVdouble SubstitutionMappingTools::getCountsPerSitePerBranch(
  const ProbabilisticSubstitutionMapping& counts,
  const vector<uint>& ids)
{
  const Vuint idc(ids.size() == 0 ? counts.getAllEdgesIndexes() : ids);
  size_t nbBr = idc.size();

  size_t nbSites = counts.getNumberOfSites();

  VVdouble result;
  VectorTools::resize2(result, nbSites, nbBr);

  for (size_t k = 0; k < nbSites; ++k)
  {
    vector<double> countsf(SubstitutionMappingTools::getCountsForSitePerBranch(counts, k));
    Vdouble* resS = &result[k];

    for (size_t i = 0; i < nbBr; ++i)
    {
      (*resS)[i] = countsf[idc[i]];
    }
  }
  return result;
}

VVdouble SubstitutionMappingTools::getCountsPerSitePerBranch(
  const ProbabilisticSubstitutionMapping& counts,
  const ProbabilisticSubstitutionMapping& factors,
  const vector<uint>& ids)
{
  const Vuint idc(ids.size() == 0 ? counts.getAllEdgesIndexes() : ids);
  size_t nbBr = idc.size();

  size_t nbSites = counts.getNumberOfSites();

  VVdouble result;
  VectorTools::resize2(result, nbSites, nbBr);

  for (size_t k = 0; k < nbSites; ++k)
  {
    vector<double> countsf(SubstitutionMappingTools::getCountsForSitePerBranch(counts, factors, k));
    Vdouble* resS = &result[k];

    for (size_t i = 0; i < nbBr; ++i)
    {
      (*resS)[i] = countsf[idc[i]];
    }
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
  shared_ptr<PhyloBranchMapping> br = counts.getEdge(branchId);

  for (size_t i = 0; i < nbSites; ++i)
  {
    for (size_t t = 0; t < nbTypes; ++t)
    {
      v[t] += br->getSiteTypeCount(counts.getSiteIndex(i), t);
    }
  }

  return v;
}

/**************************************************************************************************/

VVdouble SubstitutionMappingTools::getCountsPerBranchPerType(
  const ProbabilisticSubstitutionMapping& counts,
  const vector<uint>& ids)
{
  VVdouble result;
  const Vuint idc(ids.size() == 0 ? counts.getAllEdgesIndexes() : ids);

  for (auto id : idc)
  {
    result.push_back(SubstitutionMappingTools::getCountsForBranchPerType(counts, id));
  }

  return result;
}

/**************************************************************************************************/

VVdouble SubstitutionMappingTools::getCountsPerTypePerBranch(
  const ProbabilisticSubstitutionMapping& counts,
  const vector<uint>& ids)
{
  const Vuint idc(ids.size() == 0 ? counts.getAllEdgesIndexes() : ids);
  size_t nbBr = idc.size();

  size_t nbTypes = counts.getNumberOfSubstitutionTypes();

  VVdouble result;
  VectorTools::resize2(result, nbTypes, nbBr);

  for (size_t i = 0; i < idc.size(); i++)
  {
    Vdouble cou = SubstitutionMappingTools::getCountsForBranchPerType(counts, idc[i]);
    for (size_t nbt = 0; nbt < nbTypes; nbt++)
    {
      result[nbt][i] = cou[nbt];
    }
  }

  return result;
}


/************************************************************/


VVdouble SubstitutionMappingTools::computeCountsPerTypePerBranch(
  LikelihoodCalculationSingleProcess& rltc,
  const vector<uint>& ids,
  shared_ptr<const SubstitutionRegisterInterface> reg,
  std::shared_ptr<const AlphabetIndex2> weights,
  std::shared_ptr<const AlphabetIndex2> distances,
  short unresolvedOption,
  double threshold,
  bool verbose)
{
  auto psm = computeCounts(rltc, ids, reg, weights, distances, unresolvedOption, threshold, verbose);

  VVdouble result = getCountsPerTypePerBranch(*psm, ids);

  auto creg = dynamic_pointer_cast<const CategorySubstitutionRegister>(reg);

  if (creg && !creg->isStationary())
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
  const Vuint idc(ids.size() == 0 ? counts.getAllEdgesIndexes() : ids);

  size_t siteIndex = counts.getSiteIndex(site);

  size_t nbTypes    = counts.getNumberOfSubstitutionTypes();
  Vdouble v(nbTypes, 0);

  for (auto id : idc)
  {
    shared_ptr<PhyloBranchMapping> br = counts.getEdge(id);
    for (size_t t = 0; t < nbTypes; ++t)
    {
      v[t] += (*br)(siteIndex, t);
    }
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

  Vdouble v(nbTypes, 0);
  Vdouble n(nbTypes, 0);

  double lg = (!perTimeUnit ? 0 : 1);
  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt = counts.allEdgesIterator();
  for ( ; !brIt->end(); brIt->next())
  {
    for (size_t t = 0; t < nbTypes; ++t)
    {
      v[t] += (***brIt)(siteIndex, t);
    }
    if (!perTimeUnit)
      lg += (**brIt)->getLength();
  }

  double slg = lg / siteSize;

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> nrIt = factors.allEdgesIterator();
  for ( ; !nrIt->end(); nrIt->next())
  {
    for (size_t t = 0; t < nbTypes; ++t)
    {
      n[t] += (***nrIt)(siteIndex, t);
    }
  }

  for (size_t t = 0; t < nbTypes; ++t)
  {
    v[t] = v[t] / n[t] * slg;
  }

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
  const Vuint idc(ids.size() == 0 ? counts.getAllEdgesIndexes() : ids);

  size_t siteIndex = counts.getSiteIndex(site);
  size_t nbTypes   = counts.getNumberOfSubstitutionTypes();

  Vdouble v(nbTypes, 0);
  Vdouble n(nbTypes, 0);

  double lg = (!perTimeUnit ? 0 : 1);
  for (auto id : idc)
  {
    shared_ptr<PhyloBranchMapping> br = counts.getEdge(id);
    for (size_t t = 0; t < nbTypes; ++t)
    {
      v[t] += (*br)(siteIndex, t);
    }
    if (!perTimeUnit)
      lg += br->getLength();
  }

  double slg = lg / siteSize;

  for (auto id : idc)
  {
    shared_ptr<PhyloBranchMapping> br = factors.getEdge(id);
    for (size_t t = 0; t < nbTypes; ++t)
    {
      n[t] += (*br)(siteIndex, t);
    }
  }

  for (size_t t = 0; t < nbTypes; ++t)
  {
    v[t] = v[t] / n[t] * slg;
  }

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
  {
    result[k] = SubstitutionMappingTools::getCountsForSitePerType(counts, factors, k, perTimeUnit, siteSize);
  }

  return result;
}

/**************************************************************************************************/

VVdouble SubstitutionMappingTools::getCountsPerSitePerType(const ProbabilisticSubstitutionMapping& counts, const vector<uint>& ids)
{
  const Vuint idc(ids.size() == 0 ? counts.getAllEdgesIndexes() : ids);

  size_t nbSites = counts.getNumberOfSites();
  size_t nbTypes = counts.getNumberOfSubstitutionTypes();
  VVdouble result;
  VectorTools::resize2(result, nbSites, nbTypes);

  for (size_t k = 0; k < nbSites; ++k)
  {
    result[k] = SubstitutionMappingTools::getCountsForSitePerType(counts, k, idc);
  }

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
  const Vuint idc(ids.size() == 0 ? counts.getAllEdgesIndexes() : ids);

  size_t nbSites = counts.getNumberOfSites();
  size_t nbTypes = counts.getNumberOfSubstitutionTypes();
  VVdouble result;
  VectorTools::resize2(result, nbSites, nbTypes);

  for (size_t k = 0; k < nbSites; ++k)
  {
    result[k] = SubstitutionMappingTools::getCountsForSitePerType(counts, factors, k, idc, perTimeUnit, siteSize);
  }

  return result;
}

/**************************************************************************************************/

double SubstitutionMappingTools::getNormForSite(const ProbabilisticSubstitutionMapping& counts, size_t site)
{
  size_t siteIndex = counts.getSiteIndex(site);
  double sumSquare = 0;

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt = counts.allEdgesIterator();

  size_t nbTypes    = counts.getNumberOfSubstitutionTypes();
  for ( ; !brIt->end(); brIt->next())
  {
    double sum = 0;
    for (size_t t = 0; t < nbTypes; ++t)
    {
      sum += (***brIt)(siteIndex, t);
    }

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
  const AlignmentDataInterface& sites,
  const VVdouble& counts)
{
  size_t nbSites = counts.size();
  if (nbSites == 0)
    return;
  size_t nbBr = counts[0].size();

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
    const Vdouble& countS = counts[k];
    file << sites.site(k).getCoordinate();
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
  const SubstitutionRegisterInterface& reg,
  const AlignmentDataInterface& sites,
  const VVdouble& counts)
{
  size_t nbSites = counts.size();
  if (nbSites == 0)
    return;
  size_t nbTypes = counts[0].size();

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
    file << sites.site(k).getCoordinate();
    const Vdouble& resS = counts[k];
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
  const SubstitutionRegisterInterface& reg,
  const AlignmentDataInterface& sites,
  const VVVdouble& counts)
{
  size_t nbSites = counts.size();
  if (nbSites == 0)
    return;
  size_t nbBr = counts[0].size();
  if (nbBr == 0)
    return;
  size_t nbTypes = counts[0][0].size();

  ofstream file;

  for (size_t i = 0; i < nbTypes; ++i)
  {
    string name = reg.getTypeName(i + 1);
    if (name == "")
      name = TextTools::toString(i + 1);

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
      const VVdouble& resS = counts[j];

      file << sites.site(j).getCoordinate();

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
  const AlignmentDataInterface& sites,
  size_t type,
  ostream& out)
{
  if (!out)
    throw IOException("SubstitutionMappingTools::writeToFile. Can't write to stream.");
  out << "Branches";
  out << "\tMean";
  for (size_t i = 0; i < substitutions.getNumberOfSites(); ++i)
  {
    out << "\tSite" << sites.site(i).getCoordinate();
  }
  out << endl;

  unique_ptr<ProbabilisticSubstitutionMapping::mapTree::EdgeIterator> brIt = substitutions.allEdgesIterator();

  for ( ; !brIt->end(); brIt->next())
  {
    const shared_ptr<PhyloBranchMapping> br = **brIt;

    out << substitutions.getEdgeIndex(br) << "\t" << br->getLength();

    for (size_t i = 0; i < substitutions.getNumberOfSites(); i++)
    {
      out << "\t" << (*br)(substitutions.getSiteIndex(i), type);
    }

    out << endl;
  }
}

/**************************************************************************************************/

void SubstitutionMappingTools::readFromStream(istream& in, ProbabilisticSubstitutionMapping& substitutions, size_t type)
{
  try
  {
    auto data = DataTable::read(in, "\t", true, -1);
    vector<string> ids = data->getColumn(0);
    data->deleteColumn(0); // Remove ids
    data->deleteColumn(0); // Remove means
    // Now parse the table:
    size_t nbSites = data->getNumberOfColumns();
    substitutions.setNumberOfSites(nbSites);

    auto brIt = substitutions.allEdgesIterator();

    for ( ; !brIt->end(); brIt->next())
    {
      const shared_ptr<PhyloBranchMapping> br = **brIt;
      uint brid = substitutions.getEdgeIndex(br);
      for (size_t j = 0; j < nbSites; j++)
      {
        (*br)(j, type) = TextTools::toDouble((*data)(brid, j));
      }
    }

    // Parse the header:
    for (size_t i = 0; i < nbSites; ++i)
    {
      string siteTxt = data->getColumnName(i);
      int site = 0;
      if (siteTxt.substr(0, 4) == "Site")
        site = TextTools::to<int>(siteTxt.substr(4));
      else
        site = TextTools::to<int>(siteTxt);
      substitutions.setSitePosition(i, site);
    }
  }
  catch (Exception& e)
  {
    throw IOException(string("Bad input file. ") + e.what());
  }
}

