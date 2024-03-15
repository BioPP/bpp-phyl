// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SubstitutionMappingTools.h"
#include "../../Mapping/UniformizationSubstitutionCount.h"
#include "../../Mapping/DecompositionReward.h"
#include "ProbabilisticRewardMapping.h"
#include "RewardMappingTools.h"
#include "../Likelihood/DRTreeLikelihoodTools.h"
#include "../Likelihood/MarginalAncestralStateReconstruction.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>

using namespace bpp;

// From the STL:
#include <iomanip>

using namespace std;

/******************************************************************************/

unique_ptr<LegacyProbabilisticSubstitutionMapping> LegacySubstitutionMappingTools::computeSubstitutionVectors(
  std::shared_ptr<const DRTreeLikelihoodInterface> drtl,
  const vector<int>& nodeIds,
  std::shared_ptr<SubstitutionCountInterface> substitutionCount,
  bool verbose)
{
  // Preamble:
  if (!drtl->isInitialized())
    throw Exception("SubstitutionMappingTools::computeSubstitutionVectors(). Likelihood object is not initialized.");

  // A few variables we'll need:

  const TreeTemplate<Node> tree(drtl->tree());
  auto& sequences = drtl->data();
  auto& rDist = drtl->rateDistribution();

  size_t nbSites         = sequences.getNumberOfSites();
  size_t nbDistinctSites = drtl->likelihoodData().getNumberOfDistinctSites();
  size_t nbStates        = sequences.alphabet().getSize();
  size_t nbClasses       = rDist.getNumberOfCategories();
  size_t nbTypes         = substitutionCount->getNumberOfSubstitutionTypes();
  vector<const Node*> nodes    = tree.getNodes();
  const auto& rootPatternLinks = drtl->likelihoodData().getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  size_t nbNodes         = nodes.size();

  // We create a new ProbabilisticSubstitutionMapping object:
  auto substitutions = make_unique<LegacyProbabilisticSubstitutionMapping>(tree, substitutionCount, nbSites);

  // Store likelihood for each rate for each site:
  VVVdouble lik;
  drtl->computeLikelihoodAtNode(tree.getRootId(), lik);
  Vdouble Lr(nbDistinctSites, 0);
  Vdouble rcProbs = rDist.getProbabilities();
  Vdouble rcRates = rDist.getCategories();
  for (size_t i = 0; i < nbDistinctSites; i++)
  {
    VVdouble* lik_i = &lik[i];
    for (size_t c = 0; c < nbClasses; c++)
    {
      Vdouble* lik_i_c = &(*lik_i)[c];
      double rc = rDist.getProbability(c);
      for (size_t s = 0; s < nbStates; s++)
      {
        Lr[i] += (*lik_i_c)[s] * rc;
      }
    }
  }

  // Compute the number of substitutions for each class and each branch in the tree:
  if (verbose)
    ApplicationTools::displayTask("Compute joint node-pairs likelihood", true);

  for (size_t l = 0; l < nbNodes; ++l)
  {
    // For each node,
    const Node* currentNode = nodes[l];
    if (nodeIds.size() > 0 && !VectorTools::contains(nodeIds, currentNode->getId()))
      continue;

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
    for (size_t i = 0; i < nbDistinctSites; i++)
    {
      VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
      likelihoodsFatherConstantPart_i->resize(nbClasses);
      for (size_t c = 0; c < nbClasses; c++)
      {
        Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
        likelihoodsFatherConstantPart_i_c->resize(nbStates);
        double rc = rDist.getProbability(c);
        for (size_t s = 0; s < nbStates; s++)
        {
          // (* likelihoodsFatherConstantPart_i_c)[s] = rc * model->freq(s);
          // freq is already accounted in the array
          (*likelihoodsFatherConstantPart_i_c)[s] = rc;
        }
      }
    }

    // First, what will remain constant:
    size_t nbSons =  father->getNumberOfSons();
    for (size_t n = 0; n < nbSons; n++)
    {
      const Node* currentSon = father->getSon(n);
      if (currentSon->getId() != currentNode->getId())
      {
        const VVVdouble& likelihoodsFather_son = drtl->likelihoodData().getLikelihoodArray(father->getId(), currentSon->getId());

        // Now iterate over all site partitions:
        unique_ptr<TreeLikelihoodInterface::ConstBranchModelIterator> mit(drtl->getNewBranchModelIterator(currentSon->getId()));
        VVVdouble pxy;
        bool first;
        while (mit->hasNext())
        {
          auto bmd = mit->next();
          unique_ptr<TreeLikelihoodInterface::SiteIterator> sit(bmd->getNewSiteIterator());
          first = true;
          while (sit->hasNext())
          {
            size_t i = sit->next();
            // We retrieve the transition probabilities for this site partition:
            if (first)
            {
              pxy = drtl->getTransitionProbabilitiesPerRateClass(currentSon->getId(), i);
              first = false;
            }
            const VVdouble& likelihoodsFather_son_i = likelihoodsFather_son[i];
            VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
            for (size_t c = 0; c < nbClasses; c++)
            {
              const Vdouble& likelihoodsFather_son_i_c = likelihoodsFather_son_i[c];
              Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
              VVdouble& pxy_c = pxy[c];
              for (size_t x = 0; x < nbStates; x++)
              {
                Vdouble& pxy_c_x = pxy_c[x];
                double likelihood = 0.;
                for (size_t y = 0; y < nbStates; y++)
                {
                  likelihood += pxy_c_x[y] * likelihoodsFather_son_i_c[y];
                }
                likelihoodsFatherConstantPart_i_c[x] *= likelihood;
              }
            }
          }
        }
      }
    }
    if (father->hasFather())
    {
      const Node* currentSon = father->getFather();
      const VVVdouble& likelihoodsFather_son = drtl->likelihoodData().getLikelihoodArray(father->getId(), currentSon->getId());
      // Now iterate over all site partitions:
      unique_ptr<TreeLikelihoodInterface::ConstBranchModelIterator> mit(drtl->getNewBranchModelIterator(father->getId()));
      VVVdouble pxy;
      bool first;
      while (mit->hasNext())
      {
        auto bmd = mit->next();
        unique_ptr<TreeLikelihoodInterface::SiteIterator> sit(bmd->getNewSiteIterator());
        first = true;
        while (sit->hasNext())
        {
          size_t i = sit->next();
          // We retrieve the transition probabilities for this site partition:
          if (first)
          {
            pxy = drtl->getTransitionProbabilitiesPerRateClass(father->getId(), i);
            first = false;
          }
          const VVdouble& likelihoodsFather_son_i = likelihoodsFather_son[i];
          VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
          for (size_t c = 0; c < nbClasses; c++)
          {
            const Vdouble& likelihoodsFather_son_i_c = likelihoodsFather_son_i[c];
            Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
            VVdouble& pxy_c = pxy[c];
            for (size_t x = 0; x < nbStates; x++)
            {
              double likelihood = 0.;
              for (size_t y = 0; y < nbStates; y++)
              {
                Vdouble& pxy_c_x = pxy_c[y];
                likelihood += pxy_c_x[x] * likelihoodsFather_son_i_c[y];
              }
              likelihoodsFatherConstantPart_i_c[x] *= likelihood;
            }
          }
        }
      }
    }
    else
    {
      // Account for root frequencies:
      for (size_t i = 0; i < nbDistinctSites; i++)
      {
        vector<double> freqs = drtl->getRootFrequencies(i);
        VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
        for (size_t c = 0; c < nbClasses; c++)
        {
          Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
          for (size_t x = 0; x < nbStates; x++)
          {
            likelihoodsFatherConstantPart_i_c[x] *= freqs[x];
          }
        }
      }
    }


    // Then, we deal with the node of interest.
    // We first average upon 'y' to save computations, and then upon 'x'.
    // ('y' is the state at 'node' and 'x' the state at 'father'.)

    // Iterate over all site partitions:
    const VVVdouble& likelihoodsFather_node = drtl->likelihoodData().getLikelihoodArray(father->getId(), currentNode->getId());
    unique_ptr<TreeLikelihoodInterface::ConstBranchModelIterator> mit(drtl->getNewBranchModelIterator(currentNode->getId()));
    VVVdouble pxy;
    bool first;
    while (mit->hasNext())
    {
      auto bmd = mit->next();
      substitutionCount->setSubstitutionModel(bmd->getSubstitutionModel());
      // compute all nxy first:
      VVVVdouble nxy(nbClasses);
      for (size_t c = 0; c < nbClasses; ++c)
      {
        VVVdouble& nxy_c = nxy[c];
        double rc = rcRates[c];
        nxy_c.resize(nbTypes);
        for (size_t t = 0; t < nbTypes; ++t)
        {
          VVdouble& nxy_c_t = nxy_c[t];
          auto nijt = substitutionCount->getAllNumbersOfSubstitutions(d * rc, t + 1);
           
          nxy_c_t.resize(nbStates);
          for (size_t x = 0; x < nbStates; ++x)
          {
            Vdouble& nxy_c_t_x = nxy_c_t[x];
            nxy_c_t_x.resize(nbStates);
            for (size_t y = 0; y < nbStates; ++y)
            {
              nxy_c_t_x[y] = (*nijt)(x, y);
            }
          }
        }
      }

      // Now loop over sites:
      unique_ptr<TreeLikelihoodInterface::SiteIterator> sit(bmd->getNewSiteIterator());
      first = true;
      while (sit->hasNext())
      {
        size_t i = sit->next();
        // We retrieve the transition probabilities and substitution counts for this site partition:
        if (first)
        {
          pxy = drtl->getTransitionProbabilitiesPerRateClass(currentNode->getId(), i);
          first = false;
        }
        const VVdouble& likelihoodsFather_node_i = likelihoodsFather_node[i];
        VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
        for (size_t c = 0; c < nbClasses; ++c)
        {
          const Vdouble& likelihoodsFather_node_i_c = likelihoodsFather_node_i[c];
          Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
          const VVdouble& pxy_c = pxy[c];
          VVVdouble& nxy_c = nxy[c];
          for (size_t x = 0; x < nbStates; ++x)
          {
            double& likelihoodsFatherConstantPart_i_c_x = likelihoodsFatherConstantPart_i_c[x];
            const Vdouble& pxy_c_x = pxy_c[x];
            for (size_t y = 0; y < nbStates; ++y)
            {
              double likelihood_cxy = likelihoodsFatherConstantPart_i_c_x
                                      * pxy_c_x[y]
                                      * likelihoodsFather_node_i_c[y];

              for (size_t t = 0; t < nbTypes; ++t)
              {
                // Now the vector computation:
                substitutionsForCurrentNode[i][t] += likelihood_cxy * nxy_c[t][x][y];
                //                                   <------------>   <------------>
                // Posterior probability                   |                |
                // for site i and rate class c *           |                |
                // likelihood for this site----------------+                |
                //                                                          |
                // Substitution function for site i and rate class c--------+
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
        (*substitutions)(l, i, t) = substitutionsForCurrentNode[rootPatternLinks[i]][t] / Lr[rootPatternLinks[i]];
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

/******************************************************************************/

unique_ptr<LegacyProbabilisticSubstitutionMapping> LegacySubstitutionMappingTools::computeSubstitutionVectors(
  std::shared_ptr<const DRTreeLikelihoodInterface> drtl,
  const SubstitutionModelSet& modelSet,
  const vector<int>& nodeIds,
  std::shared_ptr<SubstitutionCountInterface> substitutionCount,
  bool verbose)
{
  // Preamble:
  if (!drtl->isInitialized())
    throw Exception("SubstitutionMappingTools::computeSubstitutionVectors(). Likelihood object is not initialized.");

  // A few variables we'll need:

  const TreeTemplate<Node> tree(drtl->tree());
  auto& sequences = drtl->data();
  auto& rDist = drtl->rateDistribution();

  size_t nbSites         = sequences.getNumberOfSites();
  size_t nbDistinctSites = drtl->likelihoodData().getNumberOfDistinctSites();
  size_t nbStates        = sequences.alphabet().getSize();
  size_t nbClasses       = rDist.getNumberOfCategories();
  size_t nbTypes         = substitutionCount->getNumberOfSubstitutionTypes();
  vector<const Node*> nodes    = tree.getNodes();
  auto& rootPatternLinks = drtl->likelihoodData().getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  size_t nbNodes         = nodes.size();

  // We create a new ProbabilisticSubstitutionMapping object:
  auto substitutions = make_unique<LegacyProbabilisticSubstitutionMapping>(tree, substitutionCount, nbSites);

  // Store likelihood for each rate for each site:
  VVVdouble lik;
  drtl->computeLikelihoodAtNode(tree.getRootId(), lik);
  Vdouble Lr(nbDistinctSites, 0);
  Vdouble rcProbs = rDist.getProbabilities();
  Vdouble rcRates = rDist.getCategories();
  for (size_t i = 0; i < nbDistinctSites; i++)
  {
    VVdouble* lik_i = &lik[i];
    for (size_t c = 0; c < nbClasses; c++)
    {
      Vdouble* lik_i_c = &(*lik_i)[c];
      double rc = rDist.getProbability(c);
      for (size_t s = 0; s < nbStates; s++)
      {
        Lr[i] += (*lik_i_c)[s] * rc;
      }
    }
  }

  // Compute the number of substitutions for each class and each branch in the tree:
  if (verbose)
    ApplicationTools::displayTask("Compute joint node-pairs likelihood", true);

  for (size_t l = 0; l < nbNodes; ++l)
  {
    // For each node,
    const Node* currentNode = nodes[l];
    if (nodeIds.size() > 0 && !VectorTools::contains(nodeIds, currentNode->getId()))
      continue;

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
    for (size_t i = 0; i < nbDistinctSites; i++)
    {
      VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
      likelihoodsFatherConstantPart_i.resize(nbClasses);
      for (size_t c = 0; c < nbClasses; c++)
      {
        Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
        likelihoodsFatherConstantPart_i_c.resize(nbStates);
        double rc = rDist.getProbability(c);
        for (size_t s = 0; s < nbStates; s++)
        {
          // (* likelihoodsFatherConstantPart_i_c)[s] = rc * model->freq(s);
          // freq is already accounted in the array
          likelihoodsFatherConstantPart_i_c[s] = rc;
        }
      }
    }

    // First, what will remain constant:
    size_t nbSons =  father->getNumberOfSons();
    for (size_t n = 0; n < nbSons; n++)
    {
      const Node* currentSon = father->getSon(n);
      if (currentSon->getId() != currentNode->getId())
      {
        const VVVdouble& likelihoodsFather_son = drtl->likelihoodData().getLikelihoodArray(father->getId(), currentSon->getId());

        // Now iterate over all site partitions:
        unique_ptr<TreeLikelihoodInterface::ConstBranchModelIterator> mit(drtl->getNewBranchModelIterator(currentSon->getId()));
        VVVdouble pxy;
        bool first;
        while (mit->hasNext())
        {
          auto bmd = mit->next();
          unique_ptr<TreeLikelihoodInterface::SiteIterator> sit(bmd->getNewSiteIterator());
          first = true;
          while (sit->hasNext())
          {
            size_t i = sit->next();
            // We retrieve the transition probabilities for this site partition:
            if (first)
            {
              pxy = drtl->getTransitionProbabilitiesPerRateClass(currentSon->getId(), i);
              first = false;
            }
            const VVdouble& likelihoodsFather_son_i = likelihoodsFather_son[i];
            VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
            for (size_t c = 0; c < nbClasses; c++)
            {
              const Vdouble& likelihoodsFather_son_i_c = likelihoodsFather_son_i[c];
              Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
              VVdouble& pxy_c = pxy[c];
              for (size_t x = 0; x < nbStates; x++)
              {
                Vdouble& pxy_c_x = pxy_c[x];
                double likelihood = 0.;
                for (size_t y = 0; y < nbStates; y++)
                {
                  likelihood += pxy_c_x[y] * likelihoodsFather_son_i_c[y];
                }
                likelihoodsFatherConstantPart_i_c[x] *= likelihood;
              }
            }
          }
        }
      }
    }
    if (father->hasFather())
    {
      const Node* currentSon = father->getFather();
      const VVVdouble& likelihoodsFather_son = drtl->likelihoodData().getLikelihoodArray(father->getId(), currentSon->getId());
      // Now iterate over all site partitions:
      unique_ptr<TreeLikelihoodInterface::ConstBranchModelIterator> mit(drtl->getNewBranchModelIterator(father->getId()));
      VVVdouble pxy;
      bool first;
      while (mit->hasNext())
      {
        auto bmd = mit->next();
        unique_ptr<TreeLikelihoodInterface::SiteIterator> sit(bmd->getNewSiteIterator());
        first = true;
        while (sit->hasNext())
        {
          size_t i = sit->next();
          // We retrieve the transition probabilities for this site partition:
          if (first)
          {
            pxy = drtl->getTransitionProbabilitiesPerRateClass(father->getId(), i);
            first = false;
          }
          const VVdouble& likelihoodsFather_son_i = likelihoodsFather_son[i];
          VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
          for (size_t c = 0; c < nbClasses; c++)
          {
            const Vdouble& likelihoodsFather_son_i_c = likelihoodsFather_son_i[c];
            Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
            VVdouble& pxy_c = pxy[c];
            for (size_t x = 0; x < nbStates; x++)
            {
              double likelihood = 0.;
              for (size_t y = 0; y < nbStates; y++)
              {
                Vdouble& pxy_c_x = pxy_c[y];
                likelihood += pxy_c_x[x] * likelihoodsFather_son_i_c[y];
              }
              likelihoodsFatherConstantPart_i_c[x] *= likelihood;
            }
          }
        }
      }
    }
    else
    {
      // Account for root frequencies:
      for (size_t i = 0; i < nbDistinctSites; i++)
      {
        vector<double> freqs = drtl->getRootFrequencies(i);
        VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
        for (size_t c = 0; c < nbClasses; c++)
        {
          Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
          for (size_t x = 0; x < nbStates; x++)
          {
            likelihoodsFatherConstantPart_i_c[x] *= freqs[x];
          }
        }
      }
    }


    // Then, we deal with the node of interest.
    // We first average upon 'y' to save computations, and then upon 'x'.
    // ('y' is the state at 'node' and 'x' the state at 'father'.)

    // Iterate over all site partitions:
    const VVVdouble& likelihoodsFather_node = drtl->likelihoodData().getLikelihoodArray(father->getId(), currentNode->getId());
    unique_ptr<TreeLikelihoodInterface::ConstBranchModelIterator> mit(drtl->getNewBranchModelIterator(currentNode->getId()));
    VVVdouble pxy;
    bool first;
    while (mit->hasNext())
    {
      auto bmd = mit->next();
      substitutionCount->setSubstitutionModel(modelSet.getSubstitutionModelForNode(currentNode->getId()));

      // compute all nxy first:
      VVVVdouble nxy(nbClasses);
      for (size_t c = 0; c < nbClasses; ++c)
      {
        VVVdouble& nxy_c = nxy[c];
        double rc = rcRates[c];
        nxy_c.resize(nbTypes);
        for (size_t t = 0; t < nbTypes; ++t)
        {
          VVdouble& nxy_c_t = nxy_c[t];
          auto nijt = substitutionCount->getAllNumbersOfSubstitutions(d * rc, t + 1);
          nxy_c_t.resize(nbStates);
          for (size_t x = 0; x < nbStates; ++x)
          {
            Vdouble& nxy_c_t_x = nxy_c_t[x];
            nxy_c_t_x.resize(nbStates);
            for (size_t y = 0; y < nbStates; ++y)
            {
              nxy_c_t_x[y] = (*nijt)(x, y);
            }
          }
        }
      }

      // Now loop over sites:
      unique_ptr<TreeLikelihoodInterface::SiteIterator> sit(bmd->getNewSiteIterator());
      first = true;
      while (sit->hasNext())
      {
        size_t i = sit->next();
        // We retrieve the transition probabilities and substitution counts for this site partition:
        if (first)
        {
          pxy = drtl->getTransitionProbabilitiesPerRateClass(currentNode->getId(), i);
          first = false;
        }
        const VVdouble& likelihoodsFather_node_i = likelihoodsFather_node[i];
        VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
        for (size_t c = 0; c < nbClasses; ++c)
        {
          const Vdouble& likelihoodsFather_node_i_c = likelihoodsFather_node_i[c];
          Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
          const VVdouble& pxy_c = pxy[c];
          VVVdouble& nxy_c = nxy[c];
          for (size_t x = 0; x < nbStates; ++x)
          {
            double& likelihoodsFatherConstantPart_i_c_x = likelihoodsFatherConstantPart_i_c[x];
            const Vdouble& pxy_c_x = pxy_c[x];
            for (size_t y = 0; y < nbStates; ++y)
            {
              double likelihood_cxy = likelihoodsFatherConstantPart_i_c_x
                                      * pxy_c_x[y]
                                      * likelihoodsFather_node_i_c[y];

              for (size_t t = 0; t < nbTypes; ++t)
              {
                // Now the vector computation:
                substitutionsForCurrentNode[i][t] += likelihood_cxy * nxy_c[t][x][y];
                //                                   <------------>   <------------>
                // Posterior probability                   |                 |
                // for site i and rate class c *           |                 |
                // likelihood for this site----------------+                 |
                //                                                           |
                // Substitution function for site i and rate class c---------+
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
        (*substitutions)(l, i, t) = substitutionsForCurrentNode[rootPatternLinks[i]][t] / Lr[rootPatternLinks[i]];
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

/**************************************************************************************************/

unique_ptr<LegacyProbabilisticSubstitutionMapping> LegacySubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(
  std::shared_ptr<const DRTreeLikelihoodInterface> drtl,
  std::shared_ptr<SubstitutionCountInterface> substitutionCount,
  bool verbose)
{
  // Preamble:
  if (!drtl->isInitialized())
    throw Exception("SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(). Likelihood object is not initialized.");

  // A few variables we'll need:
  const TreeTemplate<Node> tree(drtl->tree());
  auto& sequences = drtl->data();
  auto& rDist = drtl->rateDistribution();

  size_t nbSites         = sequences.getNumberOfSites();
  size_t nbDistinctSites = drtl->likelihoodData().getNumberOfDistinctSites();
  size_t nbStates        = sequences.alphabet().getSize();
  size_t nbClasses       = rDist.getNumberOfCategories();
  size_t nbTypes         = substitutionCount->getNumberOfSubstitutionTypes();
  vector<const Node*> nodes              = tree.getNodes();
  const vector<size_t>& rootPatternLinks = drtl->likelihoodData().getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  size_t nbNodes = nodes.size();

  // We create a new ProbabilisticSubstitutionMapping object:
  auto substitutions = make_unique<LegacyProbabilisticSubstitutionMapping>(tree, substitutionCount, nbSites);

  Vdouble rcRates = rDist.getCategories();

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
      VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
      likelihoodsFatherConstantPart_i.resize(nbClasses);
      for (size_t c = 0; c < nbClasses; ++c)
      {
        Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
        likelihoodsFatherConstantPart_i_c.resize(nbStates);
        double rc = rDist.getProbability(c);
        for (size_t s = 0; s < nbStates; ++s)
        {
          likelihoodsFatherConstantPart_i_c[s] = rc;
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
        const VVVdouble& likelihoodsFather_son = drtl->likelihoodData().getLikelihoodArray(father->getId(), currentSon->getId());

        // Now iterate over all site partitions:
        unique_ptr<TreeLikelihoodInterface::ConstBranchModelIterator> mit(drtl->getNewBranchModelIterator(currentSon->getId()));
        VVVdouble pxy;
        bool first;
        while (mit->hasNext())
        {
          auto bmd = mit->next();
          unique_ptr<TreeLikelihoodInterface::SiteIterator> sit(bmd->getNewSiteIterator());
          first = true;
          while (sit->hasNext())
          {
            size_t i = sit->next();
            // We retrieve the transition probabilities for this site partition:
            if (first)
            {
              pxy = drtl->getTransitionProbabilitiesPerRateClass(currentSon->getId(), i);
              first = false;
            }
            const VVdouble& likelihoodsFather_son_i = likelihoodsFather_son[i];
            VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
            for (size_t c = 0; c < nbClasses; ++c)
            {
              const Vdouble& likelihoodsFather_son_i_c = likelihoodsFather_son_i[c];
              Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
              VVdouble& pxy_c = pxy[c];
              for (size_t x = 0; x < nbStates; ++x)
              {
                Vdouble& pxy_c_x = pxy_c[x];
                double likelihood = 0.;
                for (size_t y = 0; y < nbStates; ++y)
                {
                  likelihood += pxy_c_x[y] * likelihoodsFather_son_i_c[y];
                }
                likelihoodsFatherConstantPart_i_c[x] *= likelihood;
              }
            }
          }
        }
      }
    }
    if (father->hasFather())
    {
      const Node* currentSon = father->getFather();
      const VVVdouble& likelihoodsFather_son = drtl->likelihoodData().getLikelihoodArray(father->getId(), currentSon->getId());
      // Now iterate over all site partitions:
      unique_ptr<TreeLikelihoodInterface::ConstBranchModelIterator> mit(drtl->getNewBranchModelIterator(father->getId()));
      VVVdouble pxy;
      bool first;
      while (mit->hasNext())
      {
        auto bmd = mit->next();
        unique_ptr<TreeLikelihoodInterface::SiteIterator> sit(bmd->getNewSiteIterator());
        first = true;
        while (sit->hasNext())
        {
          size_t i = sit->next();
          // We retrieve the transition probabilities for this site partition:
          if (first)
          {
            pxy = drtl->getTransitionProbabilitiesPerRateClass(father->getId(), i);
            first = false;
          }
          const VVdouble& likelihoodsFather_son_i = likelihoodsFather_son[i];
          VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
          for (size_t c = 0; c < nbClasses; ++c)
          {
            const Vdouble& likelihoodsFather_son_i_c = likelihoodsFather_son_i[c];
            Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
            VVdouble& pxy_c = pxy[c];
            for (size_t x = 0; x < nbStates; ++x)
            {
              double likelihood = 0.;
              for (size_t y = 0; y < nbStates; ++y)
              {
                Vdouble& pxy_c_x = pxy_c[y];
                likelihood += pxy_c_x[x] * likelihoodsFather_son_i_c[y];
              }
              likelihoodsFatherConstantPart_i_c[x] *= likelihood;
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
        vector<double> freqs = drtl->getRootFrequencies(i);
        VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
        for (size_t c = 0; c < nbClasses; ++c)
        {
          Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
          for (size_t x = 0; x < nbStates; ++x)
          {
            likelihoodsFatherConstantPart_i_c[x] *= freqs[x];
          }
        }
      }
    }

    // Then, we deal with the node of interest.
    // We first average uppon 'y' to save computations, and then uppon 'x'.
    // ('y' is the state at 'node' and 'x' the state at 'father'.)

    // Iterate over all site partitions:
    const VVVdouble& likelihoodsFather_node = drtl->likelihoodData().getLikelihoodArray(father->getId(), currentNode->getId());
    unique_ptr<TreeLikelihoodInterface::ConstBranchModelIterator> mit(drtl->getNewBranchModelIterator(currentNode->getId()));
    VVVdouble pxy;
    bool first;
    while (mit->hasNext())
    {
      auto bmd = mit->next();
      substitutionCount->setSubstitutionModel(bmd->getSubstitutionModel());
      // compute all nxy first:
      VVVVdouble nxy(nbClasses);
      for (size_t c = 0; c < nbClasses; ++c)
      {
        double rc = rcRates[c];
        VVVdouble& nxy_c = nxy[c];
        nxy_c.resize(nbTypes);
        for (size_t t = 0; t < nbTypes; ++t)
        {
          VVdouble& nxy_c_t = nxy_c[t];
          nxy_c_t.resize(nbStates);
          auto nijt = substitutionCount->getAllNumbersOfSubstitutions(d * rc, t + 1);
          for (size_t x = 0; x < nbStates; ++x)
          {
            Vdouble& nxy_c_t_x = nxy_c_t[x];
            nxy_c_t_x.resize(nbStates);
            for (size_t y = 0; y < nbStates; ++y)
            {
              nxy_c_t_x[y] = (*nijt)(x, y);
            }
          }
        }
      }

      // Now loop over sites:
      unique_ptr<TreeLikelihoodInterface::SiteIterator> sit(bmd->getNewSiteIterator());
      first = true;
      while (sit->hasNext())
      {
        size_t i = sit->next();
        // We retrieve the transition probabilities and substitution counts for this site partition:
        if (first)
        {
          pxy = drtl->getTransitionProbabilitiesPerRateClass(currentNode->getId(), i);
          first = false;
        }
        const VVdouble& likelihoodsFather_node_i = likelihoodsFather_node[i];
        VVdouble& likelihoodsFatherConstantPart_i = likelihoodsFatherConstantPart[i];
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
          const Vdouble& likelihoodsFather_node_i_c = likelihoodsFather_node_i[c];
          Vdouble& likelihoodsFatherConstantPart_i_c = likelihoodsFatherConstantPart_i[c];
          const VVdouble& pxy_c = pxy[c];
          VVVdouble& nxy_c = nxy[c];
          for (size_t x = 0; x < nbStates; ++x)
          {
            double& likelihoodsFatherConstantPart_i_c_x = likelihoodsFatherConstantPart_i_c[x];
            const Vdouble& pxy_c_x = pxy_c[x];
            for (size_t y = 0; y < nbStates; ++y)
            {
              double likelihood_cxy = likelihoodsFatherConstantPart_i_c_x
                                      * pxy_c_x[y]
                                      * likelihoodsFather_node_i_c[y];
              pairProbabilities(x, y) += likelihood_cxy; // Sum over all rate classes.
              for (size_t t = 0; t < nbTypes; ++t)
              {
                subsCounts[x][y][t] += likelihood_cxy * nxy_c[t][x][y];
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
        (*substitutions)(l, i, t) = substitutionsForCurrentNode[rootPatternLinks[i]][t];
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

/**************************************************************************************************/

unique_ptr<LegacyProbabilisticSubstitutionMapping> LegacySubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(
  shared_ptr<const DRTreeLikelihoodInterface> drtl,
  std::shared_ptr<SubstitutionCountInterface> substitutionCount,
  bool verbose)
{
  // Preamble:
  if (!drtl->isInitialized())
    throw Exception("SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(). Likelihood object is not initialized.");

  // A few variables we'll need:

  const TreeTemplate<Node> tree(drtl->tree());
  auto& sequences = drtl->data();
  auto& rDist = drtl->rateDistribution();
  auto alpha = sequences.getAlphabet();

  size_t nbSites         = sequences.getNumberOfSites();
  size_t nbDistinctSites = drtl->likelihoodData().getNumberOfDistinctSites();
  size_t nbStates        = alpha->getSize();
  size_t nbTypes         = substitutionCount->getNumberOfSubstitutionTypes();
  vector<const Node*> nodes    = tree.getNodes();
  const vector<size_t>& rootPatternLinks = drtl->likelihoodData().getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  size_t nbNodes = nodes.size();

  // We create a new ProbabilisticSubstitutionMapping object:
  auto substitutions = make_unique<LegacyProbabilisticSubstitutionMapping>(tree, substitutionCount, nbSites);

  // Compute the whole likelihood of the tree according to the specified model:

  Vdouble rcRates = rDist.getCategories();

  // Compute the number of substitutions for each class and each branch in the tree:
  if (verbose)
    ApplicationTools::displayTask("Compute marginal ancestral states");
  LegacyMarginalAncestralStateReconstruction masr(drtl);
  map<int, vector<size_t>> ancestors = masr.getAllAncestralStates();

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
    unique_ptr<TreeLikelihoodInterface::ConstBranchModelIterator> mit(drtl->getNewBranchModelIterator(currentNode->getId()));
    while (mit->hasNext())
    {
      auto bmd = mit->next();
      substitutionCount->setSubstitutionModel(bmd->getSubstitutionModel());
      // compute all nxy first:
      VVVdouble nxyt(nbTypes);
      for (size_t t = 0; t < nbTypes; ++t)
      {
        nxyt[t].resize(nbStates);
        auto nxy = substitutionCount->getAllNumbersOfSubstitutions(d, t + 1);
        for (size_t x = 0; x < nbStates; ++x)
        {
          nxyt[t][x].resize(nbStates);
          for (size_t y = 0; y < nbStates; ++y)
          {
            nxyt[t][x][y] = (*nxy)(x, y);
          }
        }
      }
      // Now loop over sites:
      unique_ptr<TreeLikelihoodInterface::SiteIterator> sit(bmd->getNewSiteIterator());
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
        (*substitutions)(l, i, t) = substitutionsForCurrentNode[rootPatternLinks[i]][t];
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

/**************************************************************************************************/

unique_ptr<LegacyProbabilisticSubstitutionMapping> LegacySubstitutionMappingTools::computeSubstitutionVectorsMarginal(
  std::shared_ptr<const DRTreeLikelihoodInterface> drtl,
  std::shared_ptr<SubstitutionCountInterface> substitutionCount,
  bool verbose)
{
  // Preamble:
  if (!drtl->isInitialized())
    throw Exception("SubstitutionMappingTools::computeSubstitutionVectorsMarginal(). Likelihood object is not initialized.");

  // A few variables we'll need:

  const TreeTemplate<Node> tree(drtl->tree());
  auto& sequences = drtl->data();
  auto& rDist = drtl->rateDistribution();

  size_t nbSites         = sequences.getNumberOfSites();
  size_t nbDistinctSites = drtl->likelihoodData().getNumberOfDistinctSites();
  size_t nbStates        = sequences.alphabet().getSize();
  size_t nbClasses       = rDist.getNumberOfCategories();
  size_t nbTypes         = substitutionCount->getNumberOfSubstitutionTypes();
  vector<const Node*> nodes    = tree.getNodes();
  const vector<size_t>& rootPatternLinks = drtl->likelihoodData().getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  size_t nbNodes = nodes.size();

  // We create a new ProbabilisticSubstitutionMapping object:
  auto substitutions = make_unique<LegacyProbabilisticSubstitutionMapping>(tree, substitutionCount, nbSites);

  // Compute the whole likelihood of the tree according to the specified model:

  Vdouble rcProbs = rDist.getProbabilities();
  Vdouble rcRates = rDist.getCategories();

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
    VVVdouble probsNode   = DRTreeLikelihoodTools::getPosteriorProbabilitiesPerStatePerRate(*drtl, currentNode->getId());
    VVVdouble probsFather = DRTreeLikelihoodTools::getPosteriorProbabilitiesPerStatePerRate(*drtl, father->getId());

    // Iterate over all site partitions:
    unique_ptr<TreeLikelihoodInterface::ConstBranchModelIterator> mit(drtl->getNewBranchModelIterator(currentNode->getId()));
    while (mit->hasNext())
    {
      auto bmd = mit->next();
      substitutionCount->setSubstitutionModel(bmd->getSubstitutionModel());
      // compute all nxy first:
      VVVVdouble nxy(nbClasses);
      for (size_t c = 0; c < nbClasses; ++c)
      {
        VVVdouble& nxy_c = nxy[c];
        double rc = rcRates[c];
        nxy_c.resize(nbTypes);
        for (size_t t = 0; t < nbTypes; ++t)
        {
          VVdouble& nxy_c_t = nxy_c[t];
          auto nijt = substitutionCount->getAllNumbersOfSubstitutions(d * rc, t + 1);
          nxy_c_t.resize(nbStates);
          for (size_t x = 0; x < nbStates; ++x)
          {
            Vdouble& nxy_c_t_x = nxy_c_t[x];
            nxy_c_t_x.resize(nbStates);
            for (size_t y = 0; y < nbStates; ++y)
            {
              nxy_c_t_x[y] = (*nijt)(x, y);
            }
          }
        }
      }

      // Now loop over sites:
      unique_ptr<TreeLikelihoodInterface::SiteIterator> sit(bmd->getNewSiteIterator());
      while (sit->hasNext())
      {
        size_t i = sit->next();
        VVdouble& probsNode_i   = probsNode[i];
        VVdouble& probsFather_i = probsFather[i];
        for (size_t c = 0; c < nbClasses; ++c)
        {
          Vdouble& probsNode_i_c   = probsNode_i[c];
          Vdouble& probsFather_i_c = probsFather_i[c];
          VVVdouble& nxy_c = nxy[c];
          for (size_t x = 0; x < nbStates; ++x)
          {
            for (size_t y = 0; y < nbStates; ++y)
            {
              double prob_cxy = probsFather_i_c[x] * probsNode_i_c[y];
              // Now the vector computation:
              for (size_t t = 0; t < nbTypes; ++t)
              {
                substitutionsForCurrentNode[i][t] += prob_cxy * nxy_c[t][x][y];
                //                                   <------>   <------------>
                // Posterior probability                 |             |
                // for site i and rate class c *         |             |
                // likelihood for this site--------------+             |
                //                                                     |
                // Substitution function for site i and rate class c---+
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
        (*substitutions)(l, i, t) = substitutionsForCurrentNode[rootPatternLinks[i]][t];
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

/**************************************************************************************************/

void LegacySubstitutionMappingTools::writeToStream(
  const LegacyProbabilisticSubstitutionMapping& substitutions,
  const SiteContainerInterface& sites,
  size_t type,
  ostream& out)
{
  if (!out)
    throw IOException("LegacySubstitutionMappingTools::writeToFile. Can't write to stream.");
  out << "Branches";
  out << "\tMean";
  for (size_t i = 0; i < substitutions.getNumberOfSites(); i++)
  {
    out << "\tSite" << sites.site(i).getCoordinate();
  }
  out << endl;

  for (size_t j = 0; j < substitutions.getNumberOfBranches(); j++)
  {
    out << substitutions.getNode(j)->getId() << "\t" << substitutions.getNode(j)->getDistanceToFather();
    for (size_t i = 0; i < substitutions.getNumberOfSites(); i++)
    {
      out << "\t" << substitutions(j, i, type);
    }
    out << endl;
  }
}

/**************************************************************************************************/

void LegacySubstitutionMappingTools::readFromStream(
  istream& in, 
  LegacyProbabilisticSubstitutionMapping& substitutions,
  size_t type)
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
    size_t nbBranches = data->getNumberOfRows();
    for (size_t i = 0; i < nbBranches; i++)
    {
      int id = TextTools::toInt(ids[i]);
      size_t br = substitutions.getNodeIndex(id);
      for (size_t j = 0; j < nbSites; j++)
      {
        substitutions(br, j, type) = TextTools::toDouble((*data)(i, j));
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
      substitutions.setSitePosition(i, site);
    }
  }
  catch (Exception& e)
  {
    throw IOException(string("Bad input file. ") + e.what());
  }
}

/**************************************************************************************************/

vector<double> LegacySubstitutionMappingTools::computeTotalSubstitutionVectorForSitePerBranch(
  const LegacySubstitutionMappingInterface& smap,
  size_t siteIndex)
{
  size_t nbBranches = smap.getNumberOfBranches();
  size_t nbTypes    = smap.getNumberOfSubstitutionTypes();
  Vdouble v(nbBranches);
  for (size_t l = 0; l < nbBranches; ++l)
  {
    v[l] = 0;
    for (size_t t = 0; t < nbTypes; ++t)
    {
      v[l] += smap(l, siteIndex, t);
    }
  }
  return v;
}

/**************************************************************************************************/

vector<double> LegacySubstitutionMappingTools::computeTotalSubstitutionVectorForSitePerType(
  const LegacySubstitutionMappingInterface& smap,
  size_t siteIndex)
{
  size_t nbBranches = smap.getNumberOfBranches();
  size_t nbTypes    = smap.getNumberOfSubstitutionTypes();
  Vdouble v(nbTypes);
  for (size_t t = 0; t < nbTypes; ++t)
  {
    v[t] = 0;
    for (size_t l = 0; l < nbBranches; ++l)
      v[t] += smap(l, siteIndex, t);

  }
  return v;
}

/**************************************************************************************************/

double LegacySubstitutionMappingTools::computeNormForSite(
  const LegacySubstitutionMappingInterface& smap,
  size_t siteIndex)
{
  double sumSquare = 0;
  for (size_t l = 0; l < smap.getNumberOfBranches(); ++l)
  {
    double sum = 0;
    for (size_t t = 0; t < smap.getNumberOfSubstitutionTypes(); ++t)
    {
      sum += smap(l, siteIndex, t);
    }
    sumSquare += sum * sum;
  }
  return sqrt(sumSquare);
}

/**************************************************************************************************/

vector<double> LegacySubstitutionMappingTools::computeSumForBranch(
  const LegacySubstitutionMappingInterface& smap,
  size_t branchIndex)
{
  size_t nbSites = smap.getNumberOfSites();
  size_t nbTypes = smap.getNumberOfSubstitutionTypes();
  Vdouble v(nbTypes, 0);
  for (size_t i = 0; i < nbSites; ++i)
  {
    for (size_t t = 0; t < nbTypes; ++t)
    {
      v[t] += smap(branchIndex, i, t);
    }
  }
  return v;
}

/**************************************************************************************************/

vector<double> LegacySubstitutionMappingTools::computeSumForSite(
  const LegacySubstitutionMappingInterface& smap,
  size_t siteIndex)
{
  size_t nbBranches = smap.getNumberOfBranches();
  size_t nbTypes = smap.getNumberOfSubstitutionTypes();
  Vdouble v(nbTypes, 0);
  for (size_t i = 0; i < nbBranches; ++i)
  {
    for (size_t t = 0; t < nbTypes; ++t)
    {
      v[t] += smap(i, siteIndex, t);
    }
  }
  return v;
}

/**************************************************************************************************/

vector<vector<double>> LegacySubstitutionMappingTools::getCountsPerBranch(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  std::shared_ptr<SubstitutionModelInterface> model,
  std::shared_ptr<const SubstitutionRegisterInterface> reg,
  double threshold,
  bool verbose)
{
  auto count = make_shared<UniformizationSubstitutionCount>(model, reg);
  auto mapping = LegacySubstitutionMappingTools::computeSubstitutionVectors(drtl, ids, count, false);

  vector< vector<double> > counts(ids.size());
  size_t nbSites = mapping->getNumberOfSites();
  size_t nbTypes = mapping->getNumberOfSubstitutionTypes();

  for (size_t k = 0; k < ids.size(); ++k)
  {
    vector<double> countsf(nbTypes, 0);
    vector<double> tmp(nbTypes, 0);
    size_t nbIgnored = 0;
    bool error = false;
    for (size_t i = 0; !error && i < nbSites; ++i)
    {
      double s = 0;
      for (size_t t = 0; t < nbTypes; ++t)
      {
        tmp[t] = (*mapping)(mapping->getNodeIndex(ids[k]), i, t);
        error = std::isnan(tmp[t]);
        if (error)
          goto ERROR;
        s += tmp[t];
      }
      if (threshold >= 0)
      {
        if (s <= threshold)
          countsf += tmp;
        else
        {
          nbIgnored++;
        }
      }
      else
      {
        countsf += tmp;
      }
    }

ERROR:
    if (error)
    {
      // We do nothing. This happens for small branches.
      if (verbose)
        ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", counts could not be computed.");
      for (size_t t = 0; t < nbTypes; ++t)
      {
        countsf[t] = 0;
      }
    }
    else
    {
      if (nbIgnored > 0)
      {
        if (verbose)
          ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", " + TextTools::toString(nbIgnored) + " sites (" + TextTools::toString(ceil(static_cast<double>(nbIgnored * 100) / static_cast<double>(nbSites))) + "%) have been ignored because they are presumably saturated.");
      }
    }

    counts[k].resize(countsf.size());
    for (size_t j = 0; j < countsf.size(); ++j)
    {
      counts[k][j] = countsf[j];
    }
  }

  return counts;
}

/**************************************************************************************************/

vector<vector<double>> LegacySubstitutionMappingTools::getCountsPerBranch(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  const SubstitutionModelSet& modelSet,
  std::shared_ptr<const SubstitutionRegisterInterface> reg,
  double threshold,
  bool verbose)
{
  auto count = make_shared<UniformizationSubstitutionCount>(modelSet.getSubstitutionModel(0), reg);

  auto mapping = LegacySubstitutionMappingTools::computeSubstitutionVectors(drtl, modelSet, ids, count, false);

  vector< vector<double> > counts(ids.size());
  size_t nbSites = mapping->getNumberOfSites();
  size_t nbTypes = mapping->getNumberOfSubstitutionTypes();

  for (size_t k = 0; k < ids.size(); ++k)
  {
    vector<double> countsf(nbTypes, 0);
    vector<double> tmp(nbTypes, 0);
    size_t nbIgnored = 0;
    bool error = false;
    for (size_t i = 0; !error && i < nbSites; ++i)
    {
      double s = 0;
      for (size_t t = 0; t < nbTypes; ++t)
      {
        tmp[t] = (*mapping)(mapping->getNodeIndex(ids[k]), i, t);
        error = std::isnan(tmp[t]);
        if (error)
          goto ERROR;
        s += tmp[t];
      }
      if (threshold >= 0)
      {
        if (s <= threshold)
          countsf += tmp;
        else
        {
          nbIgnored++;
        }
      }
      else
      {
        countsf += tmp;
      }
    }

ERROR:
    if (error)
    {
      // We do nothing. This happens for small branches.
      if (verbose)
        ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", counts could not be computed.");
      for (size_t t = 0; t < nbTypes; ++t)
      {
        countsf[t] = 0;
      }
    }
    else
    {
      if (nbIgnored > 0)
      {
        if (verbose)
          ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", " + TextTools::toString(nbIgnored) + " sites (" + TextTools::toString(ceil(static_cast<double>(nbIgnored * 100) / static_cast<double>(nbSites))) + "%) have been ignored because they are presumably saturated.");
      }
    }

    counts[k].resize(countsf.size());
    for (size_t j = 0; j < countsf.size(); ++j)
    {
      counts[k][j] = countsf[j];
    }
  }

  return counts;
}

/**************************************************************************************************/

vector<vector<double>> LegacySubstitutionMappingTools::getNormalizationsPerBranch(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  std::shared_ptr<const SubstitutionModelInterface> nullModel,
  const SubstitutionRegisterInterface& reg,
  bool verbose)
{
  size_t nbTypes = reg.getNumberOfSubstitutionTypes();
  size_t nbStates = nullModel->alphabet().getSize();
  size_t nbSites = drtl->getNumberOfSites();
  vector<int> supportedStates = nullModel->getAlphabetStates();

  // compute the AlphabetIndex for each substitutionType
  vector<shared_ptr<UserAlphabetIndex1>> usai(nbTypes);

  for (size_t nbt = 0; nbt < nbTypes; nbt++) {
    usai[nbt] = make_shared<UserAlphabetIndex1>(nullModel->getAlphabet());
    for (size_t i = 0; i < nbStates; i++)
      usai[nbt]->setIndex(supportedStates[i], 0);
  }

  for (size_t i = 0; i < nbStates; i++)
  {
    for (size_t j = 0; j < nbStates; j++)
    {
      if (i != j)
      {
        size_t nbt = reg.getType(i, j);
        if (nbt != 0)
          usai[nbt - 1]->setIndex(supportedStates[i], usai[nbt - 1]->getIndex(supportedStates[i]) + nullModel->Qij(i, j));
      }
    }
  }

  // compute the normalization for each substitutionType
  vector<vector<double>> rewards(ids.size());

  for (size_t k = 0; k < ids.size(); ++k)
  {
    rewards[k].resize(nbTypes);
  }

  for (size_t nbt = 0; nbt < nbTypes; nbt++)
  {
    auto reward = make_shared<DecompositionReward>(nullModel, usai[nbt]);
    auto mapping = LegacyRewardMappingTools::computeRewardVectors(drtl, ids, reward, false);

    for (size_t k = 0; k < ids.size(); ++k)
    {
      double s = 0;
      for (size_t i = 0; i < nbSites; ++i)
      {
        double tmp = (*mapping)(k, i);
        if (std::isnan(tmp))
        {
          if (verbose)
            ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", reward for type " + reg.getTypeName(nbt + 1) + " could not be computed.");
          s = 0;
          break;
        }
        s += tmp;
      }
      rewards[k][nbt] = s;
    }
  }
  return rewards;
}

/**************************************************************************************************/

vector<vector<double>> LegacySubstitutionMappingTools::getNormalizationsPerBranch(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  std::shared_ptr<const SubstitutionModelSet> nullModelSet,
  const SubstitutionRegisterInterface& reg,
  bool verbose)
{
  size_t nbTypes = reg.getNumberOfSubstitutionTypes();
  size_t nbStates = nullModelSet->alphabet().getSize();
  size_t nbSites = drtl->getNumberOfSites();
  size_t nbModels = nullModelSet->getNumberOfModels();

  // compute the AlphabetIndex for each substitutionType
  // compute the normalization for each substitutionType
  vector<vector<double>> rewards(ids.size());

  for (size_t k = 0; k < ids.size(); ++k)
  {
    rewards[k].resize(nbTypes);
  }

  vector<std::shared_ptr<UserAlphabetIndex1>> usai(nbTypes);
  for (size_t i = 0; i < nbTypes; ++i)
    usai[i] = make_shared<UserAlphabetIndex1>(nullModelSet->getAlphabet());

  for (size_t nbm = 0; nbm < nbModels; nbm++)
  {
    vector<int> mids = VectorTools::vectorIntersection(ids, nullModelSet->getNodesWithModel(nbm));
    
    if (mids.size() > 0)
    {
      auto modn = nullModelSet->getSubstitutionModel(nbm);
      vector<int> supportedStates = modn->getAlphabetStates();

      for (size_t nbt = 0; nbt < nbTypes; nbt++)
        for (size_t i = 0; i < nbStates; i++)
          usai[nbt]->setIndex(supportedStates[i], 0);

      for (size_t i = 0; i < nbStates; i++)
      {
        for (size_t j = 0; j < nbStates; j++)
        {
          if (i != j)
          {
            size_t nbt = reg.getType(i, j);
            if (nbt != 0)
              usai[nbt - 1]->setIndex(supportedStates[i], usai[nbt - 1]->getIndex(supportedStates[i]) + modn->Qij(i, j));
          }
        }
      }

      for (size_t nbt = 0; nbt < nbTypes; nbt++)
      {
        auto reward = make_shared<DecompositionReward>(nullModelSet->getSubstitutionModel(nbm), usai[nbt]);       
        auto mapping = LegacyRewardMappingTools::computeRewardVectors(drtl, mids, reward, false);
        
        for (size_t k = 0; k < mids.size(); k++)
        {
          double s = 0;
          for (size_t i = 0; i < nbSites; ++i)
          {
            double tmp = (*mapping)(mapping->getNodeIndex(mids[k]), i);
            if (std::isnan(tmp))
            {
              if (verbose)
                ApplicationTools::displayWarning("On branch " + TextTools::toString(mids[k]) + ", reward for type " + reg.getTypeName(nbt + 1) + " could not be computed.");
              s = 0;
              break;
            }
            else
              s += tmp;
          }
          
          rewards[VectorTools::which(ids, mids[k])][nbt] = s;
        }
      }
    }
  }

  return rewards;
}

/**************************************************************************************************/

void LegacySubstitutionMappingTools::computeCountsPerTypePerBranch(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  std::shared_ptr<SubstitutionModelInterface> model,
  std::shared_ptr<SubstitutionModelInterface> nullModel,
  std::shared_ptr<const SubstitutionRegisterInterface> reg,
  VVdouble& result,
  bool perTime,
  bool perWord,
  bool verbose)
{
  vector<vector<double>> factors;

  result = getCountsPerBranch(drtl, ids, model, reg, -1, verbose);
  factors = getNormalizationsPerBranch(drtl, ids, nullModel, *reg, verbose);


  // check if perWord

  auto wAlp = dynamic_pointer_cast<const CoreWordAlphabet>(nullModel->getAlphabet());

  float sizeWord = float((wAlp && !perWord) ? wAlp->getLength() : 1);
  

  size_t nbTypes = result[0].size();

  for (size_t k = 0; k < ids.size(); ++k)
  {
    for (size_t t = 0; t < nbTypes; ++t)
    {
      if (factors[k][t] != 0)
        result[k][t] /= (factors[k][t] * sizeWord);
    }
  }

  // Multiply by the lengths of the branches of the input tree

  if (!perTime)
  {
    const TreeTemplate<Node> tree(drtl->tree());

    for (size_t k = 0; k < ids.size(); ++k)
    {
      double l = tree.getNode(ids[k])->getDistanceToFather();
      for (size_t t = 0; t < nbTypes; ++t)
      {
        result[k][t] *= l;
      }
    }
  }
  
}

/**************************************************************************************************/

void LegacySubstitutionMappingTools::computeCountsPerTypePerBranch(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  std::shared_ptr<SubstitutionModelSet> modelSet,
  std::shared_ptr<SubstitutionModelSet> nullModelSet,
  std::shared_ptr<const SubstitutionRegisterInterface> reg,
  VVdouble& result,
  bool perTime,
  bool perWord,
  bool verbose)
{
  vector<vector<double>> factors;

  result = getCountsPerBranch(drtl, ids, *modelSet, reg, -1, verbose);
  factors = getNormalizationsPerBranch(drtl, ids, nullModelSet, *reg, verbose);

  size_t nbTypes = result[0].size();

  // check if perWord

  auto wAlp = dynamic_pointer_cast<const CoreWordAlphabet>(nullModelSet->getAlphabet());

  float sizeWord = float((wAlp && !perWord) ? wAlp->getLength() : 1);
  

  for (size_t k = 0; k < ids.size(); ++k)
  {
    for (size_t t = 0; t < nbTypes; ++t)
    {
      if (factors[k][t] != 0)
        result[k][t] /= (factors[k][t]*sizeWord);
    }
  }

  // Multiply by the lengths of the branches of the input tree

  if (!perTime)
  {
    const TreeTemplate<Node> tree(drtl->tree());

    for (size_t k = 0; k < ids.size(); ++k)
    {
      double l = tree.getNode(ids[k])->getDistanceToFather();
      for (size_t t = 0; t < nbTypes; ++t)
      {
        result[k][t] *= l;
      }
    }
  }
}

/**************************************************************************************************/

void LegacySubstitutionMappingTools::computeCountsPerTypePerBranch(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  std::shared_ptr<SubstitutionModelInterface> model,
  std::shared_ptr<const SubstitutionRegisterInterface> reg,
  VVdouble& result,
  double threshold,
  bool verbose)
{
  result = getCountsPerBranch(drtl, ids, model, reg, threshold, verbose);

  auto creg = dynamic_pointer_cast<const CategorySubstitutionRegister>(reg);

  if (creg && !creg->isStationary())
  {
    size_t nbTypes = result[0].size();

    for (size_t k = 0; k < ids.size(); ++k)
    {
      vector<double> freqs = DRTreeLikelihoodTools::getPosteriorStateFrequencies(*drtl, ids[k]);
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
}

/**************************************************************************************************/

void LegacySubstitutionMappingTools::computeCountsPerSitePerBranch(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  std::shared_ptr<SubstitutionModelInterface> model,
  std::shared_ptr<const SubstitutionRegisterInterface> reg,
  VVdouble& result)
{
  auto count = make_shared<UniformizationSubstitutionCount>(model, reg);
  auto smap = LegacySubstitutionMappingTools::computeSubstitutionVectors(drtl, ids, count, false);

  size_t nbSites = smap->getNumberOfSites();
  size_t nbBr = ids.size();

  VectorTools::resize2(result, nbSites, nbBr);

  vector<size_t> sdi(nbBr);  // reverse of ids
  for (size_t i = 0; i < nbBr; ++i)
  {
    sdi[i] = smap->getNodeIndex(ids[i]);
  }

  for (size_t k = 0; k < nbSites; ++k)
  {
    vector<double> countsf = LegacySubstitutionMappingTools::computeTotalSubstitutionVectorForSitePerBranch(*smap, k);
    Vdouble* resS=&result[k];
    
    for (size_t i = 0; i < nbBr; ++i)
      (*resS)[i]= countsf[sdi[i]];
  }
}


/**************************************************************************************************/

void LegacySubstitutionMappingTools::computeCountsPerSitePerType(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  std::shared_ptr<SubstitutionModelInterface> model,
  std::shared_ptr<const SubstitutionRegisterInterface> reg,
  VVdouble& result)
{
  auto count = make_shared<UniformizationSubstitutionCount>(model, reg);
  auto smap = LegacySubstitutionMappingTools::computeSubstitutionVectors(drtl, ids, count, false);

  size_t nbSites = smap->getNumberOfSites();
  size_t nbTypes = smap->getNumberOfSubstitutionTypes();

  VectorTools::resize2(result, nbSites, nbTypes);

  vector<unique_ptr<LegacyProbabilisticRewardMapping>> rmap;
  
  for (size_t k = 0; k < nbSites; ++k)
  {
    vector<double> countsf = LegacySubstitutionMappingTools::computeTotalSubstitutionVectorForSitePerType(*smap, k);

    Vdouble* resS=&result[k];
    
    for (size_t i = 0; i < nbTypes; ++i)
      (*resS)[i] = countsf[i];
  }
}


/**************************************************************************************************/

void LegacySubstitutionMappingTools::computeCountsPerSitePerType(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  std::shared_ptr<SubstitutionModelInterface> model,
  std::shared_ptr<SubstitutionModelInterface> nullModel,
  std::shared_ptr<const SubstitutionRegisterInterface> reg,
  VVdouble& result,
  bool perTime,
  bool perWord)
{
  
  auto count = make_shared<UniformizationSubstitutionCount>(model, reg);
  auto smap = LegacySubstitutionMappingTools::computeSubstitutionVectors(drtl, ids, count, false);

  size_t nbSites = smap->getNumberOfSites();
  size_t nbTypes = smap->getNumberOfSubstitutionTypes();
  size_t nbStates = nullModel->alphabet().getSize();
  vector<int> supportedStates = nullModel->getAlphabetStates();

  VectorTools::resize2(result, nbSites, nbTypes);

  // compute the AlphabetIndex for each substitutionType
  vector<std::shared_ptr<UserAlphabetIndex1>> usai(nbTypes);

  for (size_t nbt = 0; nbt < nbTypes; nbt++) {
    usai[nbt] = make_shared<UserAlphabetIndex1>(nullModel->getAlphabet());
    for (size_t i = 0; i < nbStates; i++)
      usai[nbt]->setIndex(supportedStates[i], 0);
  }

  for (size_t i = 0; i < nbStates; i++)
  {
    for (size_t j = 0; j < nbStates; j++)
    {
      if (i != j)
      {
        size_t nbt = reg->getType(i, j);
        if (nbt != 0)
          usai[nbt - 1]->setIndex(supportedStates[i], usai[nbt - 1]->getIndex(supportedStates[i]) + nullModel->Qij(i, j));
      }
    }
  }

  // compute the normalization for each site for each substitutionType
  vector<vector<double>> rewards(nbSites);

  for (size_t k = 0; k < nbSites; ++k)
  {
    rewards[k].resize(nbTypes);
  }
  
  for (size_t nbt = 0; nbt < nbTypes; nbt++)
  {
    auto reward = make_shared<DecompositionReward>(nullModel, usai[nbt]);
    auto mapping = LegacyRewardMappingTools::computeRewardVectors(drtl, ids, reward, false);

    for (size_t i = 0; i < nbSites; ++i)
    {
      double s = 0;
      for (size_t k = 0; k < ids.size(); ++k)
      {
        double tmp = (*mapping)(mapping->getNodeIndex(ids[k]), i);
        if (std::isnan(tmp))
        {
          s = 0;
          break;
        }
        s += tmp;
      }
      rewards[i][nbt] = s;
    }
  }

  // Compute the sum of lengths of concerned branchs
  
  double brlen = 0;

  if (!perTime)
  {
    const TreeTemplate<Node> tree(drtl->tree());
    for (size_t k = 0; k < ids.size(); ++k)
      brlen += tree.getNode(ids[k])->getDistanceToFather();
  }
  else
    brlen = 1;

  // check if perWord

  auto wAlp = dynamic_pointer_cast<const CoreWordAlphabet>(nullModel->getAlphabet());

  float sizeWord = float((wAlp && !perWord) ? wAlp->getLength() : 1);
  
  // output
  
  for (size_t k = 0; k < nbSites; ++k)
  {
    vector<double> countsf = LegacySubstitutionMappingTools::computeTotalSubstitutionVectorForSitePerType(*smap, k);

    Vdouble& resS=result[k];
    for (size_t i = 0; i < nbTypes; ++i)
    {
      resS[i] = countsf[i] / rewards[k][i] * brlen / sizeWord;
    }
  }
}


/**************************************************************************************************/

void LegacySubstitutionMappingTools::computeCountsPerSitePerType(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  std::shared_ptr<SubstitutionModelSet> modelSet,
  std::shared_ptr<SubstitutionModelSet> nullModelSet,
  std::shared_ptr<const SubstitutionRegisterInterface> reg,
  VVdouble& result,
  bool perTime,
  bool perWord)
{
  auto count = make_shared<UniformizationSubstitutionCount>(modelSet->getSubstitutionModel(0), reg);
  auto smap = LegacySubstitutionMappingTools::computeSubstitutionVectors(drtl, ids, count, false);

  size_t nbSites = smap->getNumberOfSites();
  size_t nbTypes = smap->getNumberOfSubstitutionTypes();
  size_t nbStates = nullModelSet->alphabet().getSize();
  size_t nbModels = nullModelSet->getNumberOfModels();

  VectorTools::resize2(result, nbSites, nbTypes);
  
  // compute the normalization for each site for each substitutionType
  vector< vector<double> > rewards(nbSites);

  for (size_t k = 0; k < nbSites; ++k)
  {
    rewards[k].resize(nbTypes);
  }
  
  vector<std::shared_ptr<UserAlphabetIndex1>> usai(nbTypes);

  for (size_t nbm = 0; nbm < nbModels; nbm++)
  {
    vector<int> mids = VectorTools::vectorIntersection(ids, nullModelSet->getNodesWithModel(nbm));
    
    if (mids.size()>0)
    {
      const std::shared_ptr<SubstitutionModelInterface> modn = nullModelSet->getSubstitutionModel(nbm);
      vector<int> supportedStates = modn->getAlphabetStates();
      
      for (size_t nbt = 0; nbt < nbTypes; nbt++) {
        usai[nbt] = make_shared<UserAlphabetIndex1>(nullModelSet->getAlphabet());
        for (size_t i = 0; i < nbStates; i++)
          usai[nbt]->setIndex(supportedStates[i], 0);
      }
      
      for (size_t i = 0; i < nbStates; i++)
      {
        for (size_t j = 0; j < nbStates; j++)
        {
          if (i != j)
          {
            size_t nbt = reg->getType(i, j);
            if (nbt != 0)
              usai[nbt - 1]->setIndex(supportedStates[i], usai[nbt - 1]->getIndex(supportedStates[i]) + modn->Qij(i, j));
          }
        }
      }

      for (size_t nbt = 0; nbt < nbTypes; nbt++)
      {
        auto reward = make_shared<DecompositionReward>(nullModelSet->getSubstitutionModel(nbm), usai[nbt]);
        auto mapping = LegacyRewardMappingTools::computeRewardVectors(drtl, mids, reward, false);
        
        for (size_t i = 0; i < nbSites; ++i)
        {
          double s = 0;
          for (size_t k = 0; k < mids.size(); k++)
          {
            double tmp = (*mapping)(mapping->getNodeIndex(mids[k]), i);

            if (std::isnan(tmp))
            {
              s = 0;
              break;
            }
            else
              s += tmp;
          }
          rewards[i][nbt] += s;
        }
      }
    }
  }
  
  // Compute the sum of lengths of concerned branchs
  

  double brlen = 0;

  if (!perTime){
    const TreeTemplate<Node> tree(drtl->tree());
    for (size_t k = 0; k < ids.size(); ++k)
      brlen += tree.getNode(ids[k])->getDistanceToFather();
  }
  else
    brlen = 1.;
  
      
  // check if perWord
  auto wAlp = dynamic_pointer_cast<const CoreWordAlphabet>(nullModelSet->getAlphabet());
  float sizeWord = float((wAlp && !perWord) ? wAlp->getLength() : 1);
 
  // output
  
  for (size_t k = 0; k < nbSites; ++k)
  {
    vector<double> countsf = LegacySubstitutionMappingTools::computeTotalSubstitutionVectorForSitePerType(*smap, k);

    Vdouble& resS = result[k];
    
    for (size_t i = 0; i < nbTypes; ++i)
    {
      resS[i] = countsf[i] / rewards[k][i] * brlen / sizeWord;
    }
  }
}


/**************************************************************************************************/

void LegacySubstitutionMappingTools::computeCountsPerSitePerBranchPerType(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  std::shared_ptr<SubstitutionModelInterface> model,
  std::shared_ptr<const SubstitutionRegisterInterface> reg,
  VVVdouble& result)
{
  auto count = make_shared<UniformizationSubstitutionCount>(model, reg);
  auto smap = LegacySubstitutionMappingTools::computeSubstitutionVectors(drtl, ids, count, false);

  size_t nbSites = smap->getNumberOfSites();
  size_t nbBr = ids.size();
  size_t nbTypes= reg->getNumberOfSubstitutionTypes();

  VectorTools::resize3(result, nbSites, nbBr, nbTypes);
  
  for (size_t i = 0; i < nbTypes; ++i)
  {
    for (size_t j = 0; j < nbSites; ++j)
    {
      VVdouble& resS=result[j];
      
      for (size_t k = 0; k < nbBr; ++k)
      {
        resS[k][i] = (*smap)(smap->getNodeIndex(ids[k]), j, i);
      }
    }
  }
}

/**************************************************************************************************/

void LegacySubstitutionMappingTools::computeCountsPerSitePerBranchPerType(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  std::shared_ptr<SubstitutionModelInterface> model,
  std::shared_ptr<SubstitutionModelInterface> nullModel,
  std::shared_ptr<const SubstitutionRegisterInterface> reg,
  VVVdouble& result,
  bool perTime,
  bool perWord)
{
  auto count = make_shared<UniformizationSubstitutionCount>(model, reg);
  auto smap = LegacySubstitutionMappingTools::computeSubstitutionVectors(drtl, ids, count, false);
  
  size_t nbSites = smap->getNumberOfSites();
  size_t nbTypes = smap->getNumberOfSubstitutionTypes();
  size_t nbStates = nullModel->alphabet().getSize();
  size_t nbBr = ids.size();

  VectorTools::resize3(result, nbSites, nbBr, nbTypes);

  // compute the normalization for each site for each substitutionType
  map<int, vector< vector<double> > > rewards;

  for (auto id : ids)
    VectorTools::resize2(rewards[id], nbSites, nbTypes);

  vector<int> supportedStates = nullModel->getAlphabetStates();

  // compute the AlphabetIndex for each substitutionType
  vector<std::shared_ptr<UserAlphabetIndex1>> usai(nbTypes);

  for (size_t nbt = 0; nbt < nbTypes; nbt++) {
    usai[nbt] = make_shared<UserAlphabetIndex1>(nullModel->getAlphabet());
    for (size_t i = 0; i < nbStates; i++)
      usai[nbt]->setIndex(supportedStates[i], 0);
  }

  for (size_t i = 0; i < nbStates; i++)
  {
    for (size_t j = 0; j < nbStates; j++)
    {
      if (i != j)
      {
        size_t nbt = reg->getType(i, j);
        if (nbt != 0)
          usai[nbt - 1]->setIndex(supportedStates[i], usai[nbt - 1]->getIndex(supportedStates[i]) + nullModel->Qij(i, j));
      }
    }
  }

  for (size_t nbt = 0; nbt < nbTypes; nbt++)
  {
    auto reward = make_shared<DecompositionReward>(nullModel, usai[nbt]);
    auto mapping = LegacyRewardMappingTools::computeRewardVectors(drtl, ids, reward, false);

    for (size_t i = 0; i < nbSites; ++i)
    {
      for (size_t k = 0; k < ids.size(); ++k)
      {
        double tmp = (*mapping)(mapping->getNodeIndex(ids[k]), i);

        if (std::isnan(tmp))
          tmp = 0;
        
        rewards[ids[k]][i][nbt] = tmp;
      }
    }
  }

  // check if perWord
  auto wAlp = dynamic_pointer_cast<const CoreWordAlphabet>(nullModel->getAlphabet());
  float sizeWord = float((wAlp && !perWord) ? wAlp->getLength() : 1);
  

  // output
  const TreeTemplate<Node>& tree = drtl->tree();

  for (size_t i = 0; i < nbTypes; ++i)
  {
    for (size_t j = 0; j < nbSites; ++j)
    {
      VVdouble& resS=result[j];
      
      for (size_t k = 0; k < nbBr; ++k)
      {
        resS[k][i] = (*smap)(smap->getNodeIndex(ids[k]), j, i)/rewards[ids[k]][j][i] * (perTime ? 1 : tree.getNode(ids[k])->getDistanceToFather()) / sizeWord;
      }
    }
  }
}

/**************************************************************************************************/

void LegacySubstitutionMappingTools::computeCountsPerSitePerBranchPerType(
  std::shared_ptr<DRTreeLikelihoodInterface> drtl,
  const vector<int>& ids,
  std::shared_ptr<SubstitutionModelSet> modelSet,
  std::shared_ptr<SubstitutionModelSet> nullModelSet,
  std::shared_ptr<const SubstitutionRegisterInterface> reg,
  VVVdouble& result,
  bool perTime,
  bool perWord)
{
  auto count = make_shared<UniformizationSubstitutionCount>(modelSet->getSubstitutionModel(0), reg);
  auto smap = LegacySubstitutionMappingTools::computeSubstitutionVectors(drtl, ids, count, false);
  
  size_t nbSites = smap->getNumberOfSites();
  size_t nbTypes = smap->getNumberOfSubstitutionTypes();
  size_t nbStates = nullModelSet->alphabet().getSize();
  size_t nbModels = nullModelSet->getNumberOfModels();
  size_t nbBr = ids.size();

  VectorTools::resize3(result, nbSites, nbBr, nbTypes);

  // compute the normalization for each site for each substitutionType
  map<int, vector<vector<double>>> rewards;

  for (auto id : ids)
    VectorTools::resize2(rewards[id], nbSites, nbTypes);
  
  vector<std::shared_ptr<UserAlphabetIndex1>> usai(nbTypes);
  for (size_t i = 0; i < nbTypes; ++i)
    usai[i] = make_shared<UserAlphabetIndex1>(nullModelSet->getAlphabet());

  for (size_t nbm = 0; nbm < nbModels; nbm++)
  {
    vector<int> mids = VectorTools::vectorIntersection(ids, nullModelSet->getNodesWithModel(nbm));
    
    if (mids.size()>0)
    {
      auto modn = nullModelSet->getSubstitutionModel(nbm);
      vector<int> supportedStates = modn->getAlphabetStates();
      
      for (size_t nbt = 0; nbt < nbTypes; nbt++)
        for (size_t i = 0; i < nbStates; i++)
          usai[nbt]->setIndex(supportedStates[i], 0);
      
      for (size_t i = 0; i < nbStates; i++)
      {
        for (size_t j = 0; j < nbStates; j++)
        {
          if (i != j)
          {
            size_t nbt = reg->getType(i, j);
            if (nbt != 0)
              usai[nbt - 1]->setIndex(supportedStates[i], usai[nbt - 1]->getIndex(supportedStates[i]) + modn->Qij(i, j));
          }
        }
      }

      for (size_t nbt = 0; nbt < nbTypes; nbt++)
      {
        auto reward = make_shared<DecompositionReward>(nullModelSet->getSubstitutionModel(nbm), usai[nbt]);
        auto mapping = LegacyRewardMappingTools::computeRewardVectors(drtl, mids, reward, false);
        
        for (size_t i = 0; i < nbSites; ++i)
        {
          for (size_t k = 0; k < mids.size(); k++)
          {
            double tmp =  (*mapping)(mapping->getNodeIndex(mids[k]), i);

            if (std::isnan(tmp))
              tmp = 0;

            rewards[mids[k]][i][nbt] = tmp;
          }
        }
      }
    }
  }

  // check if perWord
  auto wAlp = dynamic_pointer_cast<const CoreWordAlphabet>(nullModelSet->getAlphabet());
  float sizeWord = float((wAlp && !perWord) ? wAlp->getLength() : 1);
  

  // output
  const TreeTemplate<Node>& tree = drtl->tree();

  for (size_t i = 0; i < nbTypes; ++i)
  {
    for (size_t j = 0; j < nbSites; ++j)
    {
      VVdouble& resS=result[j];
      
      for (size_t k = 0; k < nbBr; ++k)
      {
        resS[k][i] = (*smap)(smap->getNodeIndex(ids[k]), j, i)/rewards[ids[k]][j][i]*(perTime?1:tree.getNode(ids[k])->getDistanceToFather())/sizeWord;
      }
    }
  }
}


/**************************************************************************************************/

void LegacySubstitutionMappingTools::outputPerSitePerBranch(
  const string& filename,
  const vector<int>& ids,
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

void LegacySubstitutionMappingTools::outputPerSitePerType(
  const string& filename,
  const SubstitutionRegisterInterface& reg,
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

void LegacySubstitutionMappingTools::outputPerSitePerBranchPerType(
  const string& filenamePrefix,
  const vector<int>& ids,
  const SubstitutionRegisterInterface& reg,
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
    string name = reg.getTypeName(i+1);
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
