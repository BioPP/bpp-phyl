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
#include "../Likelihood/MarginalAncestralStateReconstruction.h"

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
  const DRTreeLikelihood& drtl,
  const vector<int>& nodeIds,
  Reward& reward,
  bool verbose) throw (Exception)
{
  // Preamble:
  if (!drtl.isInitialized())
    throw Exception("RewardMappingTools::computeRewardVectors(). Likelihood object is not initialized.");

  // A few variables we'll need:

  const TreeTemplate<Node> tree(drtl.getTree());
  const SiteContainer*    sequences = drtl.getData();
  const DiscreteDistribution* rDist = drtl.getRateDistribution();

  size_t nbSites         = sequences->getNumberOfSites();
  size_t nbDistinctSites = drtl.getLikelihoodData()->getNumberOfDistinctSites();
  size_t nbStates        = sequences->getAlphabet()->getSize();
  size_t nbClasses       = rDist->getNumberOfCategories();
  vector<const Node*> nodes    = tree.getNodes();
  const vector<size_t>* rootPatternLinks
    = &drtl.getLikelihoodData()->getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  size_t nbNodes         = nodes.size();

  // We create a new ProbabilisticRewardMapping object:
  ProbabilisticRewardMapping* rewards = new ProbabilisticRewardMapping(tree, &reward, nbSites);

  // Store likelihood for each rate for each site:
  VVVdouble lik;
  drtl.computeLikelihoodAtNode(tree.getRootId(), lik);
  Vdouble Lr(nbDistinctSites, 0);
  Vdouble rcProbs = rDist->getProbabilities();
  Vdouble rcRates = rDist->getCategories();
  for (size_t i = 0; i < nbDistinctSites; i++)
  {
    VVdouble* lik_i = &lik[i];
    for (size_t c = 0; c < nbClasses; c++)
    {
      Vdouble* lik_i_c = &(*lik_i)[c];
      double rc = rDist->getProbability(c);
      for (size_t s = 0; s < nbStates; s++)
      {
        Lr[i] += (*lik_i_c)[s] * rc;
      }
    }
  }

  // Compute the reward for each class and each branch in the tree:
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
    Vdouble rewardsForCurrentNode(nbDistinctSites);

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
        double rc = rDist->getProbability(c);
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
        const VVVdouble* likelihoodsFather_son = &drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentSon->getId());

        // Now iterate over all site partitions:
        auto_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(currentSon->getId()));
        VVVdouble pxy;
        bool first;
        while (mit->hasNext())
        {
          TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
          auto_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
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
            for (size_t c = 0; c < nbClasses; c++)
            {
              const Vdouble* likelihoodsFather_son_i_c = &(*likelihoodsFather_son_i)[c];
              Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
              VVdouble* pxy_c = &pxy[c];
              for (size_t x = 0; x < nbStates; x++)
              {
                Vdouble* pxy_c_x = &(*pxy_c)[x];
                double likelihood = 0.;
                for (size_t y = 0; y < nbStates; y++)
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
      auto_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(father->getId()));
      VVVdouble pxy;
      bool first;
      while (mit->hasNext())
      {
        TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
        auto_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
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
          for (size_t c = 0; c < nbClasses; c++)
          {
            const Vdouble* likelihoodsFather_son_i_c = &(*likelihoodsFather_son_i)[c];
            Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
            VVdouble* pxy_c = &pxy[c];
            for (size_t x = 0; x < nbStates; x++)
            {
              double likelihood = 0.;
              for (size_t y = 0; y < nbStates; y++)
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
      for (size_t i = 0; i < nbDistinctSites; i++)
      {
        vector<double> freqs = drtl.getRootFrequencies(i);
        VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
        for (size_t c = 0; c < nbClasses; c++)
        {
          Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
          for (size_t x = 0; x < nbStates; x++)
          {
            (*likelihoodsFatherConstantPart_i_c)[x] *= freqs[x];
          }
        }
      }
    }

    // Then, we deal with the node of interest.
    // We first average upon 'y' to save computations, and then upon 'x'.
    // ('y' is the state at 'node' and 'x' the state at 'father'.)

    // Iterate over all site partitions:
    const VVVdouble* likelihoodsFather_node = &(drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentNode->getId()));
    auto_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(currentNode->getId()));
    VVVdouble pxy;
    bool first;
    while (mit->hasNext())
    {
      TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
      reward.setSubstitutionModel(bmd->getModel());
      // compute all nxy first:
      VVVdouble nxy(nbClasses);
      for (size_t c = 0; c < nbClasses; ++c)
      {
        VVdouble* nxy_c = &nxy[c];
        double rc = rcRates[c];
        Matrix<double>* nij = reward.getAllRewards(d * rc);
        nxy_c->resize(nbStates);
        for (size_t x = 0; x < nbStates; ++x)
        {
          Vdouble* nxy_c_x = &(*nxy_c)[x];
          nxy_c_x->resize(nbStates);
          for (size_t y = 0; y < nbStates; ++y)
          {
            (*nxy_c_x)[y] = (*nij)(x, y);
          }
        }
        delete nij;
      }

      // Now loop over sites:
      auto_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
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
        for (size_t c = 0; c < nbClasses; ++c)
        {
          const Vdouble* likelihoodsFather_node_i_c = &(*likelihoodsFather_node_i)[c];
          Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
          const VVdouble* pxy_c = &pxy[c];
          VVdouble* nxy_c = &nxy[c];
          for (size_t x = 0; x < nbStates; ++x)
          {
            double* likelihoodsFatherConstantPart_i_c_x = &(*likelihoodsFatherConstantPart_i_c)[x];
            const Vdouble* pxy_c_x = &(*pxy_c)[x];
            for (size_t y = 0; y < nbStates; ++y)
            {
              double likelihood_cxy = (*likelihoodsFatherConstantPart_i_c_x)
                                      * (*pxy_c_x)[y]
                                      * (*likelihoodsFather_node_i_c)[y];

              // Now the vector computation:
              rewardsForCurrentNode[i] += likelihood_cxy * (*nxy_c)[x][y];
              //                       <------------>   <--------------->
              // Posterior probability         |                 |
              // for site i and rate class c * |                 |
              // likelihood for this site------+                 |
              //                                                 |
              // Reward function for site i and rate class c------+
            }
          }
        }
      }
    }

    // Now we just have to copy the substitutions into the result vector:
    for (size_t i = 0; i < nbSites; ++i)
    {
      (*rewards)(l, i) = rewardsForCurrentNode[(*rootPatternLinks)[i]] / Lr[(*rootPatternLinks)[i]];
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
  const SiteContainer& sites,
  ostream& out)
throw (IOException)
{
  if (!out)
    throw IOException("RewardMappingTools::writeToFile. Can't write to stream.");
  out << "Branches";
  out << "\tMean";
  for (size_t i = 0; i < rewards.getNumberOfSites(); i++)
  {
    out << "\tSite" << sites.getSite(i).getPosition();
  }
  out << endl;

  for (size_t j = 0; j < rewards.getNumberOfBranches(); j++)
  {
    out << rewards.getNode(j)->getId() << "\t" << rewards.getNode(j)->getDistanceToFather();
    for (size_t i = 0; i < rewards.getNumberOfSites(); i++)
    {
      out << "\t" << rewards(j, i);
    }
    out << endl;
  }
}

/**************************************************************************************************/

void RewardMappingTools::readFromStream(istream& in, ProbabilisticRewardMapping& rewards)
throw (IOException)
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
    size_t nbBranches = data->getNumberOfRows();
    for (size_t i = 0; i < nbBranches; i++)
    {
      int id = TextTools::toInt(ids[i]);
      size_t br = rewards.getNodeIndex(id);
      for (size_t j = 0; j < nbSites; j++)
      {
        rewards(br, j) = TextTools::toDouble((*data)(i, j));
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

    delete data;
  }
  catch (Exception& e)
  {
    throw IOException(string("Bad input file. ") + e.what());
  }
}

/**************************************************************************************************/

double RewardMappingTools::computeSumForBranch(const RewardMapping& smap, size_t branchIndex)
{
  size_t nbSites = smap.getNumberOfSites();
  double v = 0;
  for (size_t i = 0; i < nbSites; ++i)
  {
    v += smap(branchIndex, i);
  }
  return v;
}

/**************************************************************************************************/

double RewardMappingTools::computeSumForSite(const RewardMapping& smap, size_t siteIndex)
{
  size_t nbBranches = smap.getNumberOfBranches();
  double v = 0;
  for (size_t i = 0; i < nbBranches; ++i)
  {
    v += smap(i, siteIndex);
  }
  return v;
}

/**************************************************************************************************/
