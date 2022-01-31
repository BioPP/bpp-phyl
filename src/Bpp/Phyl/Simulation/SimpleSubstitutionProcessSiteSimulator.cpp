//
// File: SimpleSubstitutionProcessSiteSimulator.cpp
// Authors:
//   Julien Dutheil
//   Bastien Boussau
//   Laurent GuÃÂ©guen
// Created: 2004-02-04 16:30:51
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Phyl/Likelihood/ProcessComputationTree.h>
#include <algorithm>

#include "SimpleSubstitutionProcessSiteSimulator.h"

// From SeqLib:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

#include "GivenDataSubstitutionProcessSiteSimulator.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

SimpleSubstitutionProcessSiteSimulator::SimpleSubstitutionProcessSiteSimulator(const SubstitutionProcess& process) :
  process_(&process),
  phyloTree_(process_->getParametrizablePhyloTree()),
  tree_(ProcessComputationTree(*process_)),
  qRates_(),
  qRoots_(),
  seqIndexes_(),
  seqNames_(),
  speciesNodes_(),
  nbNodes_(),
  nbClasses_(process_->getNumberOfClasses()),
  nbStates_(process_->getNumberOfStates()),
  continuousRates_(false),
  outputInternalSites_(false)
{
  init();
}

/******************************************************************************/

void SimpleSubstitutionProcessSiteSimulator::init()
{
  // Initialize sons & fathers of tree_ Nodes
  // set sequence names

  outputInternalSites(outputInternalSites_);

  // Set up cumsum rates

  const auto dRate = process_->getRateDistribution();
  qRates_ = VectorTools::cumSum(dRate->getProbabilities());

  // Initialize root frequencies

  auto cr = VectorTools::cumSum(process_->getRootFrequencies());

  qRoots_.resize(nbClasses_);
  for (auto& roots : qRoots_)
  {
    roots = cr;
  }

  // Initialize cumulative pxy for edges that have models
  auto edges = tree_.getAllEdges();

  for (auto& edge : edges)
  {
    const auto model = edge->getModel();
    if (edge->useProb())
      continue;

    const auto transmodel = dynamic_cast<const TransitionModel*>(model);
    if (!transmodel)
      throw Exception("SubstitutionProcessSiteSimulator::init : model "  + model->getName() + " on branch " + TextTools::toString(tree_.getEdgeIndex(edge)) + " is not a TransitionModel.");

    VVVdouble* cumpxy_node_ = &edge->cumpxy_;
    cumpxy_node_->resize(nbClasses_);

    for (size_t c = 0; c < nbClasses_; c++)
    {
      double brlen = dRate->getCategory(c) * phyloTree_->getEdge(edge->getSpeciesIndex())->getLength();

      VVdouble* cumpxy_node_c_ = &(*cumpxy_node_)[c];

      cumpxy_node_c_->resize(nbStates_);

      // process transition probabilities already consider rates &
      // branch length


      const Matrix<double>* P;

      const auto& vSub(edge->subModelNumbers());

      if (vSub.size() == 0)
        P = &transmodel->getPij_t(brlen);
      else
      {
        if (vSub.size() > 1)
          throw Exception("SubstitutionProcessSiteSimulator::init : only 1 submodel can be used.");

        const auto* mmodel = dynamic_cast<const MixedTransitionModel*>(transmodel);

        const auto* model2 = mmodel->getNModel(vSub[0]);

        P = &model2->getPij_t(brlen);
      }

      for (size_t x = 0; x < nbStates_; x++)
      {
        Vdouble* cumpxy_node_c_x_ = &(*cumpxy_node_c_)[x];
        cumpxy_node_c_x_->resize(nbStates_);
        (*cumpxy_node_c_x_)[0] = (*P)(x, 0);
        for (size_t y = 1; y < nbStates_; y++)
        {
          (*cumpxy_node_c_x_)[y] = (*cumpxy_node_c_x_)[y - 1] + (*P)(x, y);
        }
      }
    }
  }

  // Initialize cumulative prob for mixture nodes
  auto nodes = tree_.getAllNodes();

  for (auto node:nodes)
  {
    if (node->isMixture()) // set probas to chose
    {
      auto outEdges = tree_.getOutgoingEdges(node);
      Vdouble vprob(0);

      for (auto edge : outEdges)
      {
        auto model = dynamic_cast<const MixedTransitionModel*>(edge->getModel());
        if (!model)
          throw Exception("SubstitutionProcessSiteSimulator::init : model in edge " + TextTools::toString(tree_.getEdgeIndex(edge)) + " is not a mixture.");

        const auto& vNb(edge->subModelNumbers());

        double x = 0.;
        for (auto nb:vNb)
        {
          x += model->getNProbability(nb);
        }

        vprob.push_back(x);
        node->sons_.push_back(tree_.getSon(edge));
      }

      vprob /= VectorTools::sum(vprob);

      // Here there is no use to have one cumProb_ per class, but this
      // is used for a posteriori simulations

      node->cumProb_.resize(nbClasses_);
      for (size_t c = 0; c < nbClasses_; c++)
      {
        node->cumProb_[c] = VectorTools::cumSum(vprob);
      }
    }
  }
}

/******************************************************************************/

Site* SimpleSubstitutionProcessSiteSimulator::simulateSite() const
{
  if (continuousRates_ && process_->getRateDistribution())
  {
    double rate = process_->getRateDistribution()->randC();
    return simulateSite(rate);
  }
  else
  {
    size_t rateClass = RandomTools::pickFromCumSum(qRates_);
    return simulateSite(rateClass);
  }
}


Site* SimpleSubstitutionProcessSiteSimulator::simulateSite(double rate) const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  // Use rate class 0

  size_t initialStateIndex = RandomTools::pickFromCumSum(qRoots_[0]);

  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state_ = initialStateIndex;

  evolveInternal(root, rate);

  // Now create a Site object:
  Vint site(seqNames_.size());
  for (size_t i = 0; i < seqNames_.size(); ++i)
  {
    site[i] = process_->getStateMap().getAlphabetStateAsInt(speciesNodes_.at(seqIndexes_[i])->state_);
  }
  return new Site(site, getAlphabet());
}


Site* SimpleSubstitutionProcessSiteSimulator::simulateSite(size_t rateClass) const
{
  // Draw an initial state randomly according to equilibrum frequencies:

  size_t initialStateIndex = RandomTools::pickFromCumSum(qRoots_[rateClass]);

  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state_ = initialStateIndex;

  evolveInternal(root, rateClass);

  // Now create a Site object:
  Vint site(seqNames_.size());
  for (size_t i = 0; i < seqNames_.size(); ++i)
  {
    site[i] = process_->getStateMap().getAlphabetStateAsInt(speciesNodes_.at(seqIndexes_[i])->state_);
  }

  return new Site(site, getAlphabet());
}


Site* SimpleSubstitutionProcessSiteSimulator::simulateSite(size_t ancestralStateIndex, double rate) const
{
  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state_ = ancestralStateIndex;

  evolveInternal(root, rate);

  // Now create a Site object:
  Vint site(seqNames_.size());
  for (size_t i = 0; i < seqNames_.size(); ++i)
  {
    site[i] = process_->getStateMap().getAlphabetStateAsInt(speciesNodes_.at(seqIndexes_[i])->state_);
  }

  return new Site(site, getAlphabet());
}

/******************************************************************************/

SiteSimulationResult* SimpleSubstitutionProcessSiteSimulator::dSimulateSite() const
{
  if (continuousRates_ && process_->getRateDistribution())
  {
    double rate = process_->getRateDistribution()->randC();
    return dSimulateSite(rate);
  }
  else
  {
    size_t rateClass = RandomTools::pickFromCumSum(qRates_);
    return dSimulateSite(rateClass);
  }
}


SiteSimulationResult* SimpleSubstitutionProcessSiteSimulator::dSimulateSite(double rate) const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  // Use rate class 0

  size_t initialStateIndex = RandomTools::pickFromCumSum(qRoots_[0]);

  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state_ = initialStateIndex;

  SiteSimulationResult* ssr = new SiteSimulationResult(phyloTree_, &process_->getStateMap(), initialStateIndex);

  evolveInternal(root, rate, ssr);

  return ssr;
}

SiteSimulationResult* SimpleSubstitutionProcessSiteSimulator::dSimulateSite(size_t rateClass) const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  // Use rate class 0

  size_t initialStateIndex = RandomTools::pickFromCumSum(qRoots_[rateClass]);

  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state_ = initialStateIndex;

  SiteSimulationResult* ssr = new SiteSimulationResult(phyloTree_, &process_->getStateMap(), initialStateIndex);

  evolveInternal(root, rateClass, ssr);

  return ssr;
}

SiteSimulationResult* SimpleSubstitutionProcessSiteSimulator::dSimulateSite(size_t ancestralStateIndex, double rate) const
{
  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state_ = ancestralStateIndex;

  SiteSimulationResult* ssr = new SiteSimulationResult(phyloTree_, &process_->getStateMap(), ancestralStateIndex);

  evolveInternal(root, rate, ssr);

  return ssr;
}


/******************************************************************************/

void SimpleSubstitutionProcessSiteSimulator::evolveInternal(std::shared_ptr<SimProcessNode> node, size_t rateClass, SiteSimulationResult* ssr) const
{
  speciesNodes_[node->getSpeciesIndex()] = node;

  if (node->isSpeciation())
  {
    auto vEdge = tree_.getOutgoingEdges(node);

    for (auto edge : vEdge)
    {
      auto son = tree_.getSon(edge);

      if (edge->getModel())
      {
        if (ssr) // Detailed simulation
        {
          auto tm = dynamic_cast<const SubstitutionModel*>(edge->getModel());

          if (!tm)
            throw Exception("SimpleSubstitutionProcessSiteSimulator::EvolveInternal : detailed simulation not possible for non-markovian model on edge " + TextTools::toString(son->getSpeciesIndex()) + " for model " + edge->getModel()->getName());

          SimpleMutationProcess process(tm);

          double brlen = process_->getRateDistribution()->getCategory(rateClass) * phyloTree_->getEdge(edge->getSpeciesIndex())->getLength();

          MutationPath mp(tm->getAlphabet(), node->state_, brlen);
          if (dynamic_cast<const GivenDataSubstitutionProcessSiteSimulator*>(this) == 0)
          {
            mp = process.detailedEvolve(node->state_, brlen);
            son->state_ = mp.getFinalState();
          }
          else // First get final state
          {
            son->state_ = RandomTools::pickFromCumSum(edge->cumpxy_[rateClass][node->state_]);
            mp = process.detailedEvolve(node->state_, son->state_, brlen);
          }

          // Now append infos in ssr:
          ssr->addNode(edge->getSpeciesIndex(), mp);
        }
        else
          son->state_ = RandomTools::pickFromCumSum(edge->cumpxy_[rateClass][node->state_]);
      }
      else
        son->state_ = node->state_;

      evolveInternal(son, rateClass, ssr);
    }
  }
  else if (node->isMixture())
  {
    const auto& cumProb = node->cumProb_[rateClass];

    size_t y = RandomTools::pickFromCumSum(cumProb);
    auto son = node->sons_[y];
    son->state_ = node->state_;
    evolveInternal(son, rateClass, ssr);
  }
  else
    throw Exception("SimpleSubstitutionProcessSiteSimulator::evolveInternal : unknown property for node " + TextTools::toString(tree_.getNodeIndex(node)));
}

/******************************************************************************/

void SimpleSubstitutionProcessSiteSimulator::evolveInternal(std::shared_ptr<SimProcessNode> node, double rate, SiteSimulationResult* ssr) const
{
  speciesNodes_[node->getSpeciesIndex()] = node;

  if (node->isSpeciation())
  {
    auto vEdge = tree_.getOutgoingEdges(node);

    for (auto edge : vEdge)
    {
      auto son = tree_.getSon(edge);

      if (edge->getModel())
      {
        auto tm = dynamic_cast<const TransitionModel*>(edge->getModel());

        double brlen = rate * phyloTree_->getEdge(edge->getSpeciesIndex())->getLength();


        if (ssr) // Detailed simulation
        {
          auto sm = dynamic_cast<const SubstitutionModel*>(edge->getModel());

          if (!sm)
            throw Exception("SimpleSubstitutionProcessSiteSimulator::EvolveInternal : detailed simulation not possible for non-markovian model on edge " + TextTools::toString(son->getSpeciesIndex()) + " for model " + tm->getName());

          SimpleMutationProcess process(sm);

          MutationPath mp(tm->getAlphabet(), node->state_, brlen);

          if (dynamic_cast<const GivenDataSubstitutionProcessSiteSimulator*>(this) == 0)
          {
            mp = process.detailedEvolve(node->state_, brlen);
            son->state_ = mp.getFinalState();
          }
          else // First get final state
          {
            // Look for the rateClass where rate is (approximation)
            size_t rateClass = process_->getRateDistribution()->getCategoryIndex(rate);

            son->state_ = RandomTools::pickFromCumSum(edge->cumpxy_[rateClass][node->state_]);
            mp = process.detailedEvolve(node->state_, son->state_, brlen);
          }

          // Now append infos in ssr:
          ssr->addNode(edge->getSpeciesIndex(), mp);
        }
        else
        {
          // process transition probabilities already consider rates &
          // branch length

          const Matrix<double>* P;

          const auto& vSub(edge->subModelNumbers());

          if (vSub.size() == 0)
            P = &tm->getPij_t(brlen);
          else
          {
            if (vSub.size() > 1)
              throw Exception("SubstitutionProcessSiteSimulator::init : only 1 submodel can be used.");

            const auto* mmodel = dynamic_cast<const MixedTransitionModel*>(tm);
            const auto* model = mmodel->getNModel(vSub[0]);

            P = &model->getPij_t(brlen);
          }

          double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
          for (size_t y = 0; y < nbStates_; y++)
          {
            rand -= (*P)(node->state_, y);
            if (rand <= 0)
            {
              son->state_ = y;
              break;
            }
          }
        }
      }
      else
      {
        son->state_ = node->state_;
      }

      evolveInternal(son, rate, ssr);
    }
  }
  else if (node->isMixture())
  {
    const auto& cumProb = node->cumProb_[0]; // index 0 because it is only possible in a priori simulations, ie all class mixture probabilities are the same

    size_t y = RandomTools::pickFromCumSum(cumProb);
    auto son = node->sons_[y];
    son->state_ = node->state_;
    evolveInternal(son, rate, ssr);
  }
  else
    throw Exception("SimpleSubstitutionProcessSiteSimulator::evolveInternal : unknown property for node " + TextTools::toString(tree_.getNodeIndex(node)));
}


void SimpleSubstitutionProcessSiteSimulator::outputInternalSites(bool yn)
{
  outputInternalSites_ = yn;

  if (outputInternalSites_)
  {
    auto vCN = phyloTree_->getAllNodes();
    seqNames_.resize(vCN.size());
    seqIndexes_.resize(vCN.size());
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      auto index = phyloTree_->getNodeIndex(vCN[i]);
      seqNames_[i] = (phyloTree_->isLeaf(vCN[i])) ? vCN[i]->getName() : TextTools::toString(index);
      seqIndexes_[i] = index;
    }
  }
  else
  {
    auto vCN = phyloTree_->getAllLeaves();
    seqNames_.resize(vCN.size());
    seqIndexes_.resize(vCN.size());
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      seqNames_[i] = vCN[i]->getName();
      seqIndexes_[i] = phyloTree_->getNodeIndex(vCN[i]);
    }
  }
}
