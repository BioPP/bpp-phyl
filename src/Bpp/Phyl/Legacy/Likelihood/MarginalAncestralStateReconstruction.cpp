// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "MarginalAncestralStateReconstruction.h"
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

using namespace bpp;
using namespace std;

vector<size_t> LegacyMarginalAncestralStateReconstruction::getAncestralStatesForNode(int nodeId, VVdouble& probs, bool sample) const
{
  vector<size_t> ancestors(nbDistinctSites_);
  probs.resize(nbDistinctSites_);
  double cumProb = 0;
  double r;
  if (likelihood_->tree().isLeaf(nodeId))
  {
    VVdouble larray = likelihood_->likelihoodData().getLeafLikelihoods(nodeId);
    for (size_t i = 0; i < nbDistinctSites_; ++i)
    {
      Vdouble* probs_i = &probs[i];
      probs_i->resize(nbStates_);
      size_t j = VectorTools::whichMax(larray[i]);
      ancestors[i] = j;
      (*probs_i)[j] = 1.;
    }
  }
  else
  {
    VVVdouble larray;

    likelihood_->computeLikelihoodAtNode(nodeId, larray);
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble* larray_i = &larray[i];
      Vdouble* probs_i = &probs[i];
      probs_i->resize(nbStates_);
      for (size_t c = 0; c < nbClasses_; c++)
      {
        Vdouble* larray_i_c = &(*larray_i)[c];
        for (size_t x = 0; x < nbStates_; x++)
        {
          (*probs_i)[x] += (*larray_i_c)[x] * r_[c] / l_[i];
        }
      }
      if (sample)
      {
        cumProb = 0;
        r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
        for (size_t j = 0; j < nbStates_; j++)
        {
          cumProb += (*probs_i)[j];
          if (r <= cumProb)
          {
            ancestors[i] = j;
            break;
          }
        }
      }
      else
        ancestors[i] = VectorTools::whichMax(*probs_i);
    }
  }
  return ancestors;
}

map<int, vector<size_t>> LegacyMarginalAncestralStateReconstruction::getAllAncestralStates() const
{
  map<int, vector<size_t>> ancestors;
  // Clone the data into a AlignedSequenceContainer for more efficiency:
  auto data = make_unique<AlignedSequenceContainer>(dynamic_cast<const SiteContainerInterface&>(likelihood_->likelihoodData().shrunkData()));
  recursiveMarginalAncestralStates(tree_.getRootNode(), ancestors, *data);
  return ancestors;
}

std::unique_ptr<Sequence> LegacyMarginalAncestralStateReconstruction::getAncestralSequenceForNode(int nodeId, VVdouble* probs, bool sample) const
{
  string name = tree_.hasNodeName(nodeId) ? tree_.getNodeName(nodeId) : ("" + TextTools::toString(nodeId));
  const vector<size_t>& rootPatternLinks = likelihood_->likelihoodData().getRootArrayPositions();
  auto model = likelihood_->getModelForSite(tree_.getNodesId()[0], 0); // We assume all nodes have a model with the same number of states.
  vector<size_t> states;
  vector<int> allStates(nbSites_);
  VVdouble patternedProbs;
  if (probs)
  {
    states = getAncestralStatesForNode(nodeId, patternedProbs, sample);
    probs->resize(nbSites_);
    for (size_t i = 0; i < nbSites_; i++)
    {
      allStates[i] = model->getAlphabetStateAsInt(states[rootPatternLinks[i]]);
      (*probs)[i] = patternedProbs[rootPatternLinks[i]];
    }
  }
  else
  {
    states = getAncestralStatesForNode(nodeId, patternedProbs, sample);
    for (size_t i = 0; i < nbSites_; i++)
    {
      allStates[i] = model->getAlphabetStateAsInt(states[rootPatternLinks[i]]);
    }
  }
  std::shared_ptr<const Alphabet> alpha = alphabet_; //Copy needed because of const
  return make_unique<Sequence>(name, allStates, alpha);
}

void LegacyMarginalAncestralStateReconstruction::recursiveMarginalAncestralStates(
  const Node* node,
  map<int, vector<size_t> >& ancestors,
  AlignedSequenceContainer& data) const
{
  if (node->isLeaf())
  {
    const Sequence& seq = data.sequence(node->getName());
    vector<size_t>* v = &ancestors[node->getId()];
    v->resize(seq.size());
    // This is a tricky way to store the real sequence as an ancestral one...
    // In case of Markov Modulated models, we consider that the real sequences
    // Are all in the first category.
    auto model = likelihood_->getModelForSite(tree_.getNodesId()[0], 0); // We assume all nodes have a model with the same number of states.
    for (size_t i = 0; i < seq.size(); i++)
    {
      (*v)[i] = model->getModelStates(seq[i])[0];
    }
  }
  else
  {
    ancestors[node->getId()] = getAncestralStatesForNode(node->getId());
    for (size_t i = 0; i < node->getNumberOfSons(); i++)
    {
      recursiveMarginalAncestralStates(node->getSon(i), ancestors, data);
    }
  }
}

unique_ptr<SiteContainerInterface> LegacyMarginalAncestralStateReconstruction::getAncestralSequences(bool sample) const
{
  unique_ptr<SiteContainerInterface> asc = make_unique<AlignedSequenceContainer>(alphabet_);
  vector<int> ids = tree_.getInnerNodesId();
  for (size_t i = 0; i < ids.size(); i++)
  {
    auto seq = getAncestralSequenceForNode(ids[i], nullptr, sample);
    asc->addSequence(seq->getName(), seq);
  }
  return asc;
}

