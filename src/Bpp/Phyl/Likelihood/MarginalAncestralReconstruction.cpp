// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/VectorTools.h>

#include "MarginalAncestralReconstruction.h"

using namespace bpp;
using namespace std;

vector<size_t> MarginalAncestralReconstruction::getAncestralStatesForNode(uint nodeId, VVdouble& probs, bool sample) const
{
  vector<size_t> ancestors(nbSites_);

  auto vv = likelihood_->getLikelihoodsAtNode(nodeId)->targetValue();

  probs.resize(nbSites_);

  for (uint i = 0; i < static_cast<uint>(nbSites_); ++i)
  {
    copyEigenToBpp(vv.col(i) / vv.col(i).sum(), probs[static_cast<size_t>(i)]);
  }

  if (sample)
  {
    for (size_t i = 0; i < nbSites_; ++i)
    {
      const auto& coli = probs[i];
      double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
      for (size_t j = 0; j < nbStates_; ++j)
      {
        r -= coli[j];
        if (r < 0)
        {
          ancestors[i] = j;
          break;
        }
      }
    }
  }
  else
  {
    size_t pos;
    for (size_t i = 0; i < nbSites_; ++i)
    {
      vv.col(Eigen::Index(i)).maxCoeff(&pos);
      ancestors[i] = pos;
    }
  }
  return ancestors;
}

map<uint, vector<size_t>> MarginalAncestralReconstruction::getAllAncestralStates() const
{
  map<uint, vector<size_t>> ancestors;
  // Clone the data into a AlignedSequenceContainer for more efficiency:
  shared_ptr<AlignmentDataInterface> data = make_shared<AlignedSequenceContainer>(dynamic_cast<const SiteContainerInterface&>(likelihood_->shrunkData()));
  recursiveMarginalAncestralStates(tree_->getRoot(), ancestors, *data);
  return ancestors;
}

unique_ptr<Sequence> MarginalAncestralReconstruction::getAncestralSequenceForNode(uint nodeId, VVdouble* probs, bool sample) const
{
  string name = tree_->getNode(nodeId)->hasName() ? tree_->getNode(nodeId)->getName() : ("" + TextTools::toString(nodeId));
  vector<int> allStates(nbSites_);

  const auto& stateMap = likelihood_->stateMap();

  VVdouble patternedProbs;

  if (probs)
  {
    auto states = getAncestralStatesForNode(nodeId, *probs, sample);
    for (size_t i = 0; i < nbSites_; ++i)
    {
      allStates[i] = stateMap.getAlphabetStateAsInt(states[i]);
    }
  }
  else
  {
    auto states = getAncestralStatesForNode(nodeId, patternedProbs, sample);
    for (size_t i = 0; i < nbSites_; ++i)
    {
      allStates[i] = stateMap.getAlphabetStateAsInt(states[i]);
    }
  }

  return make_unique<Sequence>(name, allStates, alphabet_);
}

void MarginalAncestralReconstruction::recursiveMarginalAncestralStates(
    const std::shared_ptr<PhyloNode> node,
    map<uint, vector<size_t>>& ancestors,
    AlignmentDataInterface& data) const
{
  if (tree_->isLeaf(node))
  {
    try
    {
      auto& sc = dynamic_cast<const SiteContainerInterface&>(data);
      const Sequence& seq = sc.sequence(node->getName());
      vector<size_t>* v = &ancestors[tree_->getNodeIndex(node)];
      v->resize(seq.size());
      // This is a tricky way to store the real sequence as an ancestral one...
      // In case of Markov Modulated models, we consider that the real sequences
      // are all in the first category.
      const auto& stateMap = likelihood_->stateMap();
      for (size_t i = 0; i < seq.size(); ++i)
      {
        (*v)[i] = stateMap.getModelStates(seq[i])[0];
      }
    }
    catch (bad_cast&)
    {
      ancestors[tree_->getNodeIndex(node)] = getAncestralStatesForNode(tree_->getNodeIndex(node));
    }
  }
  else
  {
    ancestors[tree_->getNodeIndex(node)] = getAncestralStatesForNode(tree_->getNodeIndex(node));
    vector<shared_ptr<PhyloNode>> vsons = tree_->getSons(node);

    for (size_t i = 0; i < vsons.size(); i++)
    {
      recursiveMarginalAncestralStates(vsons[i], ancestors, data);
    }
  }
}

unique_ptr<AlignedSequenceContainer> MarginalAncestralReconstruction::getAncestralSequences(bool sample) const
{
  auto asc = make_unique<AlignedSequenceContainer>(alphabet_);
  vector<shared_ptr<PhyloNode>> inNodes = tree_->getAllInnerNodes();
  for (size_t i = 0; i < inNodes.size(); ++i)
  {
    auto seq = getAncestralSequenceForNode(tree_->getNodeIndex(inNodes[i]), nullptr, sample);
    asc->addSequence(seq->getName(), seq);
  }
  return asc;
}
