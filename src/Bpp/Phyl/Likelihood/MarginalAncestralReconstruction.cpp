//
// File: MarginalAncestralReconstruction.cpp
// Authors:
//   Julien Dutheil
// Created: 2005-07-08 13:32:00
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

#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/VectorTools.h>

#include "MarginalAncestralReconstruction.h"

using namespace bpp;
using namespace std;

vector<size_t> MarginalAncestralReconstruction::getAncestralStatesForNode(uint nodeId, VVdouble& probs, bool sample) const
{
  vector<size_t> ancestors(nbSites_);

  auto vv = likelihood_->getLikelihoodsAtNode(nodeId)->getTargetValue();

  probs.resize(nbSites_);

  for (auto i = 0; i < nbSites_; i++)
  {
    copyEigenToBpp(vv.col(i) / vv.col(i).sum(), probs[size_t(i)]);
  }

  if (sample)
  {
    for (size_t i = 0; i < nbSites_; i++)
    {
      const auto& coli = probs[i];
      double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
      for (size_t j = 0; j < nbStates_; j++)
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
    for (size_t i = 0; i < nbSites_; i++)
    {
      vv.col(Eigen::Index(i)).maxCoeff(&pos);
      ancestors[i] = pos;
    }
  }
  return ancestors;
}

map<uint, vector<size_t> > MarginalAncestralReconstruction::getAllAncestralStates() const
{
  map<uint, vector<size_t> > ancestors;
  // Clone the data into a AlignedSequenceContainer for more efficiency:
  shared_ptr<AlignedValuesContainer> data(dynamic_cast<AlignedValuesContainer*>(likelihood_->getShrunkData()->clone()));
  recursiveMarginalAncestralStates(tree_->getRoot(), ancestors, *data);
  return ancestors;
}

Sequence* MarginalAncestralReconstruction::getAncestralSequenceForNode(uint nodeId, VVdouble* probs, bool sample) const
{
  string name = tree_->getNode(nodeId)->hasName() ? tree_->getNode(nodeId)->getName() : ("" + TextTools::toString(nodeId));
  vector<int> allStates(nbSites_);

  const auto& rootPatternLinks = likelihood_->getRootArrayPositions();

  const auto& statemap = likelihood_->getStateMap();

  VVdouble patternedProbs;

  if (probs)
  {
    auto states = getAncestralStatesForNode(nodeId, patternedProbs, sample);
    probs->resize(nbSites_);
    for (size_t i = 0; i < nbSites_; i++)
    {
      allStates[i] = statemap.getAlphabetStateAsInt(states[rootPatternLinks(Eigen::Index(i))]);
      (*probs)[i] = patternedProbs[rootPatternLinks(Eigen::Index(i))];
    }
  }
  else
  {
    auto states = getAncestralStatesForNode(nodeId, patternedProbs, sample);
    for (size_t i = 0; i < nbSites_; i++)
    {
      allStates[i] = statemap.getAlphabetStateAsInt(states[rootPatternLinks(Eigen::Index(i))]);
    }
  }

  return new BasicSequence(name, allStates, alphabet_);
}

void MarginalAncestralReconstruction::recursiveMarginalAncestralStates(
  const std::shared_ptr<PhyloNode> node,
  map<uint, vector<size_t> >& ancestors,
  AlignedValuesContainer& data) const
{
  if (tree_->isLeaf(node))
  {
    const SiteContainer* sc = dynamic_cast<const SiteContainer*>(&data);
    if (sc)
    {
      const Sequence& seq = sc->getSequence(node->getName());
      vector<size_t>* v = &ancestors[tree_->getNodeIndex(node)];
      v->resize(seq.size());
      // This is a tricky way to store the real sequence as an ancestral one...
      // In case of Markov Modulated models, we consider that the real sequences
      // are all in the first category.
      const auto& statemap = likelihood_->getStateMap();
      for (size_t i = 0; i < seq.size(); i++)
      {
        (*v)[i] = statemap.getModelStates(seq[i])[0];
      }
    }
    else
    {
      ancestors[tree_->getNodeIndex(node)] = getAncestralStatesForNode(tree_->getNodeIndex(node));
    }
  }
  else
  {
    ancestors[tree_->getNodeIndex(node)] = getAncestralStatesForNode(tree_->getNodeIndex(node));
    vector<shared_ptr<PhyloNode> > vsons = tree_->getSons(node);

    for (size_t i = 0; i < vsons.size(); i++)
    {
      recursiveMarginalAncestralStates(vsons[i], ancestors, data);
    }
  }
}

AlignedSequenceContainer* MarginalAncestralReconstruction::getAncestralSequences(bool sample) const
{
  AlignedSequenceContainer* asc = new AlignedSequenceContainer(alphabet_);
  vector<shared_ptr<PhyloNode> > inNodes = tree_->getAllInnerNodes();
  for (size_t i = 0; i < inNodes.size(); i++)
  {
    Sequence* seq = getAncestralSequenceForNode(tree_->getNodeIndex(inNodes[i]), NULL, sample);
    asc->addSequence(*seq);
    delete seq;
  }
  return asc;
}
