//
// File: MarginalAncestralStateReconstruction.cpp
// Created by: Julien Dutheil
// Created on: Fri Jul 08 13:32 2005
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "MarginalAncestralStateReconstruction.h"
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

using namespace bpp;
using namespace std;

vector<size_t> MarginalAncestralStateReconstruction::getAncestralStatesForNode(int nodeId, VVdouble& probs, bool sample) const
{
  vector<size_t> ancestors(nbDistinctSites_);
  probs.resize(nbDistinctSites_);
  double cumProb = 0;
  double r;
  if (likelihood_->getTree().isLeaf(nodeId))
  {
    VVdouble larray = likelihood_->getLikelihoodData()->getLeafLikelihoods(nodeId);
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

map<int, vector<size_t> > MarginalAncestralStateReconstruction::getAllAncestralStates() const
{
  map<int, vector<size_t> > ancestors;
  // Clone the data into a AlignedValuesContainer for more efficiency:
  shared_ptr<AlignedValuesContainer> data(dynamic_cast<AlignedValuesContainer*>(likelihood_->getLikelihoodData()->getShrunkData()->clone()));
  recursiveMarginalAncestralStates(tree_.getRootNode(), ancestors, *data);
  return ancestors;
}

Sequence* MarginalAncestralStateReconstruction::getAncestralSequenceForNode(int nodeId, VVdouble* probs, bool sample) const
{
  string name = tree_.hasNodeName(nodeId) ? tree_.getNodeName(nodeId) : ("" + TextTools::toString(nodeId));
  const vector<size_t>* rootPatternLinks = &likelihood_->getLikelihoodData()->getRootArrayPositions();
  const TransitionModel* model = likelihood_->getModelForSite(tree_.getNodesId()[0], 0); // We assume all nodes have a model with the same number of states.
  vector<size_t> states;
  vector<int> allStates(nbSites_);
  VVdouble patternedProbs;
  if (probs)
  {
    states = getAncestralStatesForNode(nodeId, patternedProbs, sample);
    probs->resize(nbSites_);
    for (size_t i = 0; i < nbSites_; i++)
    {
      allStates[i] = model->getAlphabetStateAsInt(states[(*rootPatternLinks)[i]]);
      (*probs)[i] = patternedProbs[(*rootPatternLinks)[i]];
    }
  }
  else
  {
    states = getAncestralStatesForNode(nodeId, patternedProbs, sample);
    for (size_t i = 0; i < nbSites_; i++)
    {
      allStates[i] = model->getAlphabetStateAsInt(states[(*rootPatternLinks)[i]]);
    }
  }
  return new BasicSequence(name, allStates, alphabet_);
}

void MarginalAncestralStateReconstruction::recursiveMarginalAncestralStates(
  const Node* node,
  map<int, vector<size_t> >& ancestors,
  AlignedValuesContainer& data) const
{  
  if (node->isLeaf())
  {
    const SiteContainer* sc=dynamic_cast<const SiteContainer*>(&data);
    if (sc)
    {
      const Sequence& seq = sc->getSequence(node->getName());
      vector<size_t>* v = &ancestors[node->getId()];
      v->resize(seq.size());
    // This is a tricky way to store the real sequence as an ancestral one...
    // In case of Markov Modulated models, we consider that the real sequences
    // Are all in the first category.
      const TransitionModel* model = likelihood_->getModelForSite(tree_.getNodesId()[0], 0); // We assume all nodes have a model with the same number of states.
      for (size_t i = 0; i < seq.size(); i++)
      {
        (*v)[i] = model->getModelStates(seq[i])[0];
      }
    }
    else
    {
      ancestors[node->getId()] = getAncestralStatesForNode(node->getId());
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

AlignedSequenceContainer* MarginalAncestralStateReconstruction::getAncestralSequences(bool sample) const
{
  AlignedSequenceContainer* asc = new AlignedSequenceContainer(alphabet_);
  vector<int> ids = tree_.getInnerNodesId();
  for (size_t i = 0; i < ids.size(); i++)
  {
    Sequence* seq = getAncestralSequenceForNode(ids[i], NULL, sample);
    asc->addSequence(*seq);
    delete seq;
  }
  return asc;
}

