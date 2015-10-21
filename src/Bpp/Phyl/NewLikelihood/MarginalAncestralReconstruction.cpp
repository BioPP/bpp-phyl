//
// File: MarginalAncestralReconstruction.cpp
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

#include "MarginalAncestralReconstruction.h"
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

using namespace bpp;
using namespace std;

vector<size_t> MarginalAncestralReconstruction::getAncestralStatesForNode(int nodeId, VVdouble& probs, bool sample) const
{
  vector<size_t> ancestors(nbDistinctSites_);
  probs.resize(nbDistinctSites_);
  double cumProb = 0;
  double r;

  likelihood_->computeLikelihoodsAtNode(nodeId);

  if (!likelihood_->getLikelihoodData().getNodeData(nodeId, 0).usesLog())
    probs=likelihood_->getLikelihoodData().getNodeData(nodeId, 0).getLikelihoodArray(ComputingNode::D0)*likelihood_->getSubstitutionProcess()->getProbabilityForModel(0);
  else
    probs=VectorTools::exp(likelihood_->getLikelihoodData().getNodeData(nodeId, 0).getLikelihoodArray(ComputingNode::D0) + log(likelihood_->getSubstitutionProcess()->getProbabilityForModel(0)));
  
  for (size_t c = 1; c < nbClasses_; c++)
  {
    if (!likelihood_->getLikelihoodData().getNodeData(nodeId, c).usesLog())
      probs+=likelihood_->getLikelihoodData().getNodeData(nodeId, c).getLikelihoodArray(ComputingNode::D0)*likelihood_->getSubstitutionProcess()->getProbabilityForModel(c);
    else
      probs+=VectorTools::exp(likelihood_->getLikelihoodData().getNodeData(nodeId, 0).getLikelihoodArray(ComputingNode::D0) + log(likelihood_->getSubstitutionProcess()->getProbabilityForModel(0)));
  } 

  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    Vdouble* probs_i = &probs[i];
    double s=VectorTools::sum(*probs_i);
    (*probs_i)/=s;
  }
    
  if (sample)
  { 
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      Vdouble* probs_i = &probs[i];
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
  }
  else
    for (size_t i = 0; i < nbDistinctSites_; i++)
      ancestors[i] = VectorTools::whichMax(probs[i]);

  return ancestors;
}

map<int, vector<size_t> > MarginalAncestralReconstruction::getAllAncestralStates() const
{
  map<int, vector<size_t> > ancestors;
  // Clone the data into a AlignedSequenceContainer for more efficiency:
  AlignedSequenceContainer* data = new AlignedSequenceContainer(*likelihood_->getLikelihoodData().getShrunkData());
  recursiveMarginalAncestralStates(tree_.getRootNode(), ancestors, *data);
  delete data;
  return ancestors;
}

Sequence* MarginalAncestralReconstruction::getAncestralSequenceForNode(int nodeId, VVdouble* probs, bool sample) const
{
  string name = tree_.hasNodeName(nodeId) ? tree_.getNodeName(nodeId) : ("" + TextTools::toString(nodeId));
  const vector<size_t>* rootPatternLinks = &likelihood_->getLikelihoodData().getRootArrayPositions();
  const SubstitutionModel& model = likelihood_->getSubstitutionProcess()->getSubstitutionModel(tree_.getNodesId()[0], 0); // We assume all nodes have a model with the same set of states.
  vector<size_t> states;
  vector<int> allStates(nbSites_);
  VVdouble patternedProbs;
  if (probs)
  {
    states = getAncestralStatesForNode(nodeId, patternedProbs, sample);
    probs->resize(nbSites_);
    for (size_t i = 0; i < nbSites_; i++)
    {
      allStates[i] = model.getAlphabetStateAsInt(states[(*rootPatternLinks)[i]]);
      (*probs)[i] = patternedProbs[(*rootPatternLinks)[i]];
    }
  }
  else
  {
    states = getAncestralStatesForNode(nodeId, patternedProbs, sample);
    for (size_t i = 0; i < nbSites_; i++)
    {
      allStates[i] = model.getAlphabetStateAsInt(states[(*rootPatternLinks)[i]]);
    }
  }
  return new BasicSequence(name, allStates, alphabet_);
}

void MarginalAncestralReconstruction::recursiveMarginalAncestralStates(
  const Node* node,
  map<int, vector<size_t> >& ancestors,
  AlignedSequenceContainer& data) const
{
  if (node->isLeaf())
  {
    const Sequence& seq = data.getSequence(node->getName());
    vector<size_t>* v = &ancestors[node->getId()];
    v->resize(seq.size());
    // This is a tricky way to store the real sequence as an ancestral one...
    // In case of Markov Modulated models, we consider that the real sequences
    // Are all in the first category.
    const SubstitutionModel& model = likelihood_->getSubstitutionProcess()->getSubstitutionModel(tree_.getNodesId()[0], 0); // We assume all nodes have a model with the same number of states.
    for (size_t i = 0; i < seq.size(); i++)
    {
      (*v)[i] = model.getModelStates(seq[i])[0];
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

AlignedSequenceContainer* MarginalAncestralReconstruction::getAncestralSequences(bool sample) const
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

