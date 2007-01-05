//
// File: MarginalAncestralStateReconstruction.cpp
// Created by: Julien Dutheil
// Created on: Fri Jul 08 13:32 2005
//

/*
Copyright or � or Copr. CNRS, (November 16, 2004)

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
#include <NumCalc/VectorTools.h>
using namespace VectorFunctions;

vector<int> MarginalAncestralStateReconstruction::getAncestralStatesForNode(const Node * node) const
{
	vector<int> ancestors(_nDistinctSites);
	if(node -> isLeaf()) {
		VVdouble larray = _likelihood -> getLikelihoodData() -> getLeafLikelihoods(node);
		for(unsigned int i = 0; i < _nDistinctSites; i++) {
			ancestors[i] = (int)whichmax(larray[i]);
		}
	} else {
		VVVdouble larray = _likelihood -> computeLikelihoodAtNode(node);
		Vdouble freqs    = _likelihood -> getSubstitutionModel() -> getFrequencies();
		Vdouble rcProbs  = _likelihood -> getRateDistribution() -> getProbabilities(); 
		for(unsigned int i = 0; i < _nDistinctSites; i++) {
			Vdouble likelihoods(_nStates, 0);
			VVdouble * larray_i = & larray[i];
			for(unsigned int c = 0; c < _nClasses; c++) {
				Vdouble * larray_i_c = & (* larray_i)[c];
				double rcp = rcProbs[c];
				for(unsigned int x = 0; x < _nStates; x++) {
					likelihoods[x] += (* larray_i_c)[x] * freqs[x] * rcp;
				}
			}
			ancestors[i] = (int)whichmax(likelihoods);
		}
	}
	return ancestors;
}

map<const Node *, vector<int> > MarginalAncestralStateReconstruction::getAllAncestralStates() const
{
	map<const Node *, vector<int> > ancestors;
	// Clone the data into a AlignedSequenceContainer for more efficiency:
	AlignedSequenceContainer * data = new AlignedSequenceContainer(* _likelihood->getLikelihoodData()->getShrunkData());
	recursiveMarginalAncestralStates(dynamic_cast<const TreeTemplate<Node> *>(_likelihood->getTree())->getRootNode(), ancestors, *data);
	delete data;
	return ancestors;
}

Sequence * MarginalAncestralStateReconstruction::getAncestralSequenceForNode(const Node * node) const
{
	string name = node->hasName() ? node->getName() : "" + node->getId();
	vector<int> states = getAncestralStatesForNode(node);
	vector<int> allStates(_nSites);
	const vector<unsigned int> * rootPatternLinks = &_likelihood->getLikelihoodData()->getRootArrayPositions();
	for(unsigned int i = 0; i < _nSites; i++)
  {
		allStates[i] = states[(* rootPatternLinks)[i]];
	}
	return new Sequence(name, allStates, _alphabet);
}

void MarginalAncestralStateReconstruction::recursiveMarginalAncestralStates(
			const Node * node,
			map<const Node *, vector<int> > & ancestors,
			AlignedSequenceContainer & data) const
{
	if(node->isLeaf())
  {
		ancestors[node] = data.getSequence(node->getName())->getContent();
	}
  else
  {
		ancestors[node] = getAncestralStatesForNode(node);
		for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
    {
			recursiveMarginalAncestralStates(node->getSon(i), ancestors, data);
		}
	}
}

