//
// File: DetailedSequenceSimulator.h
// Created by: Julien Dutheil
// Created on: Wed Aug  24 15:20 2005
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _DETAILEDSEQUENCESIMULATOR_H_
#define _DETAILEDSEQUENCESIMULATOR_H_

#include "SequenceSimulator.h"
#include "TreeTemplate.h"
#include "MutationProcess.h"

// From the STL:
#include <map>
#include <vector>
using namespace std;

/**
 * @brief Data structure to store the result of a DetailedSequenceSimulator.
 *
 * This data structure stores each transitional state, and the time when it occured.
 */
class SequenceSimulationResult
{
	protected:
		mutable map<const Node *, unsigned int> _indexes;
		unsigned int _currentIndex;
		vector<MutationPath> _paths;
		vector<int> _ancestralStates;
		const TreeTemplate<Node> * _tree;
		vector<const Node *> _leaves;
		
	public:
		SequenceSimulationResult(const TreeTemplate<Node> * tree, int ancestralState):
			_currentIndex(0) {
			_tree = tree;
			_indexes[tree -> getRootNode()] = 0;
			_ancestralStates.push_back(ancestralState);
			_leaves = tree -> getLeaves();
		}

		virtual ~SequenceSimulationResult() {}
	
	public:
		virtual void addNode(const Node * node, MutationPath path) {
			_currentIndex++;
			_indexes[node] = _currentIndex;
			_paths.push_back(path);
			_ancestralStates.push_back(path.getFinalState());
		}

		virtual int getAncestralState(unsigned int i)    const { return _ancestralStates[i]; }

		virtual int getAncestralState(const Node * node) const { return _ancestralStates[_indexes[node]]; }

		virtual unsigned int getSubstitutionCount(unsigned int i)    const { return _paths[i].getNumberOfEvents(); }
		
		virtual unsigned int getSubstitutionCount(const Node * node) const { return _paths[_indexes[node]].getNumberOfEvents(); }
		
		virtual vector<double> getSubstitutionVector() const
		{
			unsigned int n = _paths.size();
			vector<double> counts(n);
			for(unsigned int i = 0; i < n; i++) counts[i] = (double)_paths[i].getNumberOfEvents();
			return counts;
		}

		virtual vector<int> getFinalStates() const
		{
			unsigned int n = _leaves.size(); 
			vector<int> states(n);
			for(unsigned int i = 0; i < n; i++) {
				states[i] = _ancestralStates[_indexes[_leaves[i]]];
			}
			return states;
		}

};

//---------------------------------------------------------------------------

/**
 * @brief This interface adds the dSimulate method to the SequenceSimulator interface.
 *
 * Instances of this class should be used when a detailed output of the simulation is needed.
 */
class DetailedSequenceSimulator: public virtual SequenceSimulator
{
	public:
		DetailedSequenceSimulator() {}
		virtual ~DetailedSequenceSimulator() {}
	
	public:
		/**
		 * @brief Get a detailed simulation result for one site.
		 *
		 * @return A SequenceSimulationResult object with all ancestral
		 * states for all nodes and branches.
		 */
		virtual SequenceSimulationResult * dSimulate() const = 0;
		
};

#endif // _DETAILEDSEQUENCESIMULATOR_H_

