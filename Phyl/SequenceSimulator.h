//
// File: SequenceSimulator.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Feb  4 16:30:51 2004
//

#ifndef _SEQUENCESIMULATOR_H_
#define _SEQUENCESIMULATOR_H_

#include "Tree.h"
#include "SubstitutionModel.h"
#include "MutationProcess.h"

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/Site.h>
#include <Seq/SiteContainer.h>

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/RandomTools.h>

// From the STL:
#include <map>
#include <vector>
using namespace std;


//---------------------------------------------------------------------------

class SequenceSimulationResult
{
	protected:
		mutable map<const Node *, unsigned int> _indexes;
		unsigned int _currentIndex;
		vector<MutationPath> _paths;
		vector<int> _ancestralStates;
		const Tree * _tree;
		vector<const Node *> _leaves;
		
	public:
		SequenceSimulationResult(const Tree * tree, int ancestralState):
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

class SequenceSimulator
{
	public:
		SequenceSimulator() {}
		virtual ~SequenceSimulator() {}
	
	public:
		virtual Site * simulate() const = 0;
		virtual SiteContainer * simulate(unsigned int numberOfSites) const = 0;
		virtual SequenceSimulationResult * preciseSimulate() const = 0;
	
};

//---------------------------------------------------------------------------

/**
 * @brief This interface adds the preciseSimulate method to the SequenceSimulator interface.
 *
 * Instances of this class should be used when a detailed output of the simulation is needed.
 */
class PreciseSequenceSimulator
{
	public:
		PreciseSequenceSimulator() {}
		virtual ~PreciseSequenceSimulator() {}
	
	public:
		/**
		 * @brief Get a detailed simulation result for one site.
		 *
		 * @return A SequenceSimulationResult object with all ancestral
		 * states for all nodes and branches.
		 */
		virtual SequenceSimulationResult * preciseSimulate() const = 0;
	
};

//---------------------------------------------------------------------------

class HomogeneousSequenceSimulationResult: public SequenceSimulationResult
{
	protected:
		double _rate;
		
	public:
		HomogeneousSequenceSimulationResult(const Tree * tree, int ancestralState, double rate):
			SequenceSimulationResult(tree, ancestralState),
			_rate(rate) {}

		virtual ~HomogeneousSequenceSimulationResult() {}
	
	public:
		virtual double getRate() const { return _rate; }
};

//---------------------------------------------------------------------------

class HomogeneousSequenceSimulator: public PreciseSequenceSimulator
{
	protected:
		const MutationProcess * _process;
		const Alphabet * _alphabet;
		const SubstitutionModel * _model;
		const DiscreteDistribution * _rate;
		const Tree * _tree;
	
		/**
		 * @brief This stores once for all all leaves in a given order.
		 * This order will be used during site creation.
		 */
		vector<const Node *> _leaves;
	
		vector<string> _seqNames;
	
		/**
		 * @brief Stores intermediate results.
		 */
		mutable map<const Node *, int> _states;
		mutable HomogeneousSequenceSimulationResult * _hssr;
	
	public:
		
		HomogeneousSequenceSimulator(
			const MutationProcess * process,
			const DiscreteDistribution * rate,
			const Tree * tree
		);
			
		virtual ~HomogeneousSequenceSimulator() {}

	public:
	
		/**
		 * @name The SequenceSimulator interface
		 *
		 * @{
		 */
		Site * simulate() const;
		SiteContainer * simulate(unsigned int numberOfSites) const;
		/** @} */

		/**
		 * @names the PreciseSequenceSimulator interface.
		 */
		SequenceSimulationResult * preciseSimulate() const;
		/** @} */
	
		/**
		 * @brief Get the mutation process associated to this instance.
		 *
		 * @return The MutationProcess object associated to this instance.
		 */
		const MutationProcess * getMutationProcess() const;
		
		/**
		 * @brief Get the substitution model associated to this instance.
		 *
		 * @return The SubstitutionModel object associated to this instance.
		 */
		const SubstitutionModel * getSubstitutionModel() const;
		
		/**
		 * @brief Get the rate distribution associated to this instance.
		 *
		 * @return The DiscreteDistribution object associated to this instance.
		 */
		const DiscreteDistribution * getRateDistribution() const;

		/**
		 * @brief Get the tree associated to this instance.
		 *
		 * @return The Tree object associated to this instance.
		 */
		const Tree * getTree() const;
	
	protected:
		Site * evolve(int initialState, double rate = 1.) const;
		void preciseEvolve(int initialState, double rate = 1.) const;
		void evolveInternal(const Node * node, double rate) const;
		void preciseEvolveInternal(const Node * node, double rate) const;

};

//---------------------------------------------------------------------------

#endif	//_SEQUENCESIMULATOR_H_
