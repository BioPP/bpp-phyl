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

class SequenceSimulator
{
	public:
		SequenceSimulator() {}
		virtual ~SequenceSimulator() {}
	
	public:
		virtual Site * simulate() const = 0;
		virtual SiteContainer * simulate(unsigned int numberOfSites) const = 0;
	
};

class HomogeneousSequenceSimulator: public SequenceSimulator
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
	
	public:
		
		HomogeneousSequenceSimulator(
			const MutationProcess * process,
			const DiscreteDistribution * rate,
			const Tree * tree);
			
		virtual ~HomogeneousSequenceSimulator() {}

	public:
		Site * evolve(int initialState, double rate = 1.) const;
		
		/**
		 * @name The SequenceSimulator interface
		 *
		 * @{
		 */
		Site * simulate() const;
		SiteContainer * simulate(unsigned int numberOfSites) const;
		/** @} */
	
		const MutationProcess * getMutationProcess() const;
		const SubstitutionModel * getSubstitutionModel() const;
		const DiscreteDistribution * getRateDistribution() const;
		const Tree * getTree() const;
	
	protected:
		void evolveInternal(const Node * node, double rate) const;

};

#endif	//_SEQUENCESIMULATOR_H_
