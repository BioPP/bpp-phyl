//
// File: HomogeneousSequenceSimulator.h
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

#ifndef _HOMOGENEOUSSEQUENCESIMULATOR_H_
#define _HOMOGENEOUSSEQUENCESIMULATOR_H_

#include "PreciseSequenceSimulator.h"
#include "SiteSimulator.h"
#include "TreeTemplate.h"
#include "SubstitutionModel.h"

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

class HomogeneousSequenceSimulationResult: public SequenceSimulationResult
{
	protected:
		double _rate;
		
	public:
		HomogeneousSequenceSimulationResult(const TreeTemplate<Node> * tree, int ancestralState, double rate):
			SequenceSimulationResult(tree, ancestralState),
			_rate(rate) {}

		virtual ~HomogeneousSequenceSimulationResult() {}
	
	public:
		virtual double getRate() const { return _rate; }
};

//---------------------------------------------------------------------------

/**
 * @brief Site and sequences simulation undes homogeneous models.
 *
 * Rate across sites variation is supported, using a DiscreteDistribution object or by specifying explicitely the rate of the sites to simulate.
 */
class HomogeneousSequenceSimulator: public virtual PreciseSequenceSimulator, public virtual SiteSimulator
{
	protected:
		const MutationProcess * _process;
		const Alphabet * _alphabet;
		const SubstitutionModel * _model;
		const DiscreteDistribution * _rate;
		const TreeTemplate<Node> * _tree;
	
		/**
		 * @brief This stores once for all all leaves in a given order.
		 * This order will be used during site creation.
		 */
		vector<const Node *> _leaves;
	
		vector<string> _seqNames;

		unsigned int _nbNodes;
		unsigned int _nbClasses;
		unsigned int _nbStates;
	
		/**
		 * @brief Stores intermediate results.
		 */
		mutable map<const Node *, int> _states;
		mutable map<const Node *, const Vint *> _multipleStates;
		mutable HomogeneousSequenceSimulationResult * _hssr;
		mutable map<const Node *, VVVdouble> _cumpxy;
	
	public:
		
		HomogeneousSequenceSimulator(
			const MutationProcess * process,
			const DiscreteDistribution * rate,
			const TreeTemplate<Node> * tree,
			bool verbose = true
		);
			
		virtual ~HomogeneousSequenceSimulator() {}

	public:
	
		/**
		 * @name The SequenceSimulator interface
		 *
		 * @{
		 */
		Site * simulate() const;
		Site * simulate(int ancestralState) const;
		Site * simulate(int ancestralState, double rate) const;
		Site * simulate(double rate) const;
		SiteContainer * simulate(unsigned int numberOfSites) const;
		/** @} */
		virtual Site * simulate(int ancestralState, unsigned int rateClass) const;

		/**
		 * @name the PreciseSequenceSimulator interface.
		 *
		 * @{
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
		const TreeTemplate<Node> * getTree() const;
	
	protected:
		
		/**
		 * @brief Evolve from an initial state along a branch, knowing the evolutionary rate class.
		 *
		 * This method is fast since all pijt have been computed in the constructor of the class.
		 */
		int evolve(int initialState, const Node * node, unsigned int rateClass) const;
		
		/**
		 * @brief Evolve from an initial state along a branch, knowing the evolutionary rate.
		 *
		 * This method is slower than the privious one since exponential terms must be computed.
		 */
		int evolve(int initialState, const Node * node, double rate) const;
		
		void evolveInternal(const Node * node, unsigned int rateClass) const;
		void evolveInternal(const Node * node, double rate) const;
		
		void multipleEvolve(const Vint & initialState, const Node * node, const vector<unsigned int> & rateClasses, Vint & finalStates) const;
		void multipleEvolveInternal(const Node * node, const vector<unsigned int> & rateClasses) const;
		SiteContainer * multipleEvolve(const Vint & initialStates, const vector<unsigned int> & rateClasses) const;
		
		void preciseEvolve(int initialState, unsigned int rateClass) const;
		void preciseEvolveInternal(const Node * node, unsigned int rateClass) const;

};

#endif //_HOMOGENEOUSSEQUENCESIMULATOR_H_

