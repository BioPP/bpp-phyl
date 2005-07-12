//
// File: SequenceSimulator.h
// Created by: Julien Dutheil
// Created on: Wed Feb  4 16:30:51 2004
//

/*
Copyright ou © ou Copr. CNRS, (16 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant à fournir des classes
pour l'analyse de données phylogénétiques.

Ce logiciel est régi par la licence CeCILL soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accepté les
termes.
*/

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

Julien.Dutheil@univ-montp2.fr

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

#ifndef _SEQUENCESIMULATOR_H_
#define _SEQUENCESIMULATOR_H_

#include "TreeTemplate.h"
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
class PreciseSequenceSimulator: public SequenceSimulator
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
		HomogeneousSequenceSimulationResult(const TreeTemplate<Node> * tree, int ancestralState, double rate):
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
		mutable map<const Node *, Vint *> _multipleStates;
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
		SiteContainer * simulate(unsigned int numberOfSites) const;
		/** @} */

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
		int evolve(int initialState, const Node * node, int rateClass) const;
		void multipleEvolve(Vint & initialState, const Node * node, Vint & rateClasses, Vint & finalStates) const;
		Site * evolve(int initialState, int rateClass) const;
		SiteContainer * multipleEvolve(Vint & initialStates, Vint & rateClasses) const;
		void preciseEvolve(int initialState, int rateClass) const;
		void evolveInternal(const Node * node, int rateClass) const;
		void multipleEvolveInternal(const Node * node, Vint & rateClasses) const;
		void preciseEvolveInternal(const Node * node, int rateClass) const;

};

//---------------------------------------------------------------------------

#endif	//_SEQUENCESIMULATOR_H_

