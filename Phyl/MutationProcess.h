//
// File: MutationProcess.h
// Created by: jdutheil <julien.dutheil@ens-lyon.fr>
// Created on: Wed Mar 12 16:11:44 2003
//

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
#ifndef _MUTATIONPROCESS_H_
#define _MUTATIONPROCESS_H_

#include "SubstitutionModel.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>

/**
 * @brief This class is used by MutationProcess to store detailed results of simulations.
 */
class MutationPath {
	
	protected:

		/**
		 * @brief The states taken, without intiial state.
		 */
		vector<int>    _states;

		/**
		 * @brief Times between states.
		 * The first element in array is the time between the initial state and the first state in _states.
		 */
		vector<double> _times;
		
		/**
		 * @brief The initial state.
		 */
		int _initialState;

		/**
		 * @brief Total time of evolution.
		 * Typically, this is a branch length.
		 */
		double _totalTime;

	public:

		/**
		 * @brief Builds a new MutationPath object with initial state 'initialState' and total time 'time'.
		 *
		 * @param initialState The initial state.
		 * @param time         The total time of evolution.
		 */
		MutationPath(int initialState, double time):
			_initialState(initialState), _totalTime(time) {};

		~MutationPath() {};

	public:
		
		/**
		 * @brief Add a new mutation event.
		 *
		 * @param state The new state after mutation event.
		 * @param time  The time between this mutation and previous mutation (or initial state).
		 */
		void addEvent(int state, double time) {
			_states.push_back(state);
			_times.push_back(time);
		}

		/**
		 * @brief Retrieve the initial state.
		 *
		 * @return The initial state of this path.
		 */
		int getInitialState() const { return _initialState; }

		/**
		 * @brief Retrieve the total time of evolution.
		 *
		 * @return The total time of evolution.
		 */
		double getTotalTime() const { return _totalTime; }
		
		/**
		 * @brief Retrieve the number of mutation events.
		 *
		 * @return The numbe rof mutation events, i.e. the numer of states (without initial state).
		 */
		unsigned int getNumberOfEvents() const { return _states.size(); }

		/**
		 * @brief Retrieve the final state of this path.
		 *
		 * @return The initial state if no mutation occured, otherwise sends the state after last mutation event.
		 */
		int getFinalState() const {
			if(_states.size() == 0) return _initialState;
			else return _states[_states.size() - 1];
		}
};

/**
 * @brief Base class for simulations.
 *
 * A mutation process defines the rules for mutations to occure.
 * The MutationProcess interface provides two methods, one for mutating a character in
 * state i in another character, another for achieving this task n times.
 */
class MutationProcess {

	public:
		MutationProcess() {};
		virtual ~MutationProcess() {};
	
	public:
		
    /**
     * @brief Mutate a character in state i.
		 *
		 * @param state The current state of the character.
     */
    virtual int mutate(int state) const = 0;

    /**
     * @brief Mutate a character in state i n times.
     * 
		 * @param state The current state of the character.
		 * @param n The number of mutations to perform.
     */
    virtual int mutate(int state, int n) const = 0;
	
		/**
		 * @brief Get the time before next mutation event.
		 *
		 * @param state The actual state of the chain;
		 * @return A random time before next mutation event.
		 */
		virtual double getTimeBeforeNextMutationEvent(int state) const = 0;
		
		/**
		 * @brief Simulation a character evolution during a specified time
		 * according to the given substitution model and send the final state.
		 *
		 * @param initialState The state before beginning evolution.
		 * @param time         The time during which evolution must occure.
		 * @return The resulting state after evolution is completed.
		 */
		virtual int evolve(int initialState, double time) const = 0;
	
		/**
		 * @brief Simulation a character evolution during a specified time
		 * according to the given substitution model and send the total path
		 * with all intermediate states and times between mutation events.
		 *
		 * @param initialState The state before beginning evolution.
		 * @param time         The time during which evolution must occure.
		 * @return The resulting mutation path.
		 */
		virtual MutationPath detailedEvolve(int initialState, double time) const = 0;

		/**
		 * @brief Get the substitution model associated to the mutation process.
		 *
		 * @return The SubstitutionModel associated to this instance.
		 */
		virtual const SubstitutionModel * getSubstitutionModel() const = 0;
};

/**
 * @brief Low implementation of the MutationProcess class.
 *
 * This class provides an implementation of the MutationProcess interface.
 * It assumes that there are _size states allowed for the character of interest,
 * and that the distribution of probabilities are in _repartition.
 * As a matter of facts, probabilities must be cumulative ones, so that _repartition
 * contains values of the repartition function.
 * The mutate function hence draws a random number between 0 and 1 and gives the
 * corresponding character using the bijection of the repartition function.
 *
 * All derived classes must initialize the _repartition and _size fields.
 */
class AbstractMutationProcess: public MutationProcess {
	
	protected:
		
		/**
		 * @brief The substitution model to use:
		 */
		const SubstitutionModel * _model;
	
		/**
		 * @brief The number of states allowed for the character to mutate.
		 */
    int _size;
	
		/**
		 * @brief The repartition function for states probabilities.
		 *
		 * _repartition[i][j] = probability that, being in state i at time t,
		 * we'll be in state <= j at time t+1.
		 */
    VVdouble _repartition;
	
	public:
		AbstractMutationProcess(const SubstitutionModel * model);
		virtual ~AbstractMutationProcess();
	
	public:
    int mutate(int state) const;
    int mutate(int state, int n) const;
		double getTimeBeforeNextMutationEvent(int state) const;
		int evolve(int initialState, double time) const;
		MutationPath detailedEvolve(int initialState, double time) const;
		const SubstitutionModel * getSubstitutionModel() const;
};

/**
 * @brief Generally used mutation process model.
 *
 * This builds a MutationProcess according to a given SubstitutionModel.
 * The underlying mutation process is the following:
 * <ol>
 * <li>Draw a random time @f$ t @f$ from an exponential law with parameter
 * @f$ - \lambda_i @f$,</li>
 * <li> Mutate the initial state. The probability of mutating state @f$ i @f$ 
 * to state @f$ j @f$ is:
 * @f[ \frac{Q_{i,j}}{\sum_k Q_{i,k}}. @f]</li>
 * </ol>
 */
class SimpleMutationProcess : public AbstractMutationProcess {
	
	public: // Constructor and destructor:
		
		/**
		 * @brief Build a new SimpleMutationProcess object.
		 *
		 * @param model The substitution model to use.
		 */
  	SimpleMutationProcess(const SubstitutionModel * model);
	
		~SimpleMutationProcess();

    /**
     * @brief Method redefinition for better performance.
     *
		 * @param initialState The state before beginning evolution.
		 * @param time         The time during which evolution must occure.
		 * @return The resulting state after evolution is completed.
		 */
		int evolve(int initialState, double time) const;
};

/**
 * This class is for testing purpose.
 * It allow "self" mutation of the kind i->i;
 */
class SelfMutationProcess : public AbstractMutationProcess {

  	public:
  		SelfMutationProcess(int alphabetSize);
	
			~SelfMutationProcess();
};

#endif	//_MUTATIONPROCESS_H_
