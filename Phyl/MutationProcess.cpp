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
 
#include "MutationProcess.h"

// From NumCalc:
#include <NumCalc/RandomTools.h>

// From Utils:
#include <Utils/TextTools.h>

/******************************************************************************/

AbstractMutationProcess::AbstractMutationProcess(const SubstitutionModel * model):
_model(model) {}
	
AbstractMutationProcess::~AbstractMutationProcess() {}
	
/******************************************************************************/
	
const SubstitutionModel * AbstractMutationProcess::getSubstitutionModel() const {
	return _model;
}
	
/******************************************************************************/

int AbstractMutationProcess::mutate(int state) const
{
  double alea = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
  for(int j = 0; j < _size; j++) { 
   	if(alea < _repartition[state][j]) return j;
	}
  return _size;
}

/******************************************************************************/

int AbstractMutationProcess::mutate(int state, int n) const
{
 	int s = state;
 	for(int k = 0; k < n; k++) {
   	double alea = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
   	for(int j = 1; j < _size + 1; j++) { 
   		if(alea < _repartition[s][j]) {
       	s = j;
     		break;
  		}
   	}
  }
  return s;
}

/******************************************************************************/

double AbstractMutationProcess::getTimeBeforeNextMutationEvent(int state) const {
	return RandomTools::randExponential(- _model -> Qij(state, state));
}

/******************************************************************************/

int AbstractMutationProcess::evolve(int initialState, double time) const
{
	double t = 0;
	int currentState = initialState;
  t += getTimeBeforeNextMutationEvent(currentState);
	while(t < time) {
   	currentState = mutate(currentState);
		t += getTimeBeforeNextMutationEvent(currentState);
  }
  return currentState;
}

/******************************************************************************/

MutationPath AbstractMutationProcess::detailedEvolve(int initialState, double time) const
{
	MutationPath mp(initialState, time);
	double t = 0;
	int currentState = initialState;
  t += getTimeBeforeNextMutationEvent(currentState);
	while(t < time) {
   	currentState = mutate(currentState);
		mp.addEvent(currentState, t);
		t += getTimeBeforeNextMutationEvent(currentState);
  }
  return mp;
}

/******************************************************************************/

SimpleMutationProcess::SimpleMutationProcess(const SubstitutionModel * model):
AbstractMutationProcess(model) {
	_size = model -> getAlphabet() -> getSize();
  _repartition = VVdouble(_size);
  // Each element contains the probabilities concerning each character in the alphabet.

  // We will now initiate each of these probability vector.
  for(int i = 0; i < _size; i++) {
   	_repartition[i] = Vdouble(_size);
   	double cum = 0;
   	double sum_Q = 0;
   	for(int j = 0; j < _size; j++) {
   		if(j != i) sum_Q += model -> Qij(i, j);
   	}
   	for(int j = 0; j < _size; j++) {
    	if(j != i) {
 		   	cum += model -> Qij(i, j) / sum_Q;
     		_repartition[i][j] = cum;
     	} else _repartition[i][j] = -1; // Forbiden value: does not correspond to a change.
    }
 	}
  // Note that I use cumulative probabilities in _repartition (hence the name).
  // These cumulative probabilities are useful for the 'mutate(...)' function.
}

SimpleMutationProcess::~SimpleMutationProcess() {}

/******************************************************************************/

int SimpleMutationProcess::evolve(int initialState, double time) const
{
	// Compute all cumulative pijt:
	unsigned int s = _model -> getAlphabet() -> getSize();
	Vdouble pijt(s);
	pijt[0] = _model -> Pij_t(initialState, 0, time);
	for(unsigned int i = 1; i < s; i++) {
		pijt[i] = pijt[i - 1] + _model -> Pij_t(initialState, i, time);
	}
	double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
	for(unsigned int i = 0; i < s; i++) {
		if(rand < pijt[i]) return i;
	}
	throw Exception("SimpleSimulationProcess::evolve(intialState, time): error all pijt do not sum to one (total sum = " + TextTools::toString(pijt[s - 1]) + ").");
}

/******************************************************************************/

SelfMutationProcess::SelfMutationProcess(int alphabetSize):
AbstractMutationProcess(NULL) {
  	_size = alphabetSize;
  	_repartition = VVdouble(_size);
  	// Each element contains the probabilities concerning each character in the alphabet.

  	// We will now initiate each of these probability vector.
  	for(int i = 0; i < _size; i++) {
	    _repartition[i] = Vdouble(_size);
    	for(int j = 0; j < _size; j++) {
      	_repartition[i][j] = (j+1.0)/_size;
    	}
  	}
  	// Note that I use cumulative probabilities in _repartition (hence the name).
  	// These cumulative probabilities are useful for the 'mutate(...)' function.
}

SelfMutationProcess::~SelfMutationProcess() {}
	
/******************************************************************************/
