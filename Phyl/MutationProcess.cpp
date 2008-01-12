//
// File: MutationProcess.cpp
// Created by: Julien Dutheil
// Created on: Wed Mar 12 16:11:44 2003
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
 
#include "MutationProcess.h"

// From NumCalc:
#include <NumCalc/RandomTools.h>

// From Utils:
#include <Utils/TextTools.h>

using namespace bpp;

/******************************************************************************/

AbstractMutationProcess::AbstractMutationProcess(const SubstitutionModel * model): _model(model) {}
	
AbstractMutationProcess::~AbstractMutationProcess() {}
	
/******************************************************************************/
	
const SubstitutionModel * AbstractMutationProcess::getSubstitutionModel() const
{
	return _model;
}
	
/******************************************************************************/

int AbstractMutationProcess::mutate(int state) const
{
  double alea = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
  for(unsigned int j = 0; j < _size; j++)
  { 
   	if(alea < _repartition[state][j]) return j;
	}
  return _size;
}

/******************************************************************************/

int AbstractMutationProcess::mutate(int state, unsigned int n) const
{
 	int s = state;
 	for(unsigned int k = 0; k < n; k++)
  {
   	double alea = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
   	for(unsigned int j = 1; j < _size + 1; j++)
    { 
   		if(alea < _repartition[s][j])
      {
       	s = j;
     		break;
  		}
   	}
  }
  return s;
}

/******************************************************************************/

double AbstractMutationProcess::getTimeBeforeNextMutationEvent(int state) const
{
	return RandomTools::randExponential(- _model->Qij(state, state));
}

/******************************************************************************/

int AbstractMutationProcess::evolve(int initialState, double time) const
{
	double t = 0;
	int currentState = initialState;
  t += getTimeBeforeNextMutationEvent(currentState);
	while(t < time)
  {
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
	while(t < time)
  {
   	currentState = mutate(currentState);
		mp.addEvent(currentState, t);
		t += getTimeBeforeNextMutationEvent(currentState);
  }
  return mp;
}

/******************************************************************************/

SimpleMutationProcess::SimpleMutationProcess(const SubstitutionModel * model):
  AbstractMutationProcess(model)
{
	_size = model->getNumberOfStates();
  _repartition = VVdouble(_size);
  // Each element contains the probabilities concerning each character in the alphabet.

  // We will now initiate each of these probability vector.
	RowMatrix<double> Q = model->getGenerator();
  for(unsigned int i = 0; i < _size; i++)
  {
   	_repartition[i] = Vdouble(_size);
   	double cum = 0;
   	double sum_Q = 0;
   	for(unsigned int j = 0; j < _size; j++)
    {
   		if(j != i) sum_Q += Q(i, j);
   	}
   	for(unsigned int j = 0; j < _size; j++)
    {
    	if(j != i)
      {
 		   	cum += model->Qij(i, j) / sum_Q;
     		_repartition[i][j] = cum;
     	}
      else _repartition[i][j] = -1; // Forbiden value: does not correspond to a change.
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
	Vdouble pijt(_size);
	pijt[0] = _model->Pij_t(initialState, 0, time);
	for(unsigned int i = 1; i < _size; i++)
  {
		pijt[i] = pijt[i - 1] + _model->Pij_t(initialState, i, time);
	}
	double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
	for(unsigned int i = 0; i < _size; i++)
  {
		if(rand < pijt[i]) return i;
	}
	throw Exception("SimpleSimulationProcess::evolve(intialState, time): error all pijt do not sum to one (total sum = " + TextTools::toString(pijt[_size - 1]) + ").");
}

/******************************************************************************/

SelfMutationProcess::SelfMutationProcess(int alphabetSize):
AbstractMutationProcess(NULL)
{
  	_size = alphabetSize;
  	_repartition = VVdouble(_size);
  	// Each element contains the probabilities concerning each character in the alphabet.

  	// We will now initiate each of these probability vector.
  	for(unsigned int i = 0; i < _size; i++)
    {
	    _repartition[i] = Vdouble(_size);
    	for(unsigned int j = 0; j < _size; j++)
      {
      	_repartition[i][j] = (j+1.0)/_size;
    	}
  	}
  	// Note that I use cumulative probabilities in _repartition (hence the name).
  	// These cumulative probabilities are useful for the 'mutate(...)' function.
}

SelfMutationProcess::~SelfMutationProcess() {}
	
/******************************************************************************/

