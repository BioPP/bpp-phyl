//
// File: RE08.cpp
// Created by: Julien Dutheil
// Created on: Mon Dec 29 10:15 2008
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

#include "RE08.h"

using namespace bpp;

#include <cmath>

using namespace std;

/******************************************************************************/

RE08::RE08(ReversibleSubstitutionModel *simpleModel, double lambda, double mu):
  AbstractReversibleSubstitutionModel(simpleModel->getAlphabet()),
  _simpleModel(simpleModel), _lambda(lambda), _mu(mu)
{
  _parameters.addParameter(Parameter("lambda", lambda, &Parameter::R_PLUS));  
  _parameters.addParameter(Parameter("mu", mu, &Parameter::R_PLUS));
  _parameters.addParameters(simpleModel->getParameters());
  //We need to overrired this from the AbstractSubstitutionModel constructor,
  //since the number of states in the model is no longer equal to the size of the alphabet.
  _size = simpleModel->getNumberOfStates() + 1;
  _generator.resize(_size, _size);
  _exchangeability.resize(_size, _size);
  _freq.resize(_size);
  _eigenValues.resize(_size);
  _leftEigenVectors.resize(_size, _size);
  _rightEigenVectors.resize(_size, _size);
  _p.resize(_size,_size);
  updateMatrices();
}

/******************************************************************************/
  
void RE08::updateMatrices()
{
  double f = (_lambda == 0 && _mu == 0) ? 1 : _lambda / (_lambda + _mu);
  
  // Frequencies:
  for(unsigned int i = 0; i < _size - 1; i++)
    _freq[i] = _simpleModel->freq(i) * f;

  _freq[_size-1] = (1. - f);

  _simpleGenerator = _simpleModel->getGenerator();
  _simpleExchangeabilities = _simpleModel->getExchangeabilityMatrix();

  // Generator and exchangeabilities:
  for(int i = 0; i < (int)_size - 1; i++)
  {
    for(int j = 0; j < (int)_size - 1; j++)
    {
      _generator(i, j) = _simpleGenerator(i, j);
      _exchangeability(i, j) = _simpleExchangeabilities(i, j) / f;
      if(i == j) 
      {
        _generator(i, j) -= _mu;
        _exchangeability(i, j) -= (_mu / f) / _simpleModel->freq(i);
      }
    }
    _generator(i, _size - 1) = _mu;
    _generator(_size - 1, i) = _lambda * _simpleModel->freq(i);
    _exchangeability(i, _size - 1) = _lambda + _mu;
    _exchangeability(_size - 1, i) = _lambda + _mu;
  }
  _generator(_size - 1, _size - 1) = -_lambda;
  _exchangeability(_size - 1, _size - 1) = -(_lambda + _mu); 

  //It is very likely that we are able to compute the eigen values and vector from the one of the simple model.
  //For now however, we will use a numerical diagonalization:
  AbstractSubstitutionModel::updateMatrices();
  //We do not use the one from  AbstractReversibleSubstitutionModel, since we already computed the generator.
}
  
/******************************************************************************/

double RE08::Pij_t(int i, int j, double d) const
{
  double f = (_lambda == 0 && _mu == 0) ? 1. : _lambda / (_lambda + _mu);
  if(i < (int)_size - 1 && j < (int)_size - 1)
  {
    return (_simpleModel->Pij_t(i, j, d) - _simpleModel->freq(j)) * exp(-_mu * d)
      + _freq[j] + (_simpleModel->freq(j) - _freq[j]) * exp(-(_lambda + _mu) * d);
  }
  else
  {
    if(i == (int)_size - 1)
    {
      if(j < (int)_size - 1)
      {
        return _freq[j] * (1. - exp(-(_lambda + _mu) * d));
      }
      else
      {
        return 1. - f * (1. - exp(-(_lambda + _mu) * d));
      }
    }
    else
    {  
      return _freq[j] * (1. - exp(-(_lambda + _mu) * d));
    }
  }
}

/******************************************************************************/

double RE08::dPij_dt(int i, int j, double d) const
{
  double f = (_lambda == 0 && _mu == 0) ? 1. : _lambda / (_lambda + _mu);
  if(i < (int)_size - 1 && j < (int)_size - 1)
  {
    return _simpleModel->dPij_dt(i, j, d) * exp(-_mu * d)
      - _mu * (_simpleModel->Pij_t(i, j, d) - _simpleModel->freq(j)) * exp(-_mu * d)
      - (_lambda + _mu) * (_simpleModel->freq(j) - _freq[j]) * exp(-(_lambda + _mu) * d);
  }
  else
  {
    if(i == (int)_size - 1)
    {
      if(j < (int)_size - 1)
      {
        return (_lambda + _mu) * _freq[j] * exp(-(_lambda + _mu) * d);
      }
      else
      {
        return - f * (_lambda + _mu) * exp(-(_lambda + _mu) * d);
      }
    }
    else
    {  
      return (_lambda + _mu) * _freq[j] * exp(-(_lambda + _mu) * d);
    }
  }
}

/******************************************************************************/

double RE08::d2Pij_dt2(int i, int j, double d) const
{
  double f = (_lambda == 0 && _mu == 0) ? 1. : _lambda / (_lambda + _mu);
  if(i < (int)_size - 1 && j < (int)_size - 1)
  {
    return _simpleModel->d2Pij_dt2(i, j, d) * exp(-_mu * d)
      - 2 * _mu * _simpleModel->dPij_dt(i, j, d) * exp(-_mu * d)
      + _mu * _mu * (_simpleModel->Pij_t(i, j, d) - _simpleModel->freq(j)) * exp(-_mu * d)
      + (_lambda + _mu) * (_lambda + _mu) * (_simpleModel->freq(j) - _freq[j]) * exp(-(_lambda + _mu) * d);
  }
  else
  {
    if(i == (int)_size - 1)
    {
      if(j < (int)_size - 1)
      {
        return - (_lambda + _mu) * (_lambda + _mu) * _freq[j] * exp(-(_lambda + _mu) * d);
      }
      else
      {
        return f * (_lambda + _mu) * (_lambda + _mu) * exp(-(_lambda + _mu) * d);
      }
    }
    else
    {  
      return - (_lambda + _mu) * (_lambda + _mu) * _freq[j] * exp(-(_lambda + _mu) * d);
    }
  }
}

/******************************************************************************/

RowMatrix<double> RE08::getPij_t(double d) const
{
  RowMatrix<double> simpleP = _simpleModel->getPij_t(d);
  double f = (_lambda == 0 && _mu == 0) ? 1. : _lambda / (_lambda + _mu);
  for(int i = 0; i < (int)_size - 1; i++)
  {
    for(int j = 0; j < (int)_size - 1; j++)
    {
      _p(i, j) = (simpleP(i, j) - _simpleModel->freq(j)) * exp(-_mu * d)
          + _freq[j] + (_simpleModel->freq(j) - _freq[j]) * exp(-(_lambda + _mu) * d);
    }
  }
  for(unsigned int j = 0; j < _size - 1; j++)
  {
    _p(_size - 1, j) = _freq[j] * (1. - exp(-(_lambda + _mu) * d));
  }
  _p(_size - 1, _size - 1) = 1. - f * (1. - exp(-(_lambda + _mu) * d));
  for(unsigned int i = 0; i < _size - 1; i++)
  {  
    _p(i, _size - 1) = _freq[_size - 1] * (1. - exp(-(_lambda + _mu) * d));
  }
  return _p;
}

/******************************************************************************/

RowMatrix<double> RE08::getdPij_dt(double d) const
{
  RowMatrix<double> simpleP = _simpleModel->getPij_t(d);
  RowMatrix<double> simpleDP = _simpleModel->getdPij_dt(d);
  double f = (_lambda == 0 && _mu == 0) ? 1. : _lambda / (_lambda + _mu);
  for(int i = 0; i < (int)_size - 1; i++)
  {
    for(int j = 0; j < (int)_size - 1; j++)
    {
      _p(i, j) = simpleDP(i, j) * exp(-_mu * d)
          - _mu * (simpleP(i, j) - _simpleModel->freq(j)) * exp(-_mu * d)
          - (_lambda + _mu) * (_simpleModel->freq(j) - _freq[j]) * exp(-(_lambda + _mu) * d);
    }
  }
  for(unsigned int j = 0; j < _size - 1; j++)
  {
    _p(_size - 1, j) = (_lambda + _mu) * _freq[j] * exp(-(_lambda + _mu) * d);
  }
  _p(_size - 1, _size - 1) = - f * (_lambda + _mu) * exp(-(_lambda + _mu) * d);
  for(unsigned int i = 0; i < _size - 1; i++)
  {  
    _p(i, _size - 1) = (_lambda + _mu) * _freq[_size - 1] * exp(-(_lambda + _mu) * d);
  }
  return _p;
}

/******************************************************************************/

RowMatrix<double> RE08::getd2Pij_dt2(double d) const
{
  RowMatrix<double> simpleP = _simpleModel->getPij_t(d);
  RowMatrix<double> simpleDP = _simpleModel->getdPij_dt(d);
  RowMatrix<double> simpleD2P = _simpleModel->getd2Pij_dt2(d);
  double f = (_lambda == 0 && _mu == 0) ? 1. : _lambda / (_lambda + _mu);
  for(int i = 0; i < (int)_size - 1; i++)
  {
    for(int j = 0; j < (int)_size - 1; j++)
    {
      _p(i, j) = simpleD2P(i, j) * exp(-_mu * d)
          - 2 * _mu * simpleDP(i, j) * exp(-_mu * d)
          + _mu * _mu * (simpleP(i, j) - _simpleModel->freq(j)) * exp(-_mu * d)
          + (_lambda + _mu) * (_lambda + _mu) * (_simpleModel->freq(j) - _freq[j]) * exp(-(_lambda + _mu) * d);
    }
  }
  for(unsigned int j = 0; j < _size - 1; j++)
  {
    _p(_size - 1, j) = - (_lambda + _mu) * (_lambda + _mu) * _freq[j] * exp(-(_lambda + _mu) * d);
  }
  _p(_size - 1, _size - 1) = f * (_lambda + _mu) * (_lambda + _mu) * exp(-(_lambda + _mu) * d);
  for(unsigned int i = 0; i < _size - 1; i++)
  {  
    _p(i, _size - 1) = - (_lambda + _mu) * (_lambda + _mu) * _freq[_size - 1] * exp(-(_lambda + _mu) * d);
  }
  return _p;
}

/******************************************************************************/

double RE08::getInitValue(int i, int state) const throw (BadIntException)
{
  if(i < 0 || i >= (int)_size) throw BadIntException(i, "RE08::getInitValue");
  if(state < -1 || !_alphabet->isIntInAlphabet(state)) throw BadIntException(state, "RE08::getInitValue. Character " + _alphabet->intToChar(state) + " is not allowed in model.");
  if(i == (int)_size - 1 && state == -1) return 1.;
  vector<int> states = _alphabet->getAlias(state);
  for(unsigned int j = 0; j < states.size(); j++) if(i == states[j]) return 1.;
  return 0.;
}

/******************************************************************************/

