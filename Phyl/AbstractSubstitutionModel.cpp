//
// File: AbstractSubstitutionModel.h
// Created by:  <@bogdanof>
// Created on: Tue May 27 10:31:49 2003
//

#include "AbstractSubstitutionModel.h"

// From Utils:
#include <Utils/TextTools.h>

// From NumCalc:
#include <NumCalc/MatrixTools.h>
using namespace MatrixTools;
#include <NumCalc/VectorTools.h>
using namespace VectorFunctions;
using namespace VectorOperators;

// From SeqLib:
#include <Seq/SequenceContainerTools.h>

AbstractSubstitutionModel::AbstractSubstitutionModel(const Alphabet * alpha): alphabet(alpha)
{
	_size = alpha -> getSize();
	_generator         = Matrix(_size, _size);
	_exchangeability   = Matrix(_size, _size);
	_freq              = Vector(_size);
	_eigenValues       = Vector(_size);
	_leftEigenVectors  = Matrix(_size, _size);
	_rightEigenVectors = Matrix(_size, _size);
}

/******************************************************************************/
	
const Alphabet * AbstractSubstitutionModel::getAlphabet() const { return alphabet; }

/******************************************************************************/
	
Vector AbstractSubstitutionModel::getFrequencies() const { return _freq; }

Matrix AbstractSubstitutionModel::getExchangeabilityMatrix() const { return _exchangeability; }

Matrix AbstractSubstitutionModel::getGenerator() const { return _generator; }

Vector AbstractSubstitutionModel::eigenValues() const { return _eigenValues; }

Matrix AbstractSubstitutionModel::horizontalLeftEigenVectors() const { return _leftEigenVectors; }

Matrix AbstractSubstitutionModel::verticalRightEigenVectors() const { return _rightEigenVectors; }

Matrix AbstractSubstitutionModel::getPij_t(double t) const
{
	return _rightEigenVectors * diag<Matrix, double>(exp(_eigenValues)) * _leftEigenVectors;
}

Matrix AbstractSubstitutionModel::getdPij_dt(double t) const
{
	return _rightEigenVectors * diag<Matrix, double>(_eigenValues * exp(_eigenValues)) * _leftEigenVectors;
}

Matrix AbstractSubstitutionModel::getd2Pij_dt2(double t) const
{
	return _rightEigenVectors * diag<Matrix, double>(sqr(_eigenValues) * exp(_eigenValues)) * _leftEigenVectors;
}

double AbstractSubstitutionModel::freq(int i) const { return _freq[i]; }

double AbstractSubstitutionModel::Qij(int i, int j) const { return _generator(i, j); }

double AbstractSubstitutionModel::Pij_t(int i, int j, double t) const
{
	return scalar(row<Matrix, Vector>(_rightEigenVectors, i) * exp(_eigenValues), col<Matrix, Vector>(_leftEigenVectors, j));
}

double AbstractSubstitutionModel::dPij_dt(int i, int j, double t) const
{
	return scalar(row<Matrix, Vector>(_rightEigenVectors, i) * _eigenValues * exp(_eigenValues), col<Matrix, Vector>(_leftEigenVectors, j));
}

double AbstractSubstitutionModel::d2Pij_dt2(int i, int j, double t) const
{
	return scalar(row<Matrix, Vector>(_rightEigenVectors, i) * sqr(_eigenValues) *exp(_eigenValues), col<Matrix, Vector>(_leftEigenVectors, j));
}
	

/******************************************************************************/

double AbstractSubstitutionModel::getInitValue(int i, int state) const throw (BadIntException) {
	if(i     < 0 || i     > (int)alphabet -> getSize()) throw BadIntException(i    , "AbstractSubstitutionModel::getInitValue");
	//Old method: do not care about generic characters:
	//if(state < 0 || state > (int)alphabet -> getSize()) throw BadIntException(state, "AbstractSubstitutionModel::getInitValue. Character " + alphabet -> intToChar(state) + " is not allowed in model.");
	//return i == state ? 1 : 0;
	if(state < 0 || !alphabet -> isIntInAlphabet(state)) throw BadIntException(state, "AbstractSubstitutionModel::getInitValue. Character " + alphabet -> intToChar(state) + " is not allowed in model.");
	vector<int> states = alphabet -> getAlias(state);
	for(unsigned int j = 0; j < states.size(); j++) if(i == states[j]) return 1;
	return 0;
}

/******************************************************************************/

void AbstractSubstitutionModel::setFreqFromData(const SequenceContainer & data) {
	map<int, double> freqs = SequenceContainerTools::getFrequencies(data);
	double t = 0;
	for(unsigned int i = 0; i < _size; i++) t += freqs[i];
	for(unsigned int i = 0; i < _size; i++) _freq[i] = freqs[i] / t;
	//Re-compute generator and eigen values:
	updateMatrices();
}

/******************************************************************************/

ParameterList AbstractSubstitutionModel::getParameters() const {
	return _parameters;
}

/******************************************************************************/

double AbstractSubstitutionModel::getParameter(const string & name) const
throw (ParameterNotFoundException)
{
	Parameter * p = _parameters.getParameter(name);
	if(p == NULL) throw ParameterNotFoundException("AbstractSubstitutionModel::getParameter", name);
	return p -> getValue();
}

/******************************************************************************/

void AbstractSubstitutionModel::setAllParametersValues(const ParameterList & params)
throw (ParameterNotFoundException, ConstraintException)
{
	_parameters.setAllParametersValues(params);
	updateMatrices();
}

/******************************************************************************/

void AbstractSubstitutionModel::setParameterValue(const string & name, double value)
throw (ParameterNotFoundException, ConstraintException)
{
	_parameters.setParameterValue(name, value);
	updateMatrices();
}

/******************************************************************************/

void AbstractSubstitutionModel::setParametersValues(const ParameterList & params)
throw (ParameterNotFoundException, ConstraintException)
{ 
	_parameters.setParametersValues(params);
	updateMatrices();
}

/******************************************************************************/

void AbstractSubstitutionModel::matchParametersValues(const ParameterList & params)
throw (ConstraintException)
{
	_parameters.matchParametersValues(params);		
	updateMatrices();
}

/******************************************************************************/
