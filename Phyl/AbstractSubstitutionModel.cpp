//
// File: AbstractSubstitutionModel.cpp
// Created by: Julien Dutheil
// Created on: Tue May 27 10:31:49 2003
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

#include "AbstractSubstitutionModel.h"

// From Utils:
#include <Utils/TextTools.h>

// From NumCalc:
#include <NumCalc/MatrixTools.h>
#include <NumCalc/VectorTools.h>
using namespace VectorOperators;
#include <NumCalc/EigenValue.h>

// From SeqLib:
#include <Seq/SequenceContainerTools.h>

/******************************************************************************/

AbstractSubstitutionModel::AbstractSubstitutionModel(const Alphabet * alpha): _alphabet(alpha)
{
	_size = alpha->getSize();
	_generator.resize(_size, _size);
	_freq.resize(_size);
	_eigenValues.resize(_size);
	_leftEigenVectors.resize(_size, _size);
	_rightEigenVectors.resize(_size, _size);
}

/******************************************************************************/
	
void AbstractSubstitutionModel::updateMatrices()
{
	// Compute eigen values and vectors:
	EigenValue<double> ev(_generator);
	_rightEigenVectors = ev.getV();
	_leftEigenVectors = MatrixTools::inv(_rightEigenVectors);
	_eigenValues = ev.getRealEigenValues();
}

/******************************************************************************/

RowMatrix<double> AbstractSubstitutionModel::getPij_t(double t) const
{
	if(t == 0) return MatrixTools::getId< RowMatrix<double> >(_size);
  return MatrixTools::mult(_rightEigenVectors, VectorTools::exp(_eigenValues * t), _leftEigenVectors);
}

RowMatrix<double> AbstractSubstitutionModel::getdPij_dt(double t) const
{
	return MatrixTools::mult(_rightEigenVectors, _eigenValues * VectorTools::exp(_eigenValues * t), _leftEigenVectors);
}

RowMatrix<double> AbstractSubstitutionModel::getd2Pij_dt2(double t) const
{
	return MatrixTools::mult(_rightEigenVectors, NumTools::sqr(_eigenValues) * VectorTools::exp(_eigenValues * t), _leftEigenVectors);
}

/******************************************************************************/

double AbstractSubstitutionModel::getInitValue(int i, int state) const throw (BadIntException)
{
	if(i < 0 || i > (int)_alphabet->getSize()) throw BadIntException(i, "AbstractSubstitutionModel::getInitValue");
	//Old method: do not care about generic characters:
	//if(state < 0 || state > (int)alphabet->getSize()) throw BadIntException(state, "AbstractSubstitutionModel::getInitValue. Character " + alphabet->intToChar(state) + " is not allowed in model.");
	//return i == state ? 1 : 0;
	if(state < 0 || !_alphabet->isIntInAlphabet(state)) throw BadIntException(state, "AbstractSubstitutionModel::getInitValue. Character " + _alphabet->intToChar(state) + " is not allowed in model.");
	vector<int> states = _alphabet->getAlias(state);
	for(unsigned int j = 0; j < states.size(); j++) if(i == states[j]) return 1.;
	return 0.;
}

/******************************************************************************/

void AbstractSubstitutionModel::setFreqFromData(const SequenceContainer & data)
{
	map<int, double> freqs = SequenceContainerTools::getFrequencies(data);
	double t = 0;
	for(unsigned int i = 0; i < _size; i++) t += freqs[i];
	for(unsigned int i = 0; i < _size; i++) _freq[i] = freqs[i] / t;
	//Re-compute generator and eigen values:
	updateMatrices();
}

/******************************************************************************/

double AbstractSubstitutionModel::getScale() const
{
	return -VectorTools::scalar<double, double>(MatrixTools::diag<RowMatrix<double>, double>(_generator), _freq);
}

/******************************************************************************/

AbstractReversibleSubstitutionModel::AbstractReversibleSubstitutionModel(const Alphabet * alpha):
  AbstractSubstitutionModel(alpha)
{
	_exchangeability.resize(_size, _size);
}

/******************************************************************************/
		
void AbstractReversibleSubstitutionModel::updateMatrices()
{
	RowMatrix<double> Pi = MatrixTools::diag<RowMatrix<double>, double>(_freq);
	_generator = MatrixTools::mult(_exchangeability, Pi); //Diagonal elements of the exchangability matrix will be ignored.
	// Compute diagonal elements of the generator:
	for(unsigned int i = 0; i < _size; i++)
  {
		double lambda = 0;
		for(unsigned int j = 0; j < _size; j++)
    {
			if(j!=i) lambda += _generator(i,j);
		}
		_generator(i,i) = -lambda;
	}
	// Normalization:
	double scale = getScale();
	MatrixTools::scale(_generator, 1./scale);

  // Normalize exchangeability matrix too:
	MatrixTools::scale(_exchangeability, 1./scale);
  // Compute diagonal elements of the exchangeability matrix:
	for(unsigned int i = 0; i < _size; i++)
  {
		_exchangeability(i,i) = _generator(i,i) / _freq[i];
	}
  AbstractSubstitutionModel::updateMatrices();
}

/******************************************************************************/


