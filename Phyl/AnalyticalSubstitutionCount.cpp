//
// File: AnalyticalSubstitutionCount.cpp
// Created by: Julien Dutheil
// Created on: Wed Apr 5 11:21 2006
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004, 2005, 2006)

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

#include "AnalyticalSubstitutionCount.h"

//From NumCalc:
#include <NumCalc/MatrixTools.h>

AnalyticalSubstitutionCount::AnalyticalSubstitutionCount(const SubstitutionModel * model, int cutOff):
	_model(model),
	_cuttOff(cutOff),
	_currentLength(-1.)
{
	unsigned int n = model -> getAlphabet() -> getSize();
	M = RowMatrix<double>(n, n);
	MatrixTools::fill(M, 0.);
};

/******************************************************************************/

double AnalyticalSubstitutionCount::getNumberOfSubstitutions(int initialState, int finalState, double length) const
{
  if(length == _currentLength) return M(initialState, finalState);
  if(length < 0.000001) return initialState == finalState ? 0. : 1.; //Limit case!
  // Else we need to recompute M:
  RowMatrix<double> Q = _model -> getGenerator();
  // L is the diagonal matrix with all substitution rates.
	unsigned int s = Q.nRows();
	RowMatrix<double> QL(s, s);
	for(unsigned int i = 0; i < s; i++) {
	 for(unsigned int j = 0; j < s; j++) {
	  QL(i, j) = ((i == j) ? 0. : Q(i, j)) ;
	 }
	}

	MatrixTools::fill(M, 0.);
	for(int n = 1; n < _cuttOff; n++) {
	  RowMatrix<double> M2(s, s);
	  MatrixTools::fill(M2, 0.);
	  for(int p = 0; p < n; p++) {
	    MatrixTools::add(M2, MatrixTools::mult(MatrixTools::mult(MatrixTools::pow(Q, p), QL), MatrixTools::pow(Q, n - p - 1)));
		}
		MatrixTools::scale(M2, pow(length, n) / NumTools::fact(n));
		MatrixTools::add(M, M2);
	}

	// Now we must divide by pijt:
	for(unsigned int i = 0; i < s; i++) {
	  for(unsigned int j = 0; j < s; j++) {
		  M(i, j) /= _model -> Pij_t(i, j, length);
	  }
	}

	_currentLength = length;
  return M(initialState, finalState);
}

/******************************************************************************/

Matrix<double> * AnalyticalSubstitutionCount::getAllNumbersOfSubstitutions(double length) const
{
  if(length == _currentLength) return new RowMatrix<double>(M);
	if(length < 0.000001) // Limit case!
  { 
    unsigned int s = _model->getAlphabet()->getSize();
	  for(unsigned int i = 0; i < s; i++) {
		  for(unsigned int j = 0; j < s; j++) {
			  M(i, j) = i == j ? 0. : 1.;
		  }
	  }
  } else {
	  // Else we need to recompute M:
	  RowMatrix<double> Q = _model -> getGenerator();
	  // L is the diagonal matrix with all substitution rates.
	  unsigned int s = Q.nRows();
	  RowMatrix<double> QL(s, s);
	  for(unsigned int i = 0; i < s; i++) {
		  for(unsigned int j = 0; j < s; j++) {
			  QL(i, j) = ((i == j) ? 0. : Q(i, j)) ;
		  }
	  }

	  MatrixTools::fill(M, 0.);
	  for(int n = 1; n < _cuttOff; n++) {
		  RowMatrix<double> M2(s, s);
		  MatrixTools::fill(M2, 0.);
		  for(int p = 0; p < n; p++) {
			  MatrixTools::add(M2, MatrixTools::mult(MatrixTools::mult(MatrixTools::pow(Q, p), QL), MatrixTools::pow(Q, n - p - 1)));
		  }
		  MatrixTools::scale(M2, pow(length, n) / NumTools::fact(n));
		  MatrixTools::add(M, M2);
	  }

	  // Now we must divide by pijt:
	  for(unsigned int i = 0; i < s; i++) {
		  for(unsigned int j = 0; j < s; j++) {
			  M(i, j) /= _model -> Pij_t(i, j, length);
		  }
	  }
  }

	_currentLength = length;

	return new RowMatrix<double>(M);
}

/******************************************************************************/

