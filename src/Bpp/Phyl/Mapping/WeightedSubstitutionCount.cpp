//
// File: WeightedSubstitutionCount.cpp
// Created by: Julien Dutheil
// Created on: Mon Dec 6 21:35 2010
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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

#include "WeightedSubstitutionCount.h"

using namespace bpp;
      
double WeightedSubstitutionCount::getNumberOfSubstitutions(unsigned int initialState, unsigned int finalState, double length, unsigned int type) const
{
  double weightedCount = 0;
  for (unsigned int x = 0; x < getAlphabet()->getSize(); ++x) {
    for (unsigned int y = 0; y < getAlphabet()->getSize(); ++y) {
      if (register_->getType(x, y) == type) {
        weightedCount += subCount_->getNumberOfSubstitutions(initialState, finalState, length, subCount_->getSubstitutionRegister()->getType(x, y)) * dist_->getIndex(x, y);
      }
    }
  }
  return weightedCount;
}

Matrix<double>* WeightedSubstitutionCount::getAllNumbersOfSubstitutions(double length, unsigned int type) const
{
  Matrix<double>* mat = new RowMatrix<double>(getNumberOfStates(), getNumberOfStates());
  for (unsigned int i = 0; i < mat->getNumberOfRows(); ++i)
    for (unsigned int j = 0; j < mat->getNumberOfColumns(); ++j)
      (*mat)(i, j) = getNumberOfSubstitutions(i, j, length, type);
  return mat;
}

std::vector<double> WeightedSubstitutionCount::getNumberOfSubstitutionsForEachType(unsigned int initialState, unsigned int finalState, double length) const
{
  std::vector<double> v(register_->getNumberOfSubstitutionTypes());
  for (unsigned int t = 1; t <= register_->getNumberOfSubstitutionTypes(); ++t)
    v[t - 1] = getNumberOfSubstitutions(initialState, finalState, length, t);
  return v;
}
    
