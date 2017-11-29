//
// File: NaiveSubstitutionCount.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 11:08 2006
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

#include "NaiveSubstitutionCount.h"

using namespace bpp;

Matrix<double>* NaiveSubstitutionCount::getAllNumbersOfSubstitutions(double length, size_t type) const
{ 
  size_t n = supportedChars_.size();
  RowMatrix<double>* mat = new RowMatrix<double>(n, n);
  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = 0; j < n; ++j)
    {
      (*mat)(i, j) = (register_->getType(i, j) == type ? (weights_ ? weights_->getIndex(supportedChars_[i], supportedChars_[j]) : 1.) : 0.);
    }
  }
  return mat;
}

LabelSubstitutionCount::LabelSubstitutionCount(const SubstitutionModel* model) :
  AbstractSubstitutionCount(
    new TotalSubstitutionRegister(model->getStateMap())),
  label_(model->getNumberOfStates(), model->getNumberOfStates()),
  supportedChars_(model->getAlphabetStates())
{
  size_t n = supportedChars_.size();
  double count = 0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (i == j) label_(i, j) = 0;
      else label_(i, j) = ++count;
    }
  }
}			

