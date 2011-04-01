//
// File: UniformizationSubstitutionCount.cpp
// Created by: Julien Dutheil
// Created on: Sat Mar 19 13:54 2011
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

#include "UniformizationSubstitutionCount.h"

#include "Bpp/Numeric/Matrix/MatrixTools.h"
#include "Bpp/Numeric/NumTools.h"
#include <vector>

using namespace bpp;
using namespace std;

/******************************************************************************/

UniformizationSubstitutionCount::UniformizationSubstitutionCount(const SubstitutionModel* model, SubstitutionRegister* reg) :
  AbstractSubstitutionCount(reg), model_(model),
  nbStates_(model->getNumberOfStates()),
  bMatrices_(reg->getNumberOfSubstitutionTypes()),
  power_(),
  s_(reg->getNumberOfSubstitutionTypes()),
  miu_(0),
  counts_(reg->getNumberOfSubstitutionTypes()),
  currentLength_(-1.)
{
  //Check compatiblity between model and substitution register:
  if (model->getAlphabet()->getAlphabetType() != reg->getAlphabet()->getAlphabetType())
    throw Exception("UniformizationSubstitutionCount (constructor): alphabets do not match between register and model.");

  //Initialize all B matrices according to substitution register. This is done once for all,
  //unless the number of states changes:
  for (unsigned int i = 0; i < reg->getNumberOfSubstitutionTypes(); ++i) {
    bMatrices_[i].resize(nbStates_, nbStates_);
    counts_[i].resize(nbStates_, nbStates_);
  }
  for (unsigned int j = 0; j < nbStates_; ++j) {
    for (unsigned int k = 0; k < nbStates_; ++k) {
      unsigned int i = reg->getType(static_cast<int>(j), static_cast<int>(k));
      if (i > 0) {
        bMatrices_[i - 1](j, k) = model->Qij(j, k);
      }
    }
  }

  for (unsigned int i = 0; i < nbStates_; ++i) {
    double diagQ = abs(model_->Qij(i, i));
    if (diagQ > miu_)
      miu_ = diagQ;
  }
}				
    
/******************************************************************************/

void UniformizationSubstitutionCount::computeCounts_(double length) const
{
	double lam = miu_ * length;
	RowMatrix<double> I;
  MatrixTools::getId(nbStates_, I);
  RowMatrix<double> R(model_->getGenerator());
  MatrixTools::scale(R, 1. / miu_);
  MatrixTools::add(R, I);
	
	//compute the stopping point
	//use the tail of Poission distribution
	//can be approximated by 4 + 6 * sqrt(lam) + lam
	unsigned int nMax = ceil(4 + 6 * sqrt(lam) + lam);

	//compute the powers of R
	power_.resize(nMax + 1);
	power_[0] = I;
	for (unsigned int i = 1; i < nMax + 1; ++i)
		MatrixTools::mult(power_[i - 1], R, power_[i]);

  for (unsigned int i = 0; i < register_->getNumberOfSubstitutionTypes(); ++i) {
    s_[i].resize(nMax + 1);
    MatrixTools::mult(bMatrices_[i], power_[0], s_[i][0]);
    RowMatrix<double> tmp(nbStates_, nbStates_);
	  for (unsigned int l = 1; l < nMax + 1; ++l) {
      MatrixTools::mult(R, s_[i][l - 1], s_[i][l]);
      MatrixTools::mult(bMatrices_[i], power_[l], tmp);
      MatrixTools::add(s_[i][l], tmp);
    }
    MatrixTools::fill(counts_[i], 0);
	  for (unsigned int l = 0; l < nMax + 1; ++l) {
      tmp = s_[i][l];
      MatrixTools::scale(tmp, (pow(lam, l + 1) * exp(-lam) / NumTools::fact(l + 1))/miu_);
      MatrixTools::add(counts_[i], tmp);
    }
  }

  // Now we must divide by pijt:
  RowMatrix<double> P = model_->getPij_t(length);
  for (unsigned int i = 0; i < register_->getNumberOfSubstitutionTypes(); i++) {
    for (unsigned int j = 0; j < nbStates_; j++) {
      for(unsigned int k = 0; k < nbStates_; k++) {
        counts_[i](j, k) /= P(j, k);
        if (isnan(counts_[i](j, k)))
          counts_[i](j, k) = 0;
      }
    }
  }
}

/******************************************************************************/

Matrix<double>* UniformizationSubstitutionCount::getAllNumbersOfSubstitutions(double length, unsigned int type) const
{
  if (length < 0)
    throw Exception("UniformizationSubstitutionCount::getAllNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeCounts_(length);
    currentLength_ = length;
  }
  return new RowMatrix<double>(counts_[type - 1]);
}

/******************************************************************************/

double UniformizationSubstitutionCount::getNumberOfSubstitutions(unsigned int initialState, unsigned int finalState, double length, unsigned int type) const
{
  if (length < 0)
    throw Exception("UniformizationSubstitutionCount::getNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeCounts_(length);
    currentLength_ = length;
  }
  return counts_[type - 1](initialState, finalState);
}

/******************************************************************************/

std::vector<double> UniformizationSubstitutionCount::getNumberOfSubstitutionsForEachType(unsigned int initialState, unsigned int finalState, double length) const
{
  if (length < 0)
    throw Exception("UniformizationSubstitutionCount::getNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeCounts_(length);
    currentLength_ = length;
  }
  std::vector<double> v(getNumberOfSubstitutionTypes());
  for (unsigned int t = 0; t < getNumberOfSubstitutionTypes(); ++t) {
    v[t] = counts_[t](initialState, finalState);
  }
  return v;
}
    
/******************************************************************************/

void UniformizationSubstitutionCount::setSubstitutionModel(const SubstitutionModel* model)
{
  //Check compatiblity between model and substitution register:
  if (model->getAlphabet()->getAlphabetType() != register_->getAlphabet()->getAlphabetType())
    throw Exception("UniformizationSubstitutionCount::setSubstitutionModel: alphabets do not match between register and model.");
  model_ = model;
  unsigned int n = model->getAlphabet()->getSize();
  if (n != nbStates_) {
    nbStates_ = n;
    //Re-initialize all B matrices according to substitution register.
    for (unsigned int i = 0; i < register_->getNumberOfSubstitutionTypes(); ++i) {
      bMatrices_[i].resize(nbStates_, nbStates_);
      counts_[i].resize(nbStates_, nbStates_);
    }
  }
  for (unsigned int j = 0; j < nbStates_; ++j) {
    for (unsigned int k = 0; k < nbStates_; ++k) {
      unsigned int i = register_->getType(static_cast<int>(j), static_cast<int>(k));
      if (i > 0) {
        bMatrices_[i - 1](j, k) = abs(model->Qij(j, k));
      }
    }
  }
	
  miu_ = 0;
  for (unsigned int i = 0; i < nbStates_; ++i) {
    double diagQ = abs(model_->Qij(i, i));
    if (diagQ > miu_)
      miu_ = diagQ;
  }

  //Recompute counts:
  computeCounts_(currentLength_);
}

/******************************************************************************/

