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

UniformizationSubstitutionCount::UniformizationSubstitutionCount(const SubstitutionModel* model, SubstitutionRegister* reg, const AlphabetIndex2* weights) :
  AbstractSubstitutionCount(reg),
  AbstractWeightedSubstitutionCount(weights, true),
  model_(model),
  nbStates_(model->getNumberOfStates()),
  bMatrices_(reg->getNumberOfSubstitutionTypes()),
  power_(),
  s_(reg->getNumberOfSubstitutionTypes()),
  miu_(0),
  counts_(reg->getNumberOfSubstitutionTypes()),
  currentLength_(1.)
{
  //Check compatiblity between model and substitution register:
  if (model->getAlphabet()->getAlphabetType() != reg->getAlphabet()->getAlphabetType())
    throw Exception("UniformizationSubstitutionCount (constructor): alphabets do not match between register and model.");

  //Initialize all B matrices according to substitution register. This is done once for all,
  //unless the number of states changes:
  initBMatrices_();
  fillBMatrices_();

  for (unsigned int i = 0; i < nbStates_; ++i) {
    double diagQ = abs(model_->Qij(i, i));
    if (diagQ > miu_)
      miu_ = diagQ;
  }
  
  if (miu_>10000)
    throw Exception("UniformizationSubstitutionCount::UniformizationSubstitutionCount The maximum diagonal values of generator is above 10000. Abort, chose another mapping method");
      
}        

/******************************************************************************/

void UniformizationSubstitutionCount::resetBMatrices_()
{
  size_t nbTypes = register_->getNumberOfSubstitutionTypes();
  bMatrices_.resize(nbTypes);
  counts_.resize(nbTypes);
  s_.resize(nbTypes);
}


void UniformizationSubstitutionCount::initBMatrices_()
{
  //Re-initialize all B matrices according to substitution register.
  for (size_t i = 0; i < register_->getNumberOfSubstitutionTypes(); ++i) {
    bMatrices_[i].resize(nbStates_, nbStates_);
    counts_[i].resize(nbStates_, nbStates_);
  }
}

void UniformizationSubstitutionCount::fillBMatrices_()
{
  for (size_t j = 0; j < nbStates_; ++j) {
    for (size_t k = 0; k < nbStates_; ++k) {
      size_t i = register_->getType(j, k);
      if (i > 0 && k != j) {
        //jdutheil on 25/07/14: I think this is incorrect, weights should only come at the end.
        //bMatrices_[i - 1](j, k) = model_->Qij(j, k) * (weights_ ? weights_->getIndex(fromState, toState) : 1);
        bMatrices_[i - 1](j, k) = model_->Qij(j, k);
      }
    }
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
  //use the tail of Poisson distribution
  //can be approximated by 4 + 6 * sqrt(lam) + lam
  size_t nMax = static_cast<size_t>(ceil(4 + 6 * sqrt(lam) + lam));

  //compute the powers of R
  power_.resize(nMax + 1);
  power_[0] = I;
  for (size_t i = 1; i < nMax + 1; ++i)
    MatrixTools::mult(power_[i - 1], R, power_[i]);

  for (size_t i = 0; i < register_->getNumberOfSubstitutionTypes(); ++i) {
    s_[i].resize(nMax + 1);
    MatrixTools::mult(bMatrices_[i], power_[0], s_[i][0]);
    RowMatrix<double> tmp(nbStates_, nbStates_);
    for (size_t l = 1; l < nMax + 1; ++l) {
      MatrixTools::mult(R, s_[i][l - 1], s_[i][l]);
      MatrixTools::mult(bMatrices_[i], power_[l], tmp);
      MatrixTools::add(s_[i][l], tmp);
    }
    MatrixTools::fill(counts_[i], 0);
    for (size_t l = 0; l < nMax + 1; ++l) {
      tmp = s_[i][l];
      //double f = (pow(lam, static_cast<double>(l + 1)) * exp(-lam) / static_cast<double>(NumTools::fact(l + 1))) / miu_;
      double logF = static_cast<double>(l + 1) * log(lam) - lam - log(miu_) - NumTools::logFact(static_cast<double>(l + 1));
      MatrixTools::scale(tmp, exp(logF));
      MatrixTools::add(counts_[i], tmp);
    }
  }

  // Now we must divide by pijt and account for putative weights:
  vector<int> supportedStates = model_->getAlphabetStates();
  RowMatrix<double> P = model_->getPij_t(length);
  for (size_t i = 0; i < register_->getNumberOfSubstitutionTypes(); i++) {
    for (size_t j = 0; j < nbStates_; j++) {
      for(size_t k = 0; k < nbStates_; k++) {
        counts_[i](j, k) /= P(j, k);
        if (isnan(counts_[i](j, k)) || counts_[i](j, k) < 0.)
          counts_[i](j, k) = 0;
        //Weights:
        if (weights_)
          counts_[i](j, k) *= weights_->getIndex(supportedStates[j], supportedStates[k]);
      }
    }
  }
}

/******************************************************************************/

Matrix<double>* UniformizationSubstitutionCount::getAllNumbersOfSubstitutions(double length, size_t type) const
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

double UniformizationSubstitutionCount::getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type) const
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

std::vector<double> UniformizationSubstitutionCount::getNumberOfSubstitutionsForEachType(size_t initialState, size_t finalState, double length) const
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
  size_t n = model->getNumberOfStates();
  if (n != nbStates_) {
    nbStates_ = n;
    //Re-initialize all B matrices according to substitution register.
    initBMatrices_();
  }
  fillBMatrices_();
  
  miu_ = 0;
  for (size_t i = 0; i < nbStates_; ++i) {
    double diagQ = abs(model_->Qij(i, i));
    if (diagQ > miu_)
      miu_ = diagQ;
  }

  if (miu_ > 10000)
    throw Exception("UniformizationSubstitutionCount::setSubstitutionModel(). The maximum diagonal values of generator is above 10000. Abort, chose another mapping method.");

  //Recompute counts:
  computeCounts_(currentLength_);
}

/******************************************************************************/

void UniformizationSubstitutionCount::substitutionRegisterHasChanged() throw (Exception)
{
  //Check compatiblity between model and substitution register:
  if (model_->getAlphabet()->getAlphabetType() != register_->getAlphabet()->getAlphabetType())
    throw Exception("UniformizationSubstitutionCount::substitutionRegisterHasChanged: alphabets do not match between register and model.");

  resetBMatrices_();
  initBMatrices_();
  fillBMatrices_();
  
  //Recompute counts:
  if (currentLength_ > 0)
    computeCounts_(currentLength_);
}

/******************************************************************************/

void UniformizationSubstitutionCount::weightsHaveChanged() throw (Exception)
{
  if (weights_->getAlphabet()->getAlphabetType() != register_->getAlphabet()->getAlphabetType())
    throw Exception("UniformizationSubstitutionCount::weightsHaveChanged. Incorrect alphabet type.");

  //jdutheil on 25/07/14: not necessary if weights are only accounted for in the end.
  //fillBMatrices_();
  
  //Recompute counts:
  if (currentLength_ > 0)
    computeCounts_(currentLength_);
}

/******************************************************************************/

