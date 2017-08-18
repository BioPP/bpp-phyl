//
// File: DecompositionSubstitutionCount.h
// Created by: Julien Dutheil
// Created on: Thu Mar 17 16:08 2011
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

#include "DecompositionSubstitutionCount.h"

#include "Bpp/Numeric/Matrix/MatrixTools.h"
#include <vector>
#include <typeinfo>

using namespace bpp;
using namespace std;

/******************************************************************************/

DecompositionSubstitutionCount::DecompositionSubstitutionCount(const SubstitutionModel* model, SubstitutionRegister* reg, const AlphabetIndex2* weights) :
  AbstractSubstitutionCount(reg),
  AbstractWeightedSubstitutionCount(weights, true),
  DecompositionMethods(model, reg),
  counts_(reg->getNumberOfSubstitutionTypes()),
  currentLength_(0)
{
  //Check compatiblity between model and substitution register:
  if (typeid(model->getAlphabet())!=typeid(reg->getAlphabet()))
    throw Exception("DecompositionSubstitutionCount (constructor): alphabets do not match between register and model.");

  initBMatrices_();
  initCounts_();

  fillBMatrices_();
  computeProducts_();
}

void DecompositionSubstitutionCount::initCounts_()
{
  counts_.resize(register_->getNumberOfSubstitutionTypes());
  //Re-initialize all count matrices according to substitution register.
  for (size_t i = 0; i < register_->getNumberOfSubstitutionTypes(); ++i) {
    counts_[i].resize(nbStates_, nbStates_);
  }
}

/*************************************************/

void DecompositionSubstitutionCount::fillBMatrices_()
{
  vector<int> supportedStates = model_->getAlphabetStates();
  for (size_t j = 0; j < nbStates_; ++j) {
    for (size_t k = 0; k < nbStates_; ++k) {
      size_t i = register_->getType(j, k);
      if (i > 0 && k != j) {
        bMatrices_[i - 1](j, k) = model_->Qij(j, k);
      }
    }
  }
}


/******************************************************************************/

void DecompositionSubstitutionCount::computeCounts_(double length) const
{
  computeExpectations(counts_, length);

  // Now we must divide by pijt and account for putative weights:
  vector<int> supportedStates = model_->getAlphabetStates();
  RowMatrix<double> P = model_->getPij_t(length);
  for (size_t i = 0; i < nbTypes_; i++) {
    for (size_t j = 0; j < nbStates_; j++) {
      for (size_t k = 0; k < nbStates_; k++) {
        counts_[i](j, k) /= P(j, k);
        if (std::isnan(counts_[i](j, k)) || counts_[i](j, k) < 0.) {
          counts_[i](j, k) = 0.;
          //Weights:
          if (weights_)
            counts_[i](j, k) *= weights_->getIndex(supportedStates[j], supportedStates[k]);
        }
      }
    }
  }
}

/******************************************************************************/

Matrix<double>* DecompositionSubstitutionCount::getAllNumbersOfSubstitutions(double length, size_t type) const
{
  if (length < 0)
    throw Exception("DecompositionSubstitutionCount::getAllNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeCounts_(length);
    currentLength_ = length;
  }
  return new RowMatrix<double>(counts_[type - 1]);
}

/******************************************************************************/

double DecompositionSubstitutionCount::getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type) const
{
  if (length < 0)
    throw Exception("DecompositionSubstitutionCount::getNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeCounts_(length);
    currentLength_ = length;
  }
  return counts_[type - 1](initialState, finalState);
}

/******************************************************************************/

std::vector<double> DecompositionSubstitutionCount::getNumberOfSubstitutionsForEachType(size_t initialState, size_t finalState, double length) const
{
  if (length < 0)
    throw Exception("DecompositionSubstitutionCount::getNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeCounts_(length);
    currentLength_ = length;
  }
  std::vector<double> v(getNumberOfSubstitutionTypes());
  for (size_t t = 0; t < getNumberOfSubstitutionTypes(); ++t) {
    v[t] = counts_[t](initialState, finalState);
  }
  return v;
}
    
/******************************************************************************/

void DecompositionSubstitutionCount::setSubstitutionModel(const SubstitutionModel* model)
{
  //Check compatiblity between model and substitution register:
  if (typeid(model->getAlphabet()) != typeid(register_->getAlphabet()))
    throw Exception("DecompositionMethods::setSubstitutionModel: alphabets do not match between register and model.");

  DecompositionMethods::setSubstitutionModel(model);

  initCounts_();
  
  fillBMatrices_();
  computeProducts_();
  
  //Recompute counts:
  computeCounts_(currentLength_);
}

/******************************************************************************/

void DecompositionSubstitutionCount::substitutionRegisterHasChanged() throw (Exception)
{
//Check compatiblity between model and substitution register:
  if (model_->getAlphabet()->getAlphabetType() != register_->getAlphabet()->getAlphabetType())
    throw Exception("DecompositionMethods::substitutionRegisterHasChanged: alphabets do not match between register and model.");

  initBMatrices_();
  initStates_();
  
  initCounts_();

  fillBMatrices_();
  computeProducts_();
  
  //Recompute counts:
  if (currentLength_ > 0)
    computeCounts_(currentLength_);
}

/******************************************************************************/

void DecompositionSubstitutionCount::weightsHaveChanged() throw (Exception)
{
  if (typeid(weights_->getAlphabet()) != typeid(register_->getAlphabet()))
    throw Exception("DecompositionSubstitutionCount::weightsHaveChanged. Incorrect alphabet type.");

  //Recompute counts:
  if (currentLength_ > 0)
    computeCounts_(currentLength_);
}

/******************************************************************************/

