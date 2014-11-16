//
// File: DecompositionReward.h
// Created by: Laurent Guéguen
// Created on: mercredi 27 mars 2013, à 12h 36
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

#include "DecompositionReward.h"

#include "Bpp/Numeric/Matrix/MatrixTools.h"
#include <vector>

using namespace bpp;
using namespace std;

/******************************************************************************/

DecompositionReward::DecompositionReward(const SubstitutionModel* model, AlphabetIndex1* alphIndex) :
  AbstractReward(alphIndex),
  model_(dynamic_cast<const ReversibleSubstitutionModel*>(model)),
  nbStates_(model->getNumberOfStates()),
  jMat_(nbStates_, nbStates_),
  v_(nbStates_, nbStates_),
  vInv_(nbStates_, nbStates_),
  lambda_(nbStates_, nbStates_),
  bMatrice_(nbStates_, nbStates_),
  insideProduct_(nbStates_, nbStates_),
  rewards_(nbStates_, nbStates_),
  currentLength_(-1.)
{
  //Check compatiblity between model and alphabet Index:
  if (model->getAlphabet()->getAlphabetType() != alphIndex_->getAlphabet()->getAlphabetType())
    throw Exception("DecompositionReward (constructor): alphabets do not match between alphabet index and model.");
  if (!dynamic_cast<const ReversibleSubstitutionModel*>(model))
    throw Exception("DecompositionReward::DecompositionReward. Only works with declared reversible models for now.");

  //Initialize the B matrice. This is done once for all,
  //unless the number of states changes:
  computeBMatrice_();
  computeEigen_();
  computeProducts_();
}				
    
/******************************************************************************/

void DecompositionReward::computeBMatrice_()
{
  vector<int> supportedStates = model_->getAlphabetStates();
  for (size_t j = 0; j < nbStates_; ++j) 
    bMatrice_(j, j) = getAlphabetIndex()->getIndex(supportedStates[j]);
}

void DecompositionReward::computeEigen_()
{
  v_      = model_->getColumnRightEigenVectors();
  vInv_   = model_->getRowLeftEigenVectors();
  lambda_ = model_->getEigenValues();
}

void DecompositionReward::computeProducts_()
{
  RowMatrix<double> tmp(nbStates_, nbStates_);
  MatrixTools::mult(vInv_, bMatrice_, tmp);
  MatrixTools::mult(tmp, v_, insideProduct_);
}

void DecompositionReward::resetStates_()
{
  jMat_.resize(nbStates_, nbStates_);
  v_.resize(nbStates_, nbStates_);
  vInv_.resize(nbStates_, nbStates_);
  lambda_.resize(nbStates_);
  bMatrice_.resize(nbStates_, nbStates_);
  insideProduct_.resize(nbStates_, nbStates_);
  rewards_.resize(nbStates_, nbStates_);
}

void DecompositionReward::jFunction_(const std::vector<double>& lambda, double t, RowMatrix<double>& result) const
{
  vector<double> expLam = VectorTools::exp(lambda * t);
  for (unsigned int i = 0; i < nbStates_; ++i) {
    for (unsigned int j = 0; j < nbStates_; ++j) {
      double dd = lambda[i] - lambda[j];
      if (dd == 0) {
        result(i, j) = t * expLam[i];
      } else {
        result(i, j) = (expLam[i] - expLam[j]) / dd;
      }
    }
  }
}

/******************************************************************************/

void DecompositionReward::computeRewards_(double length) const
{
  jFunction_(lambda_, length, jMat_);
  RowMatrix<double> tmp1(nbStates_, nbStates_), tmp2(nbStates_, nbStates_);
  MatrixTools::hadamardMult(jMat_, insideProduct_, tmp1);
  MatrixTools::mult(v_, tmp1, tmp2);
  MatrixTools::mult(tmp2, vInv_, rewards_);

  // Now we must divide by pijt:
  RowMatrix<double> P = model_->getPij_t(length);
  for (size_t j = 0; j < nbStates_; j++) {
    for (size_t k = 0; k < nbStates_; k++) {
      rewards_(j, k) /= P(j, k);
      if (isnan(rewards_(j, k)))
        rewards_(j, k) = 0.;
    }
  }
}

/******************************************************************************/

Matrix<double>* DecompositionReward::getAllRewards(double length) const
{
  if (length < 0)
    throw Exception("DecompositionReward::getAllRewards. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
    {
      computeRewards_(length);
      currentLength_ = length;
    }
  return new RowMatrix<double>(rewards_);
}

/******************************************************************************/

double DecompositionReward::getReward(size_t initialState, size_t finalState, double length) const
{
  if (length < 0)
    throw Exception("DecompositionReward::getRewards. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeRewards_(length);
    currentLength_ = length;
  }
  return rewards_(initialState, finalState);
}

/******************************************************************************/

void DecompositionReward::setSubstitutionModel(const SubstitutionModel* model)
{
  const ReversibleSubstitutionModel* rModel = dynamic_cast<const ReversibleSubstitutionModel*>(model);
  if (!rModel)
    throw Exception("DecompositionReward::setSubstitutionModel. Only works with reversible models for now.");

  //Check compatiblity between model and substitution register:
  if (model->getAlphabet()->getAlphabetType() != alphIndex_->getAlphabet()->getAlphabetType())
    throw Exception("DecompositionReward::setSubstitutionModel: alphabets do not match between alphabet index and model.");
  model_ = rModel;
  size_t n = model->getNumberOfStates();
  if (n != nbStates_) {
    nbStates_ = n;
    resetStates_();
  }
  computeEigen_();
  computeProducts_();

  //Recompute rewards:
  computeRewards_(currentLength_);
}

/******************************************************************************/

void DecompositionReward::alphabetIndexHasChanged() throw (Exception)
{
  //Check compatiblity between model and substitution register:
  if (model_->getAlphabet()->getAlphabetType() != alphIndex_->getAlphabet()->getAlphabetType())
    throw Exception("DecompositionReward::AlphabetIndexHasChanged: alphabets do not match between alphbaet index and model.");

  computeBMatrice_();
  computeProducts_();

  //Recompute rewards:
  if (currentLength_ > 0)
   computeRewards_(currentLength_);
}


