//
// File: MarkovModulatedSubstitutionModel.cpp
// Created by: Julien Dutheil
// Created on: Sat Aug 05 08:21 2006
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "MarkovModulatedSubstitutionModel.h"

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Matrix/EigenValue.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

MarkovModulatedSubstitutionModel::MarkovModulatedSubstitutionModel(
  const MarkovModulatedSubstitutionModel& model) :
  AbstractParameterAliasable(model),
  model_               (dynamic_cast<ReversibleSubstitutionModel*>(model.model_->clone())),
  stateMap_            (model.stateMap_),
  nbStates_            (model.nbStates_),
  nbRates_             (model.nbRates_),
  rates_               (model.rates_),
  ratesExchangeability_(model.ratesExchangeability_),
  ratesFreq_           (model.ratesFreq_),
  ratesGenerator_      (model.ratesGenerator_),
  generator_           (model.generator_),
  exchangeability_     (model.exchangeability_),
  leftEigenVectors_    (model.leftEigenVectors_),
  rightEigenVectors_   (model.rightEigenVectors_),
  eigenValues_         (model.eigenValues_),
  iEigenValues_        (model.iEigenValues_),
  eigenDecompose_      (model.eigenDecompose_),
  pijt_                (model.pijt_),
  dpijt_               (model.dpijt_),
  d2pijt_              (model.d2pijt_),
  freq_                (model.freq_),
  normalizeRateChanges_(model.normalizeRateChanges_),
  nestedPrefix_        (model.nestedPrefix_)
{}

MarkovModulatedSubstitutionModel& MarkovModulatedSubstitutionModel::operator=(
  const MarkovModulatedSubstitutionModel& model)
{
  AbstractParametrizable::operator=(model);
  model_                = dynamic_cast<ReversibleSubstitutionModel*>(model.model_->clone());
  stateMap_             = model.stateMap_;
  nbStates_             = model.nbStates_;
  nbRates_              = model.nbRates_;
  rates_                = model.rates_;
  ratesExchangeability_ = model.ratesExchangeability_;
  ratesFreq_            = model.ratesFreq_;
  ratesGenerator_       = model.ratesGenerator_;
  generator_            = model.generator_;
  exchangeability_      = model.exchangeability_;
  leftEigenVectors_     = model.leftEigenVectors_;
  rightEigenVectors_    = model.rightEigenVectors_;
  eigenValues_          = model.eigenValues_;
  iEigenValues_         = model.iEigenValues_;
  eigenDecompose_       = model.eigenDecompose_;
  pijt_                 = model.pijt_;
  dpijt_                = model.dpijt_;
  d2pijt_               = model.d2pijt_;
  freq_                 = model.freq_;
  normalizeRateChanges_ = model.normalizeRateChanges_;
  nestedPrefix_         = model.nestedPrefix_;
  return *this;
}

/******************************************************************************/

void MarkovModulatedSubstitutionModel::updateMatrices()
{
  // ratesGenerator_ and rates_ must be initialized!
  nbStates_        = model_->getNumberOfStates();
  nbRates_         = rates_.getNumberOfColumns();
  RowMatrix<double> Tmp1, Tmp2;
  MatrixTools::diag(ratesFreq_, Tmp1);
  MatrixTools::mult(ratesExchangeability_, Tmp1, ratesGenerator_);
  MatrixTools::kroneckerMult(rates_, model_->getGenerator(), generator_);

  MatrixTools::MatrixTools::getId< RowMatrix<double> >(nbStates_, Tmp1);
  MatrixTools::kroneckerMult(ratesGenerator_, Tmp1, Tmp2);
  MatrixTools::add(generator_, Tmp2);

  MatrixTools::diag(1. / ratesFreq_, Tmp1);
  MatrixTools::mult(rates_, Tmp1, Tmp2);
  MatrixTools::kroneckerMult(Tmp2, model_->getExchangeabilityMatrix(), exchangeability_);

  MatrixTools::diag(1 / model_->getFrequencies(), Tmp1);
  MatrixTools::kroneckerMult(ratesExchangeability_, Tmp1, Tmp2);
  MatrixTools::add(exchangeability_, Tmp2);
  freq_ = VectorTools::kroneckerMult(ratesFreq_, model_->getFrequencies());
  if (normalizeRateChanges_)
  {
    // Normalization:
    Vdouble Tmp;
    MatrixTools::diag(generator_, Tmp);
    double scale = -VectorTools::scalar<double, double>(Tmp, freq_);
    MatrixTools::scale(generator_, 1. / scale);

    // Normalize exchangeability matrix too:
    MatrixTools::scale(exchangeability_, 1. / scale);
  }

  // Compute eigen values and vectors:
  eigenValues_.resize(nbRates_ * nbStates_);
  iEigenValues_.resize(nbRates_ * nbStates_);
  rightEigenVectors_.resize(nbStates_ * nbRates_, nbStates_ * nbRates_);
  pijt_.resize(nbStates_ * nbRates_, nbStates_ * nbRates_);
  dpijt_.resize(nbStates_ * nbRates_, nbStates_ * nbRates_);
  d2pijt_.resize(nbStates_ * nbRates_, nbStates_ * nbRates_);

  vector<double>    modelEigenValues       = model_->getEigenValues();
  RowMatrix<double> modelRightEigenVectors = model_->getColumnRightEigenVectors();
  for (unsigned int i = 0; i < nbStates_; i++)
  {
    RowMatrix<double> tmp = rates_;
    MatrixTools::scale(tmp, modelEigenValues[i]);
    MatrixTools::add(tmp, ratesGenerator_);
    EigenValue<double> ev(tmp);
    vector<double>    values  = ev.getRealEigenValues();
    RowMatrix<double> vectors = ev.getV();
    for (size_t j = 0; j < nbRates_; j++)
    {
      size_t c = i * nbRates_ + j; // Current eigen value index.
      eigenValues_[c] = values[j];
      // Compute the Kronecker product of the jth vector and the ith modelRightEigenVector.
      for (unsigned int ii = 0; ii < nbRates_; ii++)
      {
        double vii = vectors(ii, j);
        for (unsigned int jj = 0; jj < nbStates_; jj++)
        {
          rightEigenVectors_(ii * nbStates_ + jj, c) = vii * modelRightEigenVectors(jj, i);
        }
      }
    }
  }
  // Now compute left eigen vectors by inversion:
  MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);
}

/******************************************************************************/

const Matrix<double>& MarkovModulatedSubstitutionModel::getPij_t(double t) const
{
  if (t == 0)
    MatrixTools::getId< RowMatrix<double> >(nbStates_ * nbRates_, pijt_);
  else
    MatrixTools::mult(rightEigenVectors_, VectorTools::exp(eigenValues_ * t), leftEigenVectors_, pijt_);
  return pijt_;
}

const Matrix<double>& MarkovModulatedSubstitutionModel::getdPij_dt(double t) const
{
  MatrixTools::mult(rightEigenVectors_, eigenValues_ * VectorTools::exp(eigenValues_ * t), leftEigenVectors_, dpijt_);
  return dpijt_;
}

const Matrix<double>& MarkovModulatedSubstitutionModel::getd2Pij_dt2(double t) const
{
  MatrixTools::mult(rightEigenVectors_, VectorTools::sqr(eigenValues_) * VectorTools::exp(eigenValues_ * t), leftEigenVectors_, d2pijt_);
  return d2pijt_;
}

/******************************************************************************/

double MarkovModulatedSubstitutionModel::getInitValue(size_t i, int state) const throw (IndexOutOfBoundsException, BadIntException)
{
  if (i >= (nbStates_ * nbRates_))
    throw IndexOutOfBoundsException("MarkovModulatedSubstitutionModel::getInitValue", i, 0, nbStates_ * nbRates_ - 1);
  if (state < 0 || !model_->getAlphabet()->isIntInAlphabet(state))
    throw BadIntException(state, "MarkovModulatedSubstitutionModel::getInitValue. Character " + model_->getAlphabet()->intToChar(state) + " is not allowed in model.");
  vector<int> states = model_->getAlphabet()->getAlias(state);
  for (size_t j = 0; j < states.size(); j++)
  {
    if (getAlphabetStateAsInt(i) == states[j])
      return 1.;
  }
  return 0.;
}

/******************************************************************************/

void MarkovModulatedSubstitutionModel::setNamespace(const string& prefix)
{
  AbstractParameterAliasable::setNamespace(prefix);
  // We also need to update the namespace of the nested model:
  model_->setNamespace(prefix + nestedPrefix_);
}

/******************************************************************************/

