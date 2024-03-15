// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/EigenValue.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>

#include "MarkovModulatedSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

MarkovModulatedSubstitutionModel::MarkovModulatedSubstitutionModel(
  const MarkovModulatedSubstitutionModel& model) :
  AbstractParameterAliasable(model),
  model_               (model.model_->clone()),
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
  compFreq_            (model.compFreq_),
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
  model_.reset(model.model_->clone());
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
  compFreq_             = model.compFreq_;
  pijt_                 = model.pijt_;
  dpijt_                = model.dpijt_;
  d2pijt_               = model.d2pijt_;
  freq_                 = model.freq_;
  normalizeRateChanges_ = model.normalizeRateChanges_;
  nestedPrefix_         = model.nestedPrefix_;
  return *this;
}

/******************************************************************************/

void MarkovModulatedSubstitutionModel::updateMatrices_()
{
  // ratesGenerator_ and rates_ must be initialized!
  nbStates_        = model_->getNumberOfStates();
  nbRates_         = rates_.getNumberOfColumns();
  RowMatrix<double> Tmp1, Tmp2;
  MatrixTools::diag(ratesFreq_, Tmp1);
  MatrixTools::mult(ratesExchangeability_, Tmp1, ratesGenerator_);
  MatrixTools::kroneckerMult(rates_, model_->generator(), generator_);

  MatrixTools::MatrixTools::getId< RowMatrix<double> >(nbStates_, Tmp1);
  MatrixTools::kroneckerMult(ratesGenerator_, Tmp1, Tmp2);
  MatrixTools::add(generator_, Tmp2);

  MatrixTools::diag(1. / ratesFreq_, Tmp1);
  MatrixTools::mult(rates_, Tmp1, Tmp2);
  MatrixTools::kroneckerMult(Tmp2, model_->exchangeabilityMatrix(), exchangeability_);

  MatrixTools::diag(1 / model_->getFrequencies(), Tmp1);
  MatrixTools::kroneckerMult(ratesExchangeability_, Tmp1, Tmp2);
  MatrixTools::add(exchangeability_, Tmp2);
  freq_ = VectorTools::kroneckerMult(ratesFreq_, model_->getFrequencies());
  if (normalizeRateChanges_)
  {
    // Normalization:
    double scale = getScale();
    setScale(1. / scale);

    // Normalize exchangeability matrix too:
    if (isScalable())
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

void MarkovModulatedSubstitutionModel::setDiagonal()
{
  for (size_t i = 0; i < getNumberOfStates(); i++)
  {
    double lambda = 0;
    Vdouble& row = generator_.getRow(i);

    for (size_t j = 0; j < getNumberOfStates(); j++)
    {
      if (j != i)
        lambda += row[j];
    }
    row[i] = -lambda;
  }
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

double MarkovModulatedSubstitutionModel::getInitValue(size_t i, int state) const
{
  if (i >= (nbStates_ * nbRates_))
    throw IndexOutOfBoundsException("MarkovModulatedSubstitutionModel::getInitValue", i, 0, nbStates_ * nbRates_ - 1);
  if (state < 0 || !model_->getAlphabet()->isIntInAlphabet(state))
    throw BadIntException(state, "MarkovModulatedSubstitutionModel::getInitValue. Character " + model_->getAlphabet()->intToChar(state) + " is not allowed in model.", getAlphabet().get());
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
