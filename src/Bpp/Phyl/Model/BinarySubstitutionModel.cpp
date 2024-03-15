// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "BinarySubstitutionModel.h"

// From the STL:
#include <cmath>

using namespace bpp;

#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace std;

/******************************************************************************/

BinarySubstitutionModel::BinarySubstitutionModel(
    shared_ptr<const BinaryAlphabet> alpha,
    double kappa) :
  AbstractParameterAliasable("Binary."),
  // AbstractSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "Binary."),
  AbstractReversibleSubstitutionModel(alpha, std::shared_ptr<const StateMapInterface>(new CanonicalStateMap(alpha, false)), "Binary."),
  kappa_(kappa),
  lambda_(0),
  exp_(0),
  p_(size_, size_)
{
  addParameter_(new Parameter(getNamespace() + "kappa", kappa_, Parameter::R_PLUS_STAR));
  updateMatrices_();
}

/******************************************************************************/

void BinarySubstitutionModel::updateMatrices_()
{
  kappa_ = getParameterValue("kappa"); // alpha/beta
  lambda_ = (kappa_ + 1) * (kappa_ + 1) / (2 * kappa_);

  // Frequences:
  freq_[0] = 1 / (kappa_ + 1);
  freq_[1] = kappa_ / (kappa_ + 1);

  // Generator:
  generator_(0, 0) = rate_ * -(kappa_ + 1) / 2;
  generator_(0, 1) = rate_ * (kappa_ + 1) / 2;
  generator_(1, 0) = rate_ * (kappa_ + 1) / (2 * kappa_);
  generator_(1, 1) = rate_ * -(kappa_ + 1) / (2 * kappa_);

  // Eigen values:
  eigenValues_[0] = 0;
  eigenValues_[1] = -rate_ * lambda_;

  // Eigen vectors:
  leftEigenVectors_(0, 0) = 1 / (kappa_ + 1);
  leftEigenVectors_(0, 1) = kappa_ / (kappa_ + 1);
  if (kappa_ != 1.0)
  {
    leftEigenVectors_(1, 0) = (kappa_ - 1) / (kappa_ + 1);
    leftEigenVectors_(1, 1) = -(kappa_ - 1) / (kappa_ + 1);
  }
  else
  {
    leftEigenVectors_(1, 0) = 1;
    leftEigenVectors_(1, 1) = -1;
  }

  rightEigenVectors_(0, 0) = 1;
  rightEigenVectors_(1, 0) = 1;

  if (kappa_ != 1.0)
  {
    rightEigenVectors_(0, 1) = kappa_ / (kappa_ - 1);
    rightEigenVectors_(1, 1) = -1 / (kappa_ - 1);
  }
  else
  {
    rightEigenVectors_(0, 1) = 1 / 2;
    rightEigenVectors_(1, 1) = -1 / 2;
  }
}

/******************************************************************************/

double BinarySubstitutionModel::Pij_t(size_t i, size_t j, double d) const
{
  exp_ = exp(-lambda_ * rate_ * d);

  switch (i)
  {
  case 0:
    switch (j)
    {
    case 0: return (1 + kappa_ * exp_) / (kappa_ + 1);
    case 1: return kappa_ / (kappa_ + 1) * (1 - exp_);
    default: return 0;
    }
  case 1:
    switch (j)
    {
    case 0: return (1 - exp_) / (kappa_ + 1);
    case 1: return (kappa_ + exp_) / (kappa_ + 1);
    default: return 0;
    }
  default: return 0;
  }
}

/******************************************************************************/

double BinarySubstitutionModel::dPij_dt(size_t i, size_t j, double d) const
{
  exp_ = rate_ * exp(-lambda_ * rate_ * d);

  switch (i)
  {
  case 0:
    switch (j)
    {
    case 0: return -(kappa_ + 1) / 2 * exp_;
    case 1: return (kappa_ + 1) / 2 * exp_;
    default: return 0;
    }
  case 1:
    switch (j)
    {
    case 0: return (kappa_ + 1) / (2 * kappa_) * exp_;
    case 1: return -(kappa_ + 1) / (2 * kappa_) * exp_;
    default: return 0;
    }
  default: return 0;
  }
}

/******************************************************************************/

double BinarySubstitutionModel::d2Pij_dt2(size_t i, size_t j, double d) const
{
  exp_ = rate_ * rate_ * exp(-lambda_ * rate_ * d);

  switch (i)
  {
  case 0:
    switch (j)
    {
    case 0: return lambda_ * (kappa_ + 1) / 2 * exp_;
    case 1: return -lambda_ * (kappa_ + 1) / 2 * exp_;
    default: return 0;
    }
  case 1:
    switch (j)
    {
    case 0: return -lambda_ * (kappa_ + 1) / (2 * kappa_) * exp_;
    case 1: return lambda_ * (kappa_ + 1) / (2 * kappa_) * exp_;
    default: return 0;
    }
  default: return 0;
  }
}

/******************************************************************************/

const Matrix<double>& BinarySubstitutionModel::getPij_t(double d) const
{
  exp_ = exp(-lambda_ * rate_ * d);

  p_(0, 0) = (1 + kappa_ * exp_) / (kappa_ + 1);
  p_(0, 1) = kappa_ / (kappa_ + 1) * (1 - exp_);

  p_(1, 0) =  (1 - exp_) / (kappa_ + 1);
  p_(1, 1) = (kappa_ + exp_) / (kappa_ + 1);

  return p_;
}

const Matrix<double>& BinarySubstitutionModel::getdPij_dt(double d) const
{
  exp_ = rate_ * exp(-lambda_ * rate_ * d);

  p_(0, 0) = -(kappa_ + 1) / 2 * exp_;
  p_(0, 1) = (kappa_ + 1) / 2 * exp_;

  p_(1, 0) = (kappa_ + 1) / (2 * kappa_) * exp_;
  p_(1, 1) = -(kappa_ + 1) / (2 * kappa_) * exp_;

  return p_;
}

const Matrix<double>& BinarySubstitutionModel::getd2Pij_dt2(double d) const
{
  exp_ = rate_ * rate_ * exp(-lambda_ * rate_ * d);

  p_(0, 0) = lambda_ * (kappa_ + 1) / 2 * exp_;
  p_(0, 1) = -lambda_ * (kappa_ + 1) / 2 * exp_;
  p_(1, 0) = -lambda_ * (kappa_ + 1) / (2 * kappa_) * exp_;
  p_(1, 1) = lambda_ * (kappa_ + 1) / (2 * kappa_) * exp_;

  return p_;
}

/******************************************************************************/

void BinarySubstitutionModel::setFreq(std::map<int, double>& freqs)
{
  kappa_ = freqs[1] / freqs[0];
  setParameterValue("kappa", kappa_);
  updateMatrices_();
}
