// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "TwoParameterBinarySubstitutionModel.h"

// From the STL:
#include <cmath>

using namespace bpp;

#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace std;

/******************************************************************************/

TwoParameterBinarySubstitutionModel::TwoParameterBinarySubstitutionModel(
    shared_ptr<const BinaryAlphabet> alpha,
    double mu,
    double pi0) :
  AbstractParameterAliasable("TwoParameterBinary."),
  // AbstractReversibleSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "TwoParameterBinary."),
  AbstractReversibleSubstitutionModel(alpha, std::shared_ptr<const StateMapInterface>(new CanonicalStateMap(alpha, false)), "TwoParameterBinary."),
  mu_(mu),
  pi0_(pi0),
  lambda_(0),
  exp_(0),
  p_(size_, size_)
{
  addParameter_(new Parameter(getNamespace() + "mu", mu_, std::make_shared<IntervalConstraint>(NumConstants::MILLI(), 100, false, false)));
  addParameter_(new Parameter(getNamespace() + "pi0", pi0_, std::make_shared<IntervalConstraint>(0.05, 0.95, true, true)));
  updateMatrices_();
}

/******************************************************************************/

void TwoParameterBinarySubstitutionModel::updateMatrices_()
{
  mu_ = getParameterValue("mu");
  rate_ = mu_;
  pi0_ = getParameterValue("pi0");
  lambda_ = 1;

  // Frequences:
  freq_[0] = pi0_;
  freq_[1] = 1 - pi0_;

  // Generator:
  generator_(0, 0) = -1 * rate_ * freq_[1];
  generator_(0, 1) = rate_ * freq_[1];
  generator_(1, 0) = rate_ * freq_[0];
  generator_(1, 1) = -1 * rate_ * freq_[0];

  // Eigen values:
  eigenValues_[0] = 0;
  eigenValues_[1] = -1 * mu_;

  // Eigen vectors:
  leftEigenVectors_(0, 0) = pi0_;
  leftEigenVectors_(0, 1) = 1 - pi0_;
  leftEigenVectors_(1, 0) = 1;
  leftEigenVectors_(1, 1) = -1;

  rightEigenVectors_(0, 0) = 1;
  rightEigenVectors_(1, 0) = 1;
  rightEigenVectors_(0, 1) = 1 - pi0_;
  rightEigenVectors_(1, 1) = -1 * pi0_;
}

/******************************************************************************/

double TwoParameterBinarySubstitutionModel::Pij_t(size_t i, size_t j, double d) const
{
  exp_ = exp(-lambda_ * rate_ * d);

  switch (i)
  {
  case 0:
    switch (j)
    {
    case 0: return (1 - pi0_) + pi0_ * exp_;
    case 1: return pi0_ * (1 - exp_);
    default: return 0;
    }
  case 1:
    switch (j)
    {
    case 0: return (1 - pi0_) * (1 - exp_);
    case 1: return pi0_ + (1 - pi0_) * exp_;
    default: return 0;
    }
  default: return 0;
  }
}

/******************************************************************************/

double TwoParameterBinarySubstitutionModel::dPij_dt(size_t i, size_t j, double d) const
{
  exp_ = rate_ * exp(-lambda_ * rate_ * d);

  switch (i)
  {
  case 0:
    switch (j)
    {
    case 0: return -1 * pi0_ * exp_;
    case 1: return pi0_ * exp_;
    default: return 0;
    }
  case 1:
    switch (j)
    {
    case 0: return (1 - pi0_) * exp_;
    case 1: return -1 * (1 - pi0_) * exp_;
    default: return 0;
    }
  default: return 0;
  }
}

/******************************************************************************/

double TwoParameterBinarySubstitutionModel::d2Pij_dt2(size_t i, size_t j, double d) const
{
  exp_ = rate_ * rate_ * exp(-lambda_ * rate_ * d);

  switch (i)
  {
  case 0:
    switch (j)
    {
    case 0: return pi0_ * exp_;
    case 1: return -1 * pi0_ * exp_;
    default: return 0;
    }
  case 1:
    switch (j)
    {
    case 0: return -1 * (1 - pi0_) * exp_;
    case 1: return (1 - pi0_) * exp_;
    default: return 0;
    }
  default: return 0;
  }
  return 0;
}

/******************************************************************************/

const Matrix<double>& TwoParameterBinarySubstitutionModel::getPij_t(double d) const
{
  exp_ = exp(-lambda_ * rate_ * d);

  p_(0, 0) = (1 - pi0_) + pi0_ * exp_;
  p_(0, 1) = pi0_ * (1 - exp_);

  p_(1, 0) =  (1 - pi0_) * (1 - exp_);
  p_(1, 1) = pi0_ + (1 - pi0_) * exp_;

  return p_;
}

/******************************************************************************/

const Matrix<double>& TwoParameterBinarySubstitutionModel::getdPij_dt(double d) const
{
  exp_ = rate_ * exp(-lambda_ * rate_ * d);

  p_(0, 0) = -1 * pi0_ * exp_;
  p_(0, 1) = pi0_ * exp_;

  p_(1, 0) = (1 - pi0_) * exp_;
  p_(1, 1) = -1 * (1 - pi0_) * exp_;

  return p_;
}

/******************************************************************************/

const Matrix<double>& TwoParameterBinarySubstitutionModel::getd2Pij_dt2(double d) const
{
  exp_ = rate_ * rate_ * exp(-lambda_ * rate_ * d);

  p_(0, 0) = pi0_ * exp_;
  p_(0, 1) = -1 * pi0_ * exp_;
  p_(1, 0) = -1 * (1 - pi0_) * exp_;
  p_(1, 1) = (1 - pi0_) * exp_;

  return p_;
}

/******************************************************************************/

void TwoParameterBinarySubstitutionModel::setMuBounds(double lb, double ub)
{
  std::shared_ptr<IntervalConstraint> bounds(new IntervalConstraint(lb, ub, true, true));
  getParameter_("mu").setConstraint(bounds);
}
