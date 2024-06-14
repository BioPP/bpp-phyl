// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "JCnuc.h"

using namespace bpp;

#include <cmath>

using namespace std;

/******************************************************************************/

JCnuc::JCnuc(shared_ptr<const NucleicAlphabet> alpha) :
  AbstractParameterAliasable("JC69."),
  AbstractReversibleNucleotideSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "JC69."),
  exp_(),
  p_(size_, size_)
{
  updateMatrices_();
}

/******************************************************************************/

void JCnuc::updateMatrices_()
{
  // Frequencies:
  freq_[0] = freq_[1] = freq_[2] = freq_[3] = 1. / 4.;

  // Generator and exchangeabilities:
  for (size_t i = 0; i < 4; ++i)
  {
    for (size_t j = 0; j < 4; ++j)
    {
      generator_(i, j) = (i == j) ? -1. : 1. / 3.;
      exchangeability_(i, j) = generator_(i, j) * 4.;
    }
  }

  // Eigen values:
  eigenValues_[0] = 0;
  eigenValues_[1] = eigenValues_[2] = eigenValues_[3] = -4. / 3.;

  // Eigen vectors:
  for (size_t i = 0; i < 4; i++)
  {
    leftEigenVectors_(0, i) = 1. / 4.;
  }
  for (size_t i = 1; i < 4; i++)
  {
    for (size_t j = 0; j < 4; j++)
    {
      leftEigenVectors_(i, j) = -1. / 4.;
    }
  }
  leftEigenVectors_(1, 2) = 3. / 4.;
  leftEigenVectors_(2, 1) = 3. / 4.;
  leftEigenVectors_(3, 0) = 3. / 4.;

  for (size_t i = 0; i < 4; i++)
  {
    rightEigenVectors_(i, 0) = 1.;
  }
  for (size_t i = 1; i < 4; i++)
  {
    rightEigenVectors_(3, i) = -1.;
  }
  for (size_t i = 0; i < 3; i++)
  {
    for (size_t j = 1; j < 4; j++)
    {
      rightEigenVectors_(i, j) = 0.;
    }
  }
  rightEigenVectors_(2, 1) = 1.;
  rightEigenVectors_(1, 2) = 1.;
  rightEigenVectors_(0, 3) = 1.;
}

/******************************************************************************/

double JCnuc::Pij_t(size_t i, size_t j, double d) const
{
  if (i == j)
    return 1. / 4. + 3. / 4. * exp(-rate_ * 4. / 3. * d);
  else
    return 1. / 4. - 1. / 4. * exp(-rate_ * 4. / 3. * d);
}

/******************************************************************************/

double JCnuc::dPij_dt(size_t i, size_t j, double d) const
{
  if (i == j)
    return -exp(-rate_ * 4. / 3. * d) * rate_;
  else
    return 1. / 3. * exp(-rate_ * 4. / 3. * d) * rate_;
}

/******************************************************************************/

double JCnuc::d2Pij_dt2(size_t i, size_t j, double d) const
{
  if (i == j)
    return 4. / 3. * exp(-rate_ * 4. / 3. * d) * rate_ * rate_;
  else
    return -4. / 9. * exp(-rate_ * 4. / 3. * d) * rate_ * rate_;
}

/******************************************************************************/

const Matrix<double>& JCnuc::getPij_t(double d) const
{
  exp_ = exp(-4. / 3. * d * rate_);
  for (size_t i = 0; i < size_; i++)
  {
    for (size_t j = 0; j < size_; j++)
    {
      p_(i, j) = (i == j) ? 1. / 4. + 3. / 4. * exp_ : 1. / 4. - 1. / 4. * exp_;
    }
  }
  return p_;
}

const Matrix<double>& JCnuc::getdPij_dt(double d) const
{
  exp_ = exp(-4. / 3. * d * rate_);
  for (size_t i = 0; i < size_; i++)
  {
    for (size_t j = 0; j < size_; j++)
    {
      p_(i, j) = rate_ * ((i == j) ? -exp_ : 1. / 3. * exp_);
    }
  }
  return p_;
}

const Matrix<double>& JCnuc::getd2Pij_dt2(double d) const
{
  exp_ = exp(-4. / 3. * d * rate_);
  for (size_t i = 0; i < size_; i++)
  {
    for (size_t j = 0; j < size_; j++)
    {
      p_(i, j) = rate_ * rate_ * ((i == j) ? 4. / 3. * exp_ : -4. / 9. * exp_);
    }
  }
  return p_;
}

/******************************************************************************/
