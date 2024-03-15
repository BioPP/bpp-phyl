// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "../FrequencySet/NucleotideFrequencySet.h"
#include "F81.h"

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/******************************************************************************/

F81::F81(
  shared_ptr<const NucleicAlphabet> alpha,
  double piA,
  double piC,
  double piG,
  double piT) :
  AbstractParameterAliasable("F81."),
  AbstractReversibleNucleotideSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "F81."),
  r_(),
  piA_(piA), piC_(piC), piG_(piG), piT_(piT), piY_(), piR_(),
  theta_(piG + piC), theta1_(piA / (1. - theta_)), theta2_(piG / theta_),
  exp1_(), l_(), p_(size_, size_)
{
  addParameter_(new Parameter("F81.theta", theta_, FrequencySetInterface::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter("F81.theta1", theta1_, FrequencySetInterface::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter("F81.theta2", theta2_, FrequencySetInterface::FREQUENCE_CONSTRAINT_SMALL));
  computeFrequencies(false);
  updateMatrices_();
}

/******************************************************************************/

void F81::updateMatrices_()
{
  theta_  = getParameterValue("theta");
  theta1_ = getParameterValue("theta1");
  theta2_ = getParameterValue("theta2");
  piA_ = theta1_ * (1. - theta_);
  piC_ = (1. - theta2_) * theta_;
  piG_ = theta2_ * theta_;
  piT_ = (1. - theta1_) * (1. - theta_);
  piR_   = piA_ + piG_;
  piY_   = piT_ + piC_;

  freq_[0] = piA_;
  freq_[1] = piC_;
  freq_[2] = piG_;
  freq_[3] = piT_;

  generator_(0, 0) = piA_ - 1;
  generator_(1, 1) = piC_ - 1;
  generator_(2, 2) = piG_ - 1;
  generator_(3, 3) = piT_ - 1;

  generator_(1, 0) = piA_;
  generator_(3, 0) = piA_;
  generator_(0, 1) = piC_;
  generator_(2, 1) = piC_;
  generator_(1, 2) = piG_;
  generator_(3, 2) = piG_;
  generator_(0, 3) = piT_;
  generator_(2, 3) = piT_;

  generator_(2, 0) = piA_;
  generator_(3, 1) = piC_;
  generator_(0, 2) = piG_;
  generator_(1, 3) = piT_;

  // Normalization:
  r_ = 1. / (2. * (piA_ * piC_ + piC_ * piG_ + piA_ * piT_ + piG_ * piT_ + (piC_ * piT_ + piA_ * piG_)));
  MatrixTools::scale(generator_, r_);

  // Exchangeability:
  exchangeability_(0, 0) = generator_(0, 0) / piA_;
  exchangeability_(0, 1) = generator_(0, 1) / piC_;
  exchangeability_(0, 2) = generator_(0, 2) / piG_;
  exchangeability_(0, 3) = generator_(0, 3) / piT_;

  exchangeability_(1, 0) = generator_(1, 0) / piA_;
  exchangeability_(1, 1) = generator_(1, 1) / piC_;
  exchangeability_(1, 2) = generator_(1, 2) / piG_;
  exchangeability_(1, 3) = generator_(1, 3) / piT_;

  exchangeability_(2, 0) = generator_(2, 0) / piA_;
  exchangeability_(2, 1) = generator_(2, 1) / piC_;
  exchangeability_(2, 2) = generator_(2, 2) / piG_;
  exchangeability_(2, 3) = generator_(2, 3) / piT_;

  exchangeability_(3, 0) = generator_(3, 0) / piA_;
  exchangeability_(3, 1) = generator_(3, 1) / piC_;
  exchangeability_(3, 2) = generator_(3, 2) / piG_;
  exchangeability_(3, 3) = generator_(3, 3) / piT_;

  // Eigen values:
  eigenValues_[0] = 0;
  eigenValues_[1] = -r_;
  eigenValues_[2] = -r_;
  eigenValues_[3] = -r_;

  // Eigen vectors:
  leftEigenVectors_(0, 0) = piA_;
  leftEigenVectors_(0, 1) = piC_;
  leftEigenVectors_(0, 2) = piG_;
  leftEigenVectors_(0, 3) = piT_;

  leftEigenVectors_(1, 0) = 0.;
  leftEigenVectors_(1, 1) = piT_ / piY_;
  leftEigenVectors_(1, 2) = 0.;
  leftEigenVectors_(1, 3) = -piT_ / piY_;

  leftEigenVectors_(2, 0) = piG_ / piR_;
  leftEigenVectors_(2, 1) = 0.;
  leftEigenVectors_(2, 2) = -piG_ / piR_;
  leftEigenVectors_(2, 3) = 0.;

  leftEigenVectors_(3, 0) = piA_ * piY_ / piR_;
  leftEigenVectors_(3, 1) = -piC_;
  leftEigenVectors_(3, 2) = piG_ * piY_ / piR_;
  leftEigenVectors_(3, 3) = -piT_;

  rightEigenVectors_(0, 0) = 1.;
  rightEigenVectors_(0, 1) = 0.;
  rightEigenVectors_(0, 2) = 1.;
  rightEigenVectors_(0, 3) = 1.;

  rightEigenVectors_(1, 0) = 1.;
  rightEigenVectors_(1, 1) = 1.;
  rightEigenVectors_(1, 2) = 0.;
  rightEigenVectors_(1, 3) = -piR_ / piY_;

  rightEigenVectors_(2, 0) = 1.;
  rightEigenVectors_(2, 1) = 0.;
  rightEigenVectors_(2, 2) = -piA_ / piG_;
  rightEigenVectors_(2, 3) = 1.;

  rightEigenVectors_(3, 0) = 1.;
  rightEigenVectors_(3, 1) = -piC_ / piT_;
  rightEigenVectors_(3, 2) = 0.;
  rightEigenVectors_(3, 3) = -piR_ / piY_;
}

/******************************************************************************/

double F81::Pij_t(size_t i, size_t j, double d) const
{
  l_     = rate_ * r_ * d;
  exp1_  = exp(-l_);

  switch (i)
  {
  // A
  case 0: {
    switch (j)
    {
    case 0: return piA_ * (1. + (piY_ / piR_) * exp1_) + (piG_ / piR_) * exp1_; // A
    case 1: return piC_ * (1. -               exp1_);                       // C
    case 2: return piG_ * (1. + (piY_ / piR_) * exp1_) - (piG_ / piR_) * exp1_; // G
    case 3: return piT_ * (1. -               exp1_);                       // T, U
    }
  }
  // C
  case 1: {
    switch (j)
    {
    case 0: return piA_ * (1. -               exp1_);                       // A
    case 1: return piC_ * (1. + (piR_ / piY_) * exp1_) + (piT_ / piY_) * exp1_; // C
    case 2: return piG_ * (1. -               exp1_);                       // G
    case 3: return piT_ * (1. + (piR_ / piY_) * exp1_) - (piT_ / piY_) * exp1_; // T, U
    }
  }
  // G
  case 2: {
    switch (j)
    {
    case 0: return piA_ * (1. + (piY_ / piR_) * exp1_) - (piA_ / piR_) * exp1_; // A
    case 1: return piC_ * (1. -               exp1_);                       // C
    case 2: return piG_ * (1. + (piY_ / piR_) * exp1_) + (piA_ / piR_) * exp1_; // G
    case 3: return piT_ * (1. -               exp1_);                       // T, U
    }
  }
  // T, U
  case 3: {
    switch (j)
    {
    case 0: return piA_ * (1. -               exp1_);                       // A
    case 1: return piC_ * (1. + (piR_ / piY_) * exp1_) - (piC_ / piY_) * exp1_; // C
    case 2: return piG_ * (1. -               exp1_);                       // G
    case 3: return piT_ * (1. + (piR_ / piY_) * exp1_) + (piC_ / piY_) * exp1_; // T, U
    }
  }
  }
  return 0;
}

/******************************************************************************/

double F81::dPij_dt(size_t i, size_t j, double d) const
{
  l_     = rate_ * r_ * d;
  exp1_  = exp(-l_);

  switch (i)
  {
  // A
  case 0: {
    switch (j)
    {
    case 0: return rate_ * r_ * (piA_ * -(piY_ / piR_) * exp1_ - (piG_ / piR_) * exp1_); // A
    case 1: return rate_ * r_ * (piC_ *                exp1_);                             // C
    case 2: return rate_ * r_ * (piG_ * -(piY_ / piR_) * exp1_ + (piG_ / piR_) * exp1_); // G
    case 3: return rate_ * r_ * (piT_ *                exp1_);                             // T, U
    }
  }
  // C
  case 1: {
    switch (j)
    {
    case 0: return rate_ * r_ * (piA_ *                exp1_);                             // A
    case 1: return rate_ * r_ * (piC_ * -(piR_ / piY_) * exp1_ - (piT_ / piY_) * exp1_); // C
    case 2: return rate_ * r_ * (piG_ *                exp1_);                             // G
    case 3: return rate_ * r_ * (piT_ * -(piR_ / piY_) * exp1_ + (piT_ / piY_) * exp1_); // T, U
    }
  }
  // G
  case 2: {
    switch (j)
    {
    case 0: return rate_ * r_ * (piA_ * -(piY_ / piR_) * exp1_ + (piA_ / piR_) * exp1_); // A
    case 1: return rate_ * r_ * (piC_ *                exp1_);                             // C
    case 2: return rate_ * r_ * (piG_ * -(piY_ / piR_) * exp1_ - (piA_ / piR_) * exp1_); // G
    case 3: return rate_ * r_ * (piT_ *                exp1_);                             // T, U
    }
  }
  // T, U
  case 3: {
    switch (j)
    {
    case 0: return rate_ * r_ * (piA_ *                exp1_);                             // A
    case 1: return rate_ * r_ * (piC_ * -(piR_ / piY_) * exp1_ + (piC_ / piY_) * exp1_); // C
    case 2: return rate_ * r_ * (piG_ *                exp1_);                             // G
    case 3: return rate_ * r_ * (piT_ * -(piR_ / piY_) * exp1_ - (piC_ / piY_) * exp1_); // T, U
    }
  }
  }
  return 0;
}

/******************************************************************************/

double F81::d2Pij_dt2(size_t i, size_t j, double d) const
{
  double r_2 = rate_ * rate_ * r_ * r_;
  l_ = rate_ * r_ * d;
  exp1_ = exp(-l_);

  switch (i)
  {
  // A
  case 0: {
    switch (j)
    {
    case 0: return r_2 * (piA_ * (piY_ / piR_) * exp1_ + (piG_ / piR_) * exp1_); // A
    case 1: return r_2 * (piC_ *             -exp1_);                              // C
    case 2: return r_2 * (piG_ * (piY_ / piR_) * exp1_ - (piG_ / piR_) * exp1_); // G
    case 3: return r_2 * (piT_ *             -exp1_);                              // T, U
    }
  }
  // C
  case 1: {
    switch (j)
    {
    case 0: return r_2 * (piA_ *             -exp1_);                              // A
    case 1: return r_2 * (piC_ * (piR_ / piY_) * exp1_ + (piT_ / piY_) * exp1_); // C
    case 2: return r_2 * (piG_ *             -exp1_);                              // G
    case 3: return r_2 * (piT_ * (piR_ / piY_) * exp1_ - (piT_ / piY_) * exp1_); // T, U
    }
  }
  // G
  case 2: {
    switch (j)
    {
    case 0: return r_2 * (piA_ * (piY_ / piR_) * exp1_ - (piA_ / piR_) * exp1_); // A
    case 1: return r_2 * (piC_ *             -exp1_);                              // C
    case 2: return r_2 * (piG_ * (piY_ / piR_) * exp1_ + (piA_ / piR_) * exp1_); // G
    case 3: return r_2 * (piT_ *             -exp1_);                              // T, U
    }
  }
  // T, U
  case 3: {
    switch (j)
    {
    case 0: return r_2 * (piA_ *             -exp1_);                             // A
    case 1: return r_2 * (piC_ * (piR_ / piY_) * exp1_ - (piC_ / piY_) * exp1_); // C
    case 2: return r_2 * (piG_ *             -exp1_);                             // G
    case 3: return r_2 * (piT_ * (piR_ / piY_) * exp1_ + (piC_ / piY_) * exp1_); // T, U
    }
  }
  }
  return 0;
}

/******************************************************************************/

const Matrix<double>& F81::getPij_t(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-l_);

  // A
  p_(0, 0) = piA_ * (1. + (piY_ / piR_) * exp1_) + (piG_ / piR_) * exp1_; // A
  p_(0, 1) = piC_ * (1. -               exp1_);                        // C
  p_(0, 2) = piG_ * (1. + (piY_ / piR_) * exp1_) - (piG_ / piR_) * exp1_; // G
  p_(0, 3) = piT_ * (1. -               exp1_);                        // T, U

  // C
  p_(1, 0) = piA_ * (1. -               exp1_);                        // A
  p_(1, 1) = piC_ * (1. + (piR_ / piY_) * exp1_) + (piT_ / piY_) * exp1_; // C
  p_(1, 2) = piG_ * (1. -               exp1_);                        // G
  p_(1, 3) = piT_ * (1. + (piR_ / piY_) * exp1_) - (piT_ / piY_) * exp1_; // T, U

  // G
  p_(2, 0) = piA_ * (1. + (piY_ / piR_) * exp1_) - (piA_ / piR_) * exp1_; // A
  p_(2, 1) = piC_ * (1. -               exp1_);                        // C
  p_(2, 2) = piG_ * (1. + (piY_ / piR_) * exp1_) + (piA_ / piR_) * exp1_; // G
  p_(2, 3) = piT_ * (1. -               exp1_);                        // T, U

  // T, U
  p_(3, 0) = piA_ * (1. -               exp1_);                        // A
  p_(3, 1) = piC_ * (1. + (piR_ / piY_) * exp1_) - (piC_ / piY_) * exp1_; // C
  p_(3, 2) = piG_ * (1. -               exp1_);                        // G
  p_(3, 3) = piT_ * (1. + (piR_ / piY_) * exp1_) + (piC_ / piY_) * exp1_; // T, U

  return p_;
}

const Matrix<double>& F81::getdPij_dt(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-l_);

  // A
  p_(0, 0) = rate_ * r_ * (piA_ * -(piY_ / piR_) * exp1_ - (piG_ / piR_) * exp1_); // A
  p_(0, 1) = rate_ * r_ * (piC_ *                exp1_);                              // C
  p_(0, 2) = rate_ * r_ * (piG_ * -(piY_ / piR_) * exp1_ + (piG_ / piR_) * exp1_); // G
  p_(0, 3) = rate_ * r_ * (piT_ *                exp1_);                              // T, U

  // C
  p_(1, 0) = rate_ * r_ * (piA_ *                exp1_);                              // A
  p_(1, 1) = rate_ * r_ * (piC_ * -(piR_ / piY_) * exp1_ - (piT_ / piY_) * exp1_); // C
  p_(1, 2) = rate_ * r_ * (piG_ *                exp1_);                              // G
  p_(1, 3) = rate_ * r_ * (piT_ * -(piR_ / piY_) * exp1_ + (piT_ / piY_) * exp1_); // T, U

  // G
  p_(2, 0) = rate_ * r_ * (piA_ * -(piY_ / piR_) * exp1_ + (piA_ / piR_) * exp1_); // A
  p_(2, 1) = rate_ * r_ * (piC_ *                exp1_);                              // C
  p_(2, 2) = rate_ * r_ * (piG_ * -(piY_ / piR_) * exp1_ - (piA_ / piR_) * exp1_); // G
  p_(2, 3) = rate_ * r_ * (piT_ *                exp1_);                              // T, U

  // T, U
  p_(3, 0) = rate_ * r_ * (piA_ *                exp1_);                              // A
  p_(3, 1) = rate_ * r_ * (piC_ * -(piR_ / piY_) * exp1_ + (piC_ / piY_) * exp1_); // C
  p_(3, 2) = rate_ * r_ * (piG_ *                exp1_);                              // G
  p_(3, 3) = rate_ * r_ * (piT_ * -(piR_ / piY_) * exp1_ - (piC_ / piY_) * exp1_); // T, U

  return p_;
}

const Matrix<double>& F81::getd2Pij_dt2(double d) const
{
  double r_2 = rate_ * rate_ * r_ * r_;
  l_ = rate_ * r_ * d;
  exp1_ = exp(-l_);

  // A
  p_(0, 0) = r_2 * (piA_ * (piY_ / piR_) * exp1_ + (piG_ / piR_) * exp1_); // A
  p_(0, 1) = r_2 * (piC_ *             -exp1_);                              // C
  p_(0, 2) = r_2 * (piG_ * (piY_ / piR_) * exp1_ - (piG_ / piR_) * exp1_); // G
  p_(0, 3) = r_2 * (piT_ *             -exp1_);                              // T, U

  // C
  p_(1, 0) = r_2 * (piA_ *             -exp1_);                              // A
  p_(1, 1) = r_2 * (piC_ * (piR_ / piY_) * exp1_ + (piT_ / piY_) * exp1_); // C
  p_(1, 2) = r_2 * (piG_ *             -exp1_);                              // G
  p_(1, 3) = r_2 * (piT_ * (piR_ / piY_) * exp1_ - (piT_ / piY_) * exp1_); // T, U

  // G
  p_(2, 0) = r_2 * (piA_ * (piY_ / piR_) * exp1_ - (piA_ / piR_) * exp1_); // A
  p_(2, 1) = r_2 * (piC_ *             -exp1_);                              // C
  p_(2, 2) = r_2 * (piG_ * (piY_ / piR_) * exp1_ + (piA_ / piR_) * exp1_); // G
  p_(2, 3) = r_2 * (piT_ *             -exp1_);                              // T, U

  // T, U
  p_(3, 0) = r_2 * (piA_ *             -exp1_);                              // A
  p_(3, 1) = r_2 * (piC_ * (piR_ / piY_) * exp1_ - (piC_ / piY_) * exp1_); // C
  p_(3, 2) = r_2 * (piG_ *             -exp1_);                              // G
  p_(3, 3) = r_2 * (piT_ * (piR_ / piY_) * exp1_ + (piC_ / piY_) * exp1_); // T, U

  return p_;
}

/******************************************************************************/

void F81::setFreq(std::map<int, double>& freqs)
{
  piA_ = freqs[0];
  piC_ = freqs[1];
  piG_ = freqs[2];
  piT_ = freqs[3];
  vector<string> thetas(3);
  thetas[0] = getNamespace() + "theta";
  thetas[1] = getNamespace() + "theta1";
  thetas[2] = getNamespace() + "theta2";
  ParameterList pl = getParameters().createSubList(thetas);
  pl[0].setValue(piC_ + piG_);
  pl[1].setValue(piA_ / (piA_ + piT_));
  pl[2].setValue(piG_ / (piC_ + piG_));
  setParametersValues(pl);
}

/******************************************************************************/
