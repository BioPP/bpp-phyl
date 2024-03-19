// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "../FrequencySet/NucleotideFrequencySet.h"
#include "L95.h"

// From the STL:
#include <cmath>

using namespace bpp;
using namespace std;

/******************************************************************************/

L95::L95(
    shared_ptr<const NucleicAlphabet> alphabet,
    double alpha, double beta, double gamma, double kappa, double theta) :
  AbstractParameterAliasable("L95."),
  AbstractNucleotideSubstitutionModel(
    alphabet,
    make_shared<CanonicalStateMap>(alphabet, false),
    "L95."),
  alpha_(alpha), beta_(beta), gamma_(gamma), kappa_(kappa), theta_(theta)
{
  addParameter_(new Parameter("L95.alpha", alpha, Parameter::PROP_CONSTRAINT_IN));
  addParameter_(new Parameter("L95.beta", beta, Parameter::PROP_CONSTRAINT_IN));
  addParameter_(new Parameter("L95.gamma", gamma, Parameter::PROP_CONSTRAINT_IN));
  addParameter_(new Parameter("L95.kappa", kappa, make_shared<IntervalConstraint>(0, 1000, false, false, NumConstants::MILLI())));
  addParameter_(new Parameter("L95.theta", theta, make_shared<IntervalConstraint>(0, 1, false, false, NumConstants::MILLI())));

  computeFrequencies(false);
  updateMatrices_();
}

/******************************************************************************/

void L95::updateMatrices_()
{
  alpha_  = getParameterValue("alpha");
  beta_   = getParameterValue("beta");
  gamma_  = getParameterValue("gamma");
  kappa_  = getParameterValue("kappa");
  theta_  = getParameterValue("theta");

  freq_[0] = (1 - theta_) / 2;
  freq_[1] = theta_ / 2;
  freq_[2] = theta_ / 2;
  freq_[3] = (1 - theta_) / 2;

  // Generator matrix:
  generator_(0, 0) = -kappa_ * theta_ - gamma_;
  generator_(0, 1) = kappa_ * beta_ * theta_;
  generator_(0, 2) = kappa_ * (1 - beta_) * theta_;
  generator_(0, 3) = gamma_;
  generator_(1, 0) = kappa_ * alpha_ * ( 1 - theta_);
  generator_(1, 1) = -kappa_ * (1 - theta_) + gamma_ - 1;
  generator_(1, 2) = 1 - gamma_;
  generator_(1, 3) = kappa_ * (1 - theta_) * (1 - alpha_);
  generator_(2, 0) = kappa_ * (1 - theta_) * (1 - alpha_);
  generator_(2, 1) = 1 - gamma_;
  generator_(2, 2) = -kappa_ * (1 - theta_) + gamma_ - 1;
  generator_(2, 3) = kappa_ * alpha_ * (1 - theta_);
  generator_(3, 0) = gamma_;
  generator_(3, 1) = kappa_ * (1 - beta_) * theta_;
  generator_(3, 2) = kappa_ * beta_ * theta_;
  generator_(3, 3) = -kappa_ * theta_ - gamma_;

  setScale(1. / (2 * kappa_ * theta_ * (1 - theta_) + gamma_ + theta_ - 2 * theta_ * gamma_));

  AbstractSubstitutionModel::updateMatrices_();
}

/******************************************************************************/

void L95::setFreq(map<int, double>& freqs)
{
  setParameterValue("theta", freqs[1] + freqs[2]);
  updateMatrices_();
}

/******************************************************************************/
