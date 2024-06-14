// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "SSR.h"

// From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

// From the STL:
#include <cmath>

using namespace bpp;
using namespace std;

/******************************************************************************/

SSR::SSR(
    shared_ptr<const NucleicAlphabet> alpha,
    double beta,
    double gamma,
    double delta,
    double theta) :
  AbstractParameterAliasable("SSR."),
  AbstractReversibleNucleotideSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "SSR."),
  beta_(beta), gamma_(gamma), delta_(delta), theta_(theta),
  piA_((1. - theta) / 2.), piC_(theta / 2.), piG_(theta / 2.), piT_((1. - theta) / 2.)
{
  addParameter_(new Parameter("SSR.beta", beta, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("SSR.gamma", gamma, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("SSR.delta", delta, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("SSR.theta", theta, Parameter::PROP_CONSTRAINT_EX));
  updateMatrices_();
}

/******************************************************************************/

void SSR::updateMatrices_()
{
  beta_  = getParameterValue("beta");
  gamma_ = getParameterValue("gamma");
  delta_ = getParameterValue("delta");
  theta_ = getParameterValue("theta");

  freq_[0] = piA_ = (1. - theta_) / 2.;
  freq_[1] = piC_ = theta_ / 2.;
  freq_[2] = piG_ = theta_ / 2;
  freq_[3] = piT_ = (1. - theta_) / 2.;

  // Exchangeability matrix:
  exchangeability_(0, 0) = -gamma_ * piT_ - piG_ - beta_ * piC_;
  exchangeability_(1, 0) = beta_;
  exchangeability_(0, 1) = beta_;
  exchangeability_(2, 0) = 1.;
  exchangeability_(0, 2) = 1.;
  exchangeability_(3, 0) = gamma_;
  exchangeability_(0, 3) = gamma_;
  exchangeability_(1, 1) = -piT_ - delta_ * piG_ - beta_ * piA_;
  exchangeability_(1, 2) = delta_;
  exchangeability_(2, 1) = delta_;
  exchangeability_(1, 3) = 1.;
  exchangeability_(3, 1) = 1.;
  exchangeability_(2, 2) = -beta_ * piT_ - delta_ * piC_ - piA_;
  exchangeability_(2, 3) = beta_;
  exchangeability_(3, 2) = beta_;
  exchangeability_(3, 3) = -beta_ * piG_ - piC_ - gamma_ * piA_;

  AbstractReversibleSubstitutionModel::updateMatrices_();
}

/******************************************************************************/

void SSR::setFreq(map<int, double>& freqs)
{
  piC_ = freqs[1];
  piG_ = freqs[2];
  setParameterValue("theta", piC_ + piG_);
  updateMatrices_();
}

/******************************************************************************/
