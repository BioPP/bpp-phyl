// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "../FrequencySet/NucleotideFrequencySet.h"
#include "GTR.h"

// From the STL:
#include <cmath>

using namespace bpp;
using namespace std;

/******************************************************************************/

GTR::GTR(
    shared_ptr<const NucleicAlphabet> alpha,
    double a,
    double b,
    double c,
    double d,
    double e,
    double piA,
    double piC,
    double piG,
    double piT) :
  AbstractParameterAliasable("GTR."),
  AbstractReversibleNucleotideSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "GTR."),
  a_(a), b_(b), c_(c), d_(d), e_(e), piA_(piA), piC_(piC), piG_(piG), piT_(piT), theta_(piG + piC), theta1_(piA / (1. - theta_)), theta2_(piG / theta_), p_()
{
  addParameter_(new Parameter("GTR.a", a, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("GTR.b", b, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("GTR.c", c, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("GTR.d", d, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("GTR.e", e, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("GTR.theta", theta_, FrequencySetInterface::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter("GTR.theta1", theta1_, FrequencySetInterface::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter("GTR.theta2", theta2_, FrequencySetInterface::FREQUENCE_CONSTRAINT_SMALL));
  updateMatrices_();
}

/******************************************************************************/

void GTR::updateMatrices_()
{
  a_ = getParameterValue("a");
  b_ = getParameterValue("b");
  c_ = getParameterValue("c");
  d_ = getParameterValue("d");
  e_ = getParameterValue("e");
  theta_  = getParameterValue("theta");
  theta1_ = getParameterValue("theta1");
  theta2_ = getParameterValue("theta2");
  piA_ = theta1_ * (1. - theta_);
  piC_ = (1. - theta2_) * theta_;
  piG_ = theta2_ * theta_;
  piT_ = (1. - theta1_) * (1. - theta_);
  p_ = 2 * (a_ * piC_ * piT_ + b_ * piA_ * piT_ + c_ * piG_ * piT_ + d_ * piA_ * piC_ + e_ * piC_ * piG_ + piA_ * piG_);

  freq_[0] = piA_;
  freq_[1] = piC_;
  freq_[2] = piG_;
  freq_[3] = piT_;

  // Exchangeability matrix:
  exchangeability_(0, 0) = (-b_ * piT_ - piG_ - d_ * piC_) / (piA_ * p_);
  exchangeability_(1, 0) = d_ / p_;
  exchangeability_(0, 1) = d_ / p_;
  exchangeability_(2, 0) = 1 / p_;
  exchangeability_(0, 2) = 1 / p_;
  exchangeability_(3, 0) = b_ / p_;
  exchangeability_(0, 3) = b_ / p_;
  exchangeability_(1, 1) = (-a_ * piT_ - e_ * piG_ - d_ * piA_) / (piC_ * p_);
  exchangeability_(1, 2) = e_ / p_;
  exchangeability_(2, 1) = e_ / p_;
  exchangeability_(1, 3) = a_ / p_;
  exchangeability_(3, 1) = a_ / p_;
  exchangeability_(2, 2) = (-c_ * piT_ - e_ * piC_ - piA_) / (piG_ * p_);
  exchangeability_(2, 3) = c_ / p_;
  exchangeability_(3, 2) = c_ / p_;
  exchangeability_(3, 3) = (-c_ * piG_ - a_ * piC_ - b_ * piA_) / (piT_ * p_);

  AbstractReversibleSubstitutionModel::updateMatrices_();
}

/******************************************************************************/

void GTR::setFreq(map<int, double>& freqs)
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
