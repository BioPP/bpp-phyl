// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "RN95s.h"

// From the STL:
#include <cmath>

using namespace bpp;
using namespace std;

/******************************************************************************/

RN95s::RN95s(
    shared_ptr<const NucleicAlphabet> alphabet,
    double alpha,
    double beta,
    double gamma,
    double delta) :
  AbstractParameterAliasable("RN95s."),
  AbstractNucleotideSubstitutionModel(alphabet, make_shared<CanonicalStateMap>(alphabet, false), "RN95s."),
  alpha_(),
  beta_(),
  gamma_(),
  delta_()
{
  double f = alpha + beta + gamma + delta;

  alpha_ = alpha / f;
  beta_  = beta  / f;
  gamma_ = gamma / f;
  delta_ = delta / f;

  addParameter_(new Parameter("RN95s.theta", alpha_+gamma_, std::make_shared<IntervalConstraint>(0.001, 0.999, true, true)));
  addParameter_(new Parameter("RN95s.alphaP", alpha_/(alpha_+gamma_), Parameter::PROP_CONSTRAINT_EX));
  addParameter_(new Parameter("RN95s.betaP", beta_/(beta_+delta_), Parameter::PROP_CONSTRAINT_EX));

  computeFrequencies(false);
  updateMatrices_();
}

/******************************************************************************/

void RN95s::updateMatrices_()
{
  double theta  = getParameterValue("theta");
  double alphaP  = getParameterValue("alphaP");
  double betaP  = getParameterValue("betaP");

  alpha_  = alphaP * theta;
  gamma_  = theta - alpha_;
  beta_   = betaP * (1-theta);
  delta_ = 1 - theta - beta_;
  
  // stationnary frequencies

  freq_[1] = freq_[2] = theta/2;
  freq_[3] = freq_[0] = (1- theta)/2;

  // Generator matrix:

  generator_(0, 1) = gamma_;
  generator_(0, 2) = alpha_;
  generator_(0, 3) = delta_;

  generator_(0, 0) = -(gamma_ + alpha_ + delta_);

  generator_(1, 0) = delta_;
  generator_(1, 2) = gamma_;
  generator_(1, 3) = beta_;

  generator_(1, 1) = -(delta_ + beta_ + gamma_);

  generator_(2, 0) = beta_;
  generator_(2, 1) = gamma_;
  generator_(2, 3) = delta_;

  generator_(2, 2) = -(gamma_ + beta_ + delta_);

  generator_(3, 0) = delta_;
  generator_(3, 1) = alpha_;
  generator_(3, 2) = gamma_;

  generator_(3, 3) = -(delta_ + alpha_ + gamma_);

  // Normalization

  double x = 0;
  for (unsigned int i = 0; i < 4; i++)
    x += generator_(i, i) * freq_[i];

  double r_ = isScalable() ? -1 / x : 1;

  setScale(r_);

  // variables for calculation purposes

  double c_1 = gamma_ + delta_ - beta_  - alpha_ ;
  double c_2 = delta_ * gamma_ - alpha_ * beta_;
  double c_3 = beta_ * gamma_ - alpha_ * delta_;

  // eigen vectors and values

  eigenValues_[0] = - 2 * (gamma_ + delta_) * r_;
  eigenValues_[1] = -r_;
  eigenValues_[2] = -r_;
  eigenValues_[3] = 0;

  rightEigenVectors_(0, 0) = 1.;
  rightEigenVectors_(1, 0) = -1.;
  rightEigenVectors_(2, 0) = 1.;
  rightEigenVectors_(3, 0) = -1.;

  rightEigenVectors_(0, 1) = c_2;
  rightEigenVectors_(1, 1) = 0;
  rightEigenVectors_(2, 1) = (beta_ - delta_) * (delta_ + beta_);
  rightEigenVectors_(3, 1) = -c_3;

  rightEigenVectors_(0, 2) = 0;
  rightEigenVectors_(1, 2) = c_2;
  rightEigenVectors_(2, 2) = c_3;
  rightEigenVectors_(3, 2) = (alpha_ - gamma_) * (alpha_ + gamma_);

  rightEigenVectors_(0, 3) = 1;
  rightEigenVectors_(1, 3) = 1;
  rightEigenVectors_(2, 3) = 1;
  rightEigenVectors_(3, 3) = 1;

// Need formula

  if (abs(c_1) < NumConstants::TINY() ||  abs(c_2) < NumConstants::TINY())
  {
    
    ApplicationTools::displayMessage("Singularity during diagonalization of RN95s. Taylor series used instead.");
    isNonSingular_ = false;
    isDiagonalizable_ = false;
    MatrixTools::Taylor(generator_, 30, vPowGen_);
  }
  else
  {
    isNonSingular_ = true;
    isDiagonalizable_ = true;

    leftEigenVectors_(0, 0) = (delta_ - beta_)/(2 * c_1);
    leftEigenVectors_(0, 1) = (alpha_ - gamma_)/(2 * c_1);
    leftEigenVectors_(0, 2) = - leftEigenVectors_(0, 1);
    leftEigenVectors_(0, 3) = - leftEigenVectors_(0, 0);

    leftEigenVectors_(1, 0) = (gamma_ * (gamma_ + delta_)  - alpha_ * (alpha_ + beta_))/(c_1 * c_2);
    leftEigenVectors_(1, 1) = c_3 / (c_1 * c_2);
    leftEigenVectors_(1, 2) = -leftEigenVectors_(1,0);
    leftEigenVectors_(1, 3) = -leftEigenVectors_(1,1);

    leftEigenVectors_(2, 0) = leftEigenVectors_(1,3);
    leftEigenVectors_(2, 1) = (delta_ * (gamma_ + delta_) - beta_ * (alpha_ + beta_))/(c_1 * c_2);
    leftEigenVectors_(2, 2) = leftEigenVectors_(1,1);
    leftEigenVectors_(2, 3) = -leftEigenVectors_(2,1);

    leftEigenVectors_(3, 0) = (1-theta)/2;
    leftEigenVectors_(3, 1) = theta/2;
    leftEigenVectors_(3, 2) = leftEigenVectors_(3,1);
    leftEigenVectors_(3, 3) = leftEigenVectors_(3,0);
  }

  // and the exchangeability_
  for (size_t i = 0; i < size_; i++)
  {
    for (size_t j = 0; j < size_; j++)
    {
      exchangeability_(i, j) = generator_(i, j) / freq_[j];
    }
  }

}


/******************************************************************************/
void RN95s::setFreq(map<int, double>& freqs)
{
  setParameterValue("thetaA", (freqs[0] + freqs[3]) / 2);
  updateMatrices_();
}

/******************************************************************************/
