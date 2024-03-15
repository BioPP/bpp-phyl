// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "RN95.h"

// From the STL:
#include <cmath>

using namespace bpp;
using namespace std;

/******************************************************************************/

RN95::RN95(
  shared_ptr<const NucleicAlphabet> alphabet,
  double alpha,
  double beta,
  double gamma,
  double delta,
  double epsilon,
  double kappa,
  double lambda,
  double sigma) :
  AbstractParameterAliasable("RN95."),
  AbstractNucleotideSubstitutionModel(alphabet, make_shared<CanonicalStateMap>(alphabet, false), "RN95."),
  alpha_(),
  beta_(),
  gamma_(),
  delta_(),
  epsilon_(),
  kappa_(),
  lambda_(),
  sigma_()
{
  double f = gamma + lambda + delta + kappa;

  alpha_   = alpha / f;
  beta_    = beta / f;
  gamma_   = gamma / f;
  delta_   = delta / f;
  epsilon_ = epsilon / f;
  kappa_   = kappa / f;
  lambda_  = lambda / f;
  sigma_   = sigma / f;

  double thetaR = delta_ + kappa_;
  double kappaP = kappa_ / thetaR;
  double gammaP = gamma_ / (1 - thetaR);

  addParameter_(new Parameter("RN95.thetaR", thetaR, Parameter::PROP_CONSTRAINT_EX));
  addParameter_(new Parameter("RN95.gammaP", gammaP, Parameter::PROP_CONSTRAINT_EX));
  addParameter_(new Parameter("RN95.kappaP", kappaP, Parameter::PROP_CONSTRAINT_EX));

  addParameter_(new Parameter("RN95.alpha"  , alpha_  , Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("RN95.sigma"  , sigma_  , Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("RN95.beta"   , beta_   , Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("RN95.epsilon", epsilon_, Parameter::R_PLUS_STAR));

  computeFrequencies(false);
  updateMatrices_();
}

/******************************************************************************/

void RN95::updateMatrices_()
{
  double thetaR  = getParameterValue("thetaR");
  double gammaP  = getParameterValue("gammaP");
  double kappaP  = getParameterValue("kappaP");

  alpha_ = getParameterValue("alpha");
  sigma_ = getParameterValue("sigma");
  beta_ = getParameterValue("beta");
  epsilon_ = getParameterValue("epsilon");

  // thetaR = delta_ + kappa_
  kappa_ = kappaP * thetaR;
  gamma_ = gammaP * (1 - thetaR);
  delta_ = thetaR - kappa_;
  lambda_ = 1 - thetaR - gamma_;

  // variables for calculation purposes
  double thetaY = 1 - thetaR;

  auto c_1 = thetaR + sigma_ + beta_;
  auto c_2 = 1 - c_1;
  auto c_3 = thetaY + alpha_ + epsilon_;
  auto c_4 = 1 - c_3;

  // stationnary frequencies

  freq_[0] = (delta_ * thetaY + epsilon_ * thetaR)/c_3;
  freq_[1] = (sigma_ * thetaY + gamma_   * thetaR)/c_1;
  freq_[2] = (kappa_ * thetaY + alpha_   * thetaR)/c_3;
  freq_[3] = (beta_  * thetaY + lambda_  * thetaR)/c_1;
  
  // Generator matrix:

  generator_(0, 1) = gamma_;
  generator_(0, 2) = alpha_;
  generator_(0, 3) = lambda_;

  generator_(0, 0) = -(gamma_ + alpha_ + lambda_);

  generator_(1, 0) = delta_;
  generator_(1, 2) = kappa_;
  generator_(1, 3) = beta_;

  generator_(1, 1) = -(delta_ + beta_ + kappa_);

  generator_(2, 0) = epsilon_;
  generator_(2, 1) = gamma_;
  generator_(2, 3) = lambda_;

  generator_(2, 2) = -(gamma_ + epsilon_ + lambda_);

  generator_(3, 0) = delta_;
  generator_(3, 1) = sigma_;
  generator_(3, 2) = kappa_;

  generator_(3, 3) = -(delta_ + sigma_ + kappa_);

  // Normalization

  double x = 0;
  for (size_t i = 0; i < 4; i++)
    x += generator_(i, i) * freq_[i];

  auto r_ = isScalable() ? -1 / x : 1;

  setScale(r_);


  // eigen vectors and values

  eigenValues_[0] = -r_;
  eigenValues_[1] = -c_3 * r_;
  eigenValues_[2] = -c_1 * r_;
  eigenValues_[3] = 0;

  rightEigenVectors_(0, 0) = thetaY;
  rightEigenVectors_(1, 0) = -thetaR;
  rightEigenVectors_(2, 0) = thetaY;
  rightEigenVectors_(3, 0) = -thetaR;

  rightEigenVectors_(0, 1) = (kappa_ - alpha_) * thetaY + alpha_ * c_4;
  rightEigenVectors_(1, 1) = alpha_ * delta_ - epsilon_ * kappa_;
  rightEigenVectors_(2, 1) = (epsilon_ - delta_) * thetaY - epsilon_ *  c_4;
  rightEigenVectors_(3, 1) = alpha_ * delta_ - epsilon_ * kappa_;

  rightEigenVectors_(0, 2) = sigma_ * lambda_ - beta_ * gamma_;
  rightEigenVectors_(1, 2) = (beta_ - lambda_) * thetaR - beta_ * c_2;
  rightEigenVectors_(2, 2) = sigma_ * lambda_ - beta_ * gamma_; 
  rightEigenVectors_(3, 2) = (gamma_ - sigma_) * thetaR + sigma_ *  c_2;

  rightEigenVectors_(0, 3) = 1.;
  rightEigenVectors_(1, 3) = 1.;
  rightEigenVectors_(2, 3) = 1.;
  rightEigenVectors_(3, 3) = 1.;

  // Need formula

  if (abs(c_2) < NumConstants::TINY() || abs(c_4) < NumConstants::TINY())
  {
    ApplicationTools::displayMessage("Singularity during diagonalization of RN95. Taylor series used instead.");
    
    isNonSingular_ = false;
    isDiagonalizable_ = false;
    MatrixTools::Taylor(generator_, 30, vPowGen_);
  }
  else
  {
    isNonSingular_ = true;
    isDiagonalizable_ = true;

    leftEigenVectors_(0, 0) = (delta_ - epsilon_)/c_4;
    leftEigenVectors_(0, 1) = (sigma_ - gamma_)/c_2;
    leftEigenVectors_(0, 2) = (kappa_ - alpha_)/c_4;
    leftEigenVectors_(0, 3) = (beta_ - lambda_)/c_2;

    leftEigenVectors_(1, 0) = 1/(c_3 * c_4);
    leftEigenVectors_(1, 1) = 0;
    leftEigenVectors_(1, 2) = -1/(c_3 * c_4);
    leftEigenVectors_(1, 3) = 0;

    leftEigenVectors_(2, 0) = 0;
    leftEigenVectors_(2, 1) = -1/(c_1 * c_2);
    leftEigenVectors_(2, 2) = 0;
    leftEigenVectors_(2, 3) = 1/(c_1 * c_2);

    leftEigenVectors_(3, 0) = (epsilon_ + (delta_ - epsilon_) * thetaY)/c_3;
    leftEigenVectors_(3, 1) = (gamma_ + (sigma_ - gamma_) * thetaY)/c_1;
    leftEigenVectors_(3, 2) = (alpha_ + (kappa_ - alpha_) * thetaY)/c_3;
    leftEigenVectors_(3, 3) = (lambda_ + (beta_ - lambda_) * thetaY)/c_1;
  }

  // and the exchangeability_
  for (unsigned int i = 0; i < size_; ++i)
  {
    for (unsigned int j = 0; j < size_; ++j)
    {
      exchangeability_(i, j) = generator_(i, j) / freq_[j];
    }
  }
}


/******************************************************************************/
void RN95::setFreq(map<int, double>& freqs)
{
  auto thetaR = (freqs[0] + freqs[2]) / (freqs[0] + freqs[1] + freqs[2] + freqs[3]);
  setParameterValue("thetaR", thetaR);
  setParameterValue("gammaP", gamma_ / (1 - thetaR));
  setParameterValue("kappaP", kappa_ / thetaR);
  
  updateMatrices_();
}

/******************************************************************************/
