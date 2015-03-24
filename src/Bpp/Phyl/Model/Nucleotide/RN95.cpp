//
// File: RN95.cpp
// Created by: Laurent Guéguen
// Created on: jeudi 24 février 2011, à 20h 42
//

/*
   Copyright or © or Copr. CNRS, (November 16, 2004)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#include "RN95.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

// From the STL:
#include <cmath>

using namespace bpp;
using namespace std;

/******************************************************************************/

RN95::RN95(
  const NucleicAlphabet* alphabet,
  double alpha,
  double beta,
  double gamma,
  double delta,
  double epsilon,
  double kappa,
  double lambda,
  double sigma) :
  AbstractParameterAliasable("RN95."),
  AbstractNucleotideSubstitutionModel(alphabet, new CanonicalStateMap(alphabet, false), "RN95."),
  alpha_(),
  beta_(),
  gamma_(),
  delta_(),
  epsilon_(),
  kappa_(),
  lambda_(),
  sigma_(),
  r_(),
  c1_(),
  c2_(),
  c3_(),
  c4_(),
  c5_(),
  c6_(),
  c7_(),
  c8_(),
  c9_(),
  p_(size_, size_),
  exp1_(),
  exp3_(),
  exp6_(),
  l_()
{
  double f = gamma + lambda + delta + kappa;

  alpha_  = alpha / f;
  beta_   = beta / f;
  gamma_  = gamma / f;
  delta_  = delta / f;
  epsilon_ = epsilon / f;
  kappa_  = kappa / f;
  lambda_ = lambda / f;
  sigma_  = sigma / f;

  double thetaR = delta_ + kappa_;
  double thetaC = (gamma_ * thetaR + sigma_ * (1 - thetaR)) / (beta_ + sigma_ + thetaR) / (1 - thetaR);
  double thetaG = (alpha_ * thetaR + kappa_ * (1 - thetaR)) / (alpha_ + epsilon_ + 1 - thetaR) / thetaR;
  double kappaP = kappa_ / thetaR;
  double gammaP = gamma_ / (1 - thetaR);
  double alphaP = (alpha_ * (1 - thetaG) + (thetaG < kappaP ? thetaG : kappaP) * (1 - thetaR)) / (thetaG * (1 - thetaR));
  double sigmaP = (sigma_ * (1 - thetaC) + (thetaC < gammaP ? thetaC : gammaP) * thetaR) / (thetaC * thetaR);
  addParameter_(new Parameter("RN95.thetaR", thetaR, &Parameter::PROP_CONSTRAINT_EX));
  addParameter_(new Parameter("RN95.thetaC", thetaC, &Parameter::PROP_CONSTRAINT_EX));
  addParameter_(new Parameter("RN95.thetaG", thetaG, &Parameter::PROP_CONSTRAINT_EX));
  addParameter_(new Parameter("RN95.gammaP", gammaP, &Parameter::PROP_CONSTRAINT_EX));
  addParameter_(new Parameter("RN95.kappaP", kappaP, &Parameter::PROP_CONSTRAINT_EX));
  addParameter_(new Parameter("RN95.alphaP", alphaP, new IntervalConstraint(1, 1, false), true));
  addParameter_(new Parameter("RN95.sigmaP", sigmaP, new IntervalConstraint(1, 1, false), true));

  updateMatrices();
}

/******************************************************************************/
void RN95::updateMatrices()
{
  double alphaP  = getParameterValue("alphaP");
  double sigmaP  = getParameterValue("sigmaP");
  double thetaR  = getParameterValue("thetaR");
  double thetaC  = getParameterValue("thetaC");
  double thetaG  = getParameterValue("thetaG");
  double gammaP  = getParameterValue("gammaP");
  double kappaP  = getParameterValue("kappaP");

  kappa_ = kappaP * thetaR;
  gamma_ = gammaP * (1 - thetaR);
  delta_ = thetaR - kappa_;
  lambda_ = 1 - thetaR - gamma_;
  alpha_ = (alphaP * (1 - thetaR) * thetaG - (thetaG < kappaP ? thetaG : kappaP) * (1 - thetaR)) / (1 - thetaG);
  sigma_ = (sigmaP * thetaR * thetaC - (thetaC < gammaP ? thetaC : gammaP) * thetaR) / (1 - thetaC);
  epsilon_ = (alpha_ * thetaR + kappa_ * (1 - thetaR)) / (thetaG * thetaR) - alpha_ - (1 - thetaR);
  beta_ = (gamma_ * thetaR + sigma_ * (1 - thetaR)) / (thetaC * (1 - thetaR)) - sigma_ - thetaR;

  // stationnary frequencies

  freq_[0] = (1 - thetaG) * thetaR;
  freq_[1] = thetaC * (1 - thetaR);
  freq_[2] = thetaG * thetaR;
  freq_[3] = (1 - thetaC) * (1 - thetaR);

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
  {
    x += generator_(i, i) * freq_[i];
  }

  r_ = -1 / x;

  MatrixTools::scale(generator_, r_);
  // variables for calculation purposes

  c1_ = 1;
  c2_ = gamma_ + lambda_;
  c3_ = alpha_ + gamma_ + epsilon_ + lambda_;
  c4_ = kappa_ - alpha_;
  c5_ = delta_ + kappa_;
  c6_ = delta_ + kappa_ + sigma_ + beta_;
  c7_ = gamma_ - sigma_;
  c8_ = delta_ - epsilon_;
  c9_ = lambda_ - beta_;

  // eigen vectors and values

  eigenValues_[0] = 0;
  eigenValues_[1] = -c1_ * r_;
  eigenValues_[2] = -c3_ * r_;
  eigenValues_[3] = -c6_ * r_;

  rightEigenVectors_(0, 0) = 1.;
  rightEigenVectors_(1, 0) = 1.;
  rightEigenVectors_(2, 0) = 1.;
  rightEigenVectors_(3, 0) = 1.;

  rightEigenVectors_(0, 1) = 1.;
  rightEigenVectors_(1, 1) = -c5_ / c2_;
  rightEigenVectors_(2, 1) = 1.;
  rightEigenVectors_(3, 1) = -c5_ / c2_;

  rightEigenVectors_(0, 2) = (alpha_ * (c5_ - c3_) + kappa_ * c2_) / (delta_ * (c3_ - c2_) - epsilon_ * c5_);
  rightEigenVectors_(1, 2) = 1.;
  rightEigenVectors_(2, 2) = (-epsilon_ * (c5_ - c3_) - delta_ * c2_) / (delta_ * (c3_ - c2_) - epsilon_ * c5_);
  rightEigenVectors_(3, 2) = 1.;

  rightEigenVectors_(0, 3) = 1.;
  rightEigenVectors_(1, 3) = (beta_ * (c2_ - c6_) + lambda_ * c5_) / (gamma_ * (c6_ - c5_) - sigma_ * c2_);
  rightEigenVectors_(2, 3) = 1;
  rightEigenVectors_(3, 3) = (-sigma_ * (c2_ - c6_) - gamma_ * c5_) / (gamma_ * (c6_ - c5_) - sigma_ * c2_);

  // Need formula

  if (enableEigenDecomposition())
  {
    try
    {
      MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);
      isNonSingular_ = true;
      isDiagonalizable_ = true;
      for (size_t i = 0; i < size_ && isDiagonalizable_; i++)
      {
        if (abs(iEigenValues_[i]) > NumConstants::TINY())
          isDiagonalizable_ = false;
      }
    }
    catch (ZeroDivisionException& e)
    {
      ApplicationTools::displayMessage("Singularity during diagonalization of RN95. Taylor series used instead.");
      
      isNonSingular_ = false;
      isDiagonalizable_ = false;
      MatrixTools::Taylor(generator_, 30, vPowGen_);
    }
    
    // and the exchangeability_
    for (unsigned int i = 0; i < size_; i++)
      for (unsigned int j = 0; j < size_; j++)
        exchangeability_(i,j) = generator_(i,j) / freq_[j];
  }
}

/******************************************************************************/
double RN95::Pij_t(size_t i, size_t j, double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-c1_ * l_);
  exp3_ = exp(-c3_ * l_);
  exp6_ = exp(-c6_ * l_);

  switch (i)
  {
    {
    // A
    case 0: {
      switch (j)
      {
      case 0: return freq_[0] - c2_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (alpha_ * (c3_ - c1_) - c2_ * c4_) / (c3_ * (c3_ - c1_)) * exp3_; // A
      case 1: return freq_[1] + c2_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
      case 2: return freq_[2] - c2_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (alpha_ * (c3_ - c1_) - c2_ * c4_) / (c3_ * (c3_ - c1_)) * exp3_;  // G
      case 3: return freq_[3] + c2_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;           // T, U
      }
    }
    // C
    case 1: {
      switch (j)
      {
      case 0: return freq_[0] + c5_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // A
      case 1: return freq_[1] - c5_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (beta_ * (c6_ - c1_) - c5_ * c9_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
      case 2: return freq_[2] + c5_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // G
      case 3: return freq_[3] - c5_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (beta_ * (c6_ - c1_) - c5_ * c9_) / (c6_ * (c6_ - c1_)) * exp6_;                  // T
      }
    }
    // G
    case 2: {
      switch (j)
      {
      case 0: return freq_[0] - c2_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (c2_ * c8_ - epsilon_ * (c3_ - c1_)) / (c3_ * (c3_ - c1_)) * exp3_; // A
      case 1: return freq_[1] + c2_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
      case 2: return freq_[2] - c2_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (c2_ * c8_ - epsilon_ * (c3_ - c1_)) / (c3_ * (c3_ - c1_)) * exp3_;   // G
      case 3: return freq_[3] + c2_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;           // T, U
      }
    }
    }
  // T, U
  case 3: {
    switch (j)
    {
    case 0: return freq_[0] + c5_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // A
    case 1: return freq_[1] - c5_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (c5_ * c7_ - sigma_ * (c6_ - c1_)) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
    case 2: return freq_[2] + c5_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // G
    case 3: return freq_[3] - c5_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (c5_ * c7_ - sigma_ * (c6_ - c1_)) / (c6_ * (c6_ - c1_)) * exp6_;                  // T
    }
  }
  }
  return 0;
}

/******************************************************************************/
double RN95::dPij_dt(size_t i, size_t j, double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = -c1_* rate_* r_* exp(-c1_ * l_);
  exp3_ = -c3_* rate_* r_* exp(-c3_ * l_);
  exp6_ = -c6_* rate_* r_* exp(-c6_ * l_);

  switch (i)
  {
    {
    // A
    case 0: {
      switch (j)
      {
      case 0: return -c2_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (alpha_ * (c3_ - c1_) - c2_ * c4_) / (c3_ * (c3_ - c1_)) * exp3_; // A
      case 1: return c2_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
      case 2: return -c2_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (alpha_ * (c3_ - c1_) - c2_ * c4_) / (c3_ * (c3_ - c1_)) * exp3_;  // G
      case 3: return c2_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;           // T, U
      }
    }
    // C
    case 1: {
      switch (j)
      {
      case 0: return c5_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // A
      case 1: return -c5_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (beta_ * (c6_ - c1_) - c5_ * c9_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
      case 2: return c5_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // G
      case 3: return -c5_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (beta_ * (c6_ - c1_) - c5_ * c9_) / (c6_ * (c6_ - c1_)) * exp6_;                  // T
      }
    }
    // G
    case 2: {
      switch (j)
      {
      case 0: return -c2_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (c2_ * c8_ - epsilon_ * (c3_ - c1_)) / (c3_ * (c3_ - c1_)) * exp3_; // A
      case 1: return c2_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
      case 2: return -c2_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (c2_ * c8_ - epsilon_ * (c3_ - c1_)) / (c3_ * (c3_ - c1_)) * exp3_;   // G
      case 3: return c2_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;           // T, U
      }
    }
    }
  // T, U
  case 3: {
    switch (j)
    {
    case 0: return c5_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // A
    case 1: return -c5_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (c5_ * c7_ - sigma_ * (c6_ - c1_)) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
    case 2: return c5_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // G
    case 3: return -c5_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (c5_ * c7_ - sigma_ * (c6_ - c1_)) / (c6_ * (c6_ - c1_)) * exp6_;                  // T
    }
  }
  }
  return 0;
}

/******************************************************************************/
double RN95::d2Pij_dt2(size_t i, size_t j, double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = c1_ * rate_ * r_ * c1_ * rate_ * r_ * exp(-c1_ * l_);
  exp3_ = c3_ * rate_ * r_ * c3_ * rate_ * r_ * exp(-c3_ * l_);
  exp6_ = c6_ * rate_ * r_ * c6_ * rate_ * r_ * exp(-c6_ * l_);

  switch (i)
  {
    {
    // A
    case 0: {
      switch (j)
      {
      case 0: return -c2_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (alpha_ * (c3_ - c1_) - c2_ * c4_) / (c3_ * (c3_ - c1_)) * exp3_; // A
      case 1: return c2_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
      case 2: return -c2_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (alpha_ * (c3_ - c1_) - c2_ * c4_) / (c3_ * (c3_ - c1_)) * exp3_;  // G
      case 3: return c2_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;           // T, U
      }
    }
    // C
    case 1: {
      switch (j)
      {
      case 0: return c5_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // A
      case 1: return -c5_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (beta_ * (c6_ - c1_) - c5_ * c9_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
      case 2: return c5_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // G
      case 3: return -c5_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (beta_ * (c6_ - c1_) - c5_ * c9_) / (c6_ * (c6_ - c1_)) * exp6_;                  // T
      }
    }
    // G
    case 2: {
      switch (j)
      {
      case 0: return -c2_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (c2_ * c8_ - epsilon_ * (c3_ - c1_)) / (c3_ * (c3_ - c1_)) * exp3_; // A
      case 1: return c2_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
      case 2: return -c2_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (c2_ * c8_ - epsilon_ * (c3_ - c1_)) / (c3_ * (c3_ - c1_)) * exp3_;   // G
      case 3: return c2_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;           // T, U
      }
    }
    }
  // T, U
  case 3: {
    switch (j)
    {
    case 0: return c5_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // A
    case 1: return -c5_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (c5_ * c7_ - sigma_ * (c6_ - c1_)) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
    case 2: return c5_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // G
    case 3: return -c5_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (c5_ * c7_ - sigma_ * (c6_ - c1_)) / (c6_ * (c6_ - c1_)) * exp6_;                  // T
    }
  }
  }
  return 0;
}

/******************************************************************************/

const Matrix<double>& RN95::getPij_t(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-c1_ * l_);
  exp3_ = exp(-c3_ * l_);
  exp6_ = exp(-c6_ * l_);

  // A
  p_(0, 0) = freq_[0] - c2_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (alpha_ * (c3_ - c1_) - c2_ * c4_) / (c3_ * (c3_ - c1_)) * exp3_; // A
  p_(0, 1) = freq_[1] + c2_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
  p_(0, 2) = freq_[2] - c2_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (alpha_ * (c3_ - c1_) - c2_ * c4_) / (c3_ * (c3_ - c1_)) * exp3_;  // G
  p_(0, 3) = freq_[3] + c2_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;           // T, U
  // C
  p_(1, 0) = freq_[0] + c5_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // A
  p_(1, 1) = freq_[1] - c5_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (beta_ * (c6_ - c1_) - c5_ * c9_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
  p_(1, 2) = freq_[2] + c5_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // G
  p_(1, 3) = freq_[3] - c5_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (beta_ * (c6_ - c1_) - c5_ * c9_) / (c6_ * (c6_ - c1_)) * exp6_;                  // T
  // G
  p_(2, 0) = freq_[0] - c2_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (c2_ * c8_ - epsilon_ * (c3_ - c1_)) / (c3_ * (c3_ - c1_)) * exp3_; // A
  p_(2, 1) = freq_[1] + c2_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
  p_(2, 2) = freq_[2] - c2_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (c2_ * c8_ - epsilon_ * (c3_ - c1_)) / (c3_ * (c3_ - c1_)) * exp3_;   // G
  p_(2, 3) = freq_[3] + c2_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;           // T, U
  // T, U
  p_(3, 0) = freq_[0] + c5_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // A
  p_(3, 1) = freq_[1] - c5_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (c5_ * c7_ - sigma_ * (c6_ - c1_)) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
  p_(3, 2) = freq_[2] + c5_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // G
  p_(3, 3) = freq_[3] - c5_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (c5_ * c7_ - sigma_ * (c6_ - c1_)) / (c6_ * (c6_ - c1_)) * exp6_;                  // T

  return p_;
}

/******************************************************************************/

const Matrix<double>&  RN95::getdPij_dt(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = -c1_* rate_* r_* exp(-c1_ * l_);
  exp3_ = -c3_* rate_* r_* exp(-c3_ * l_);
  exp6_ = -c6_* rate_* r_* exp(-c6_ * l_);

  // A
  p_(0, 0) = -c2_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (alpha_ * (c3_ - c1_) - c2_ * c4_) / (c3_ * (c3_ - c1_)) * exp3_; // A
  p_(0, 1) =  c2_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
  p_(0, 2) =  -c2_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (alpha_ * (c3_ - c1_) - c2_ * c4_) / (c3_ * (c3_ - c1_)) * exp3_;  // G
  p_(0, 3) = c2_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;           // T, U
  // C
  p_(1, 0) = c5_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // A
  p_(1, 1) = -c5_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (beta_ * (c6_ - c1_) - c5_ * c9_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
  p_(1, 2) = c5_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // G
  p_(1, 3) = -c5_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (beta_ * (c6_ - c1_) - c5_ * c9_) / (c6_ * (c6_ - c1_)) * exp6_;                  // T
  // G
  p_(2, 0) = -c2_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (c2_ * c8_ - epsilon_ * (c3_ - c1_)) / (c3_ * (c3_ - c1_)) * exp3_; // A
  p_(2, 1) = c2_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
  p_(2, 2) = -c2_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (c2_ * c8_ - epsilon_ * (c3_ - c1_)) / (c3_ * (c3_ - c1_)) * exp3_;   // G
  p_(2, 3) = c2_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;           // T, U
  // T, U
  p_(3, 0) = c5_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // A
  p_(3, 1) = -c5_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (c5_ * c7_ - sigma_ * (c6_ - c1_)) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
  p_(3, 2) = c5_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // G
  p_(3, 3) = -c5_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (c5_ * c7_ - sigma_ * (c6_ - c1_)) / (c6_ * (c6_ - c1_)) * exp6_;                  // T
  return p_;
}

/******************************************************************************/

const Matrix<double>&  RN95::getd2Pij_dt2(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = c1_ * rate_ * r_ * c1_ * rate_ * r_ * exp(-c1_ * l_);
  exp3_ = c3_ * rate_ * r_ * c3_ * rate_ * r_ * exp(-c3_ * l_);
  exp6_ = c6_ * rate_ * r_ * c6_ * rate_ * r_ * exp(-c6_ * l_);

  // A
  p_(0, 0) = -c2_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (alpha_ * (c3_ - c1_) - c2_ * c4_) / (c3_ * (c3_ - c1_)) * exp3_; // A
  p_(0, 1) =  c2_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
  p_(0, 2) =  -c2_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (alpha_ * (c3_ - c1_) - c2_ * c4_) / (c3_ * (c3_ - c1_)) * exp3_;  // G
  p_(0, 3) = c2_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;           // T, U
  // C
  p_(1, 0) = c5_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // A
  p_(1, 1) = -c5_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (beta_ * (c6_ - c1_) - c5_ * c9_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
  p_(1, 2) = c5_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // G
  p_(1, 3) = -c5_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (beta_ * (c6_ - c1_) - c5_ * c9_) / (c6_ * (c6_ - c1_)) * exp6_;                  // T
  // G
  p_(2, 0) = -c2_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (c2_ * c8_ - epsilon_ * (c3_ - c1_)) / (c3_ * (c3_ - c1_)) * exp3_; // A
  p_(2, 1) = c2_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
  p_(2, 2) = -c2_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (c2_ * c8_ - epsilon_ * (c3_ - c1_)) / (c3_ * (c3_ - c1_)) * exp3_;   // G
  p_(2, 3) = c2_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (lambda_ * sigma_ - gamma_ * beta_) / (c6_ * (c6_ - c1_)) * exp6_;           // T, U
  // T, U
  p_(3, 0) = c5_ * c8_ / (c1_ * (c3_ - c1_)) * exp1_ + (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // A
  p_(3, 1) = -c5_ * c7_ / (c1_ * (c6_ - c1_)) * exp1_ + (c5_ * c7_ - sigma_ * (c6_ - c1_)) / (c6_ * (c6_ - c1_)) * exp6_;                  // C
  p_(3, 2) = c5_ * c4_ / (c1_ * (c3_ - c1_)) * exp1_ - (epsilon_ * kappa_ - delta_ * alpha_) / (c3_ * (c3_ - c1_)) * exp3_; // G
  p_(3, 3) = -c5_ * c9_ / (c1_ * (c6_ - c1_)) * exp1_ - (c5_ * c7_ - sigma_ * (c6_ - c1_)) / (c6_ * (c6_ - c1_)) * exp6_;

  return p_;
}

/******************************************************************************/
void RN95::setFreq(map<int, double>& freqs)
{
  setParameterValue("thetaR", (freqs[0] + freqs[2]) / (freqs[0] + freqs[1] + freqs[2] + freqs[3]));
  setParameterValue("thetaC", freqs[1] / (freqs[1] + freqs[3]));
  setParameterValue("thetaG", freqs[2] / (freqs[0] + freqs[2]));

  updateMatrices();
}

/******************************************************************************/

