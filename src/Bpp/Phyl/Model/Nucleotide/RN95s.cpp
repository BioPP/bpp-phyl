//
// File: RN95s.cpp
// Created by: Laurent Guéguen
// Created on: jeudi 24 février 2011, à 20h 42
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "RN95s.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

// From the STL:
#include <cmath>

using namespace bpp;
using namespace std;

/******************************************************************************/

RN95s::RN95s(const NucleicAlphabet* alphabet,
             double alpha,
             double beta,
             double gamma,
             double delta) :
  AbstractParameterAliasable("RN95s."),
  AbstractNucleotideSubstitutionModel(alphabet, new CanonicalStateMap(alphabet, false), "RN95s."),
  alpha_(),
  beta_(),
  gamma_(),
  delta_(),
  r_(),
  c3_(),
  c4_(),
  c8_(),
  p_(size_, size_),
  exp1_(),
  exp3_(),
  l_()
{
  double f = 2 * (gamma + delta);

  alpha_ = alpha / f;
  beta_ = beta / f;
  gamma_ = gamma / f;
  delta_ = delta / f;

  double piA = 0.5 * (beta_ + delta_) / (alpha_ + beta_ + 0.5);
  double alphaP = (2 * alpha_ * piA + ((0.5 - piA < gamma_) ? 0.5 - piA : gamma_)) / (0.5 - piA);
  addParameter_(new Parameter("RN95s.thetaA", piA, new IntervalConstraint(0, 0.5, false, false), true));
  addParameter_(new Parameter("RN95s.gamma", gamma_, new IntervalConstraint(0, 0.5, false, false), true));
  addParameter_(new Parameter("RN95s.alphaP", alphaP, new IntervalConstraint(1, 1, false), true));

  updateMatrices();
}

/******************************************************************************/
void RN95s::updateMatrices()
{
  freq_[0]  = getParameterValue("thetaA");
  double alphaP  = getParameterValue("alphaP");
  gamma_  = getParameterValue("gamma");
  alpha_  = (alphaP * (0.5 - freq_[0]) - ((0.5 - freq_[0] < gamma_) ? 0.5 - freq_[0] : gamma_)) / (2 * freq_[0]);
  delta_  = 0.5 - gamma_;
  beta_   = (2 * freq_[0] * (alpha_ + 0.5) - delta_) / (1 - 2 * freq_[0]);

  // stationnary frequencies

  freq_[1] = 0.5 - freq_[0];
  freq_[2] = freq_[1];
  freq_[3] = freq_[0];

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
  {
    x += generator_(i, i) * freq_[i];
  }

  r_ = -1 / x;

  MatrixTools::scale(generator_, r_);

  // variables for calculation purposes

  c3_ = alpha_ + gamma_ + beta_ + delta_;
  c4_ = gamma_ - alpha_;
  c8_ = delta_ - beta_;


  // eigen vectors and values

  eigenValues_[0] = 0;
  eigenValues_[1] = -r_;
  eigenValues_[2] = -c3_ * r_;
  eigenValues_[3] = -c3_ * r_;

  rightEigenVectors_(0, 0) = 1.;
  rightEigenVectors_(1, 0) = 1.;
  rightEigenVectors_(2, 0) = 1.;
  rightEigenVectors_(3, 0) = 1.;

  rightEigenVectors_(0, 1) = 1.;
  rightEigenVectors_(1, 1) = -1;
  rightEigenVectors_(2, 1) = 1.;
  rightEigenVectors_(3, 1) = -1;

  rightEigenVectors_(0, 2) = (alpha_ * (0.5 - c3_) + gamma_ * 0.5) / (delta_ * (c3_ - 0.5) - beta_ * 0.5);
  rightEigenVectors_(1, 2) = 1.;
  rightEigenVectors_(2, 2) = (-beta_ * (0.5 - c3_) - delta_ * 0.5) / (delta_ * (c3_ - 0.5) - beta_ * 0.5);
  rightEigenVectors_(3, 2) = 1.;

  rightEigenVectors_(0, 3) = 1.;
  rightEigenVectors_(1, 3) = (beta_ * (0.5 - c3_) + delta_ * 0.5) / (gamma_ * (c3_ - 0.5) - alpha_ * 0.5);
  rightEigenVectors_(2, 3) = 1;
  rightEigenVectors_(3, 3) = (-alpha_ * (0.5 - c3_) - gamma_ * 0.5) / (gamma_ * (c3_ - 0.5) - alpha_ * 0.5);

// Need formula

  if (enableEigenDecomposition())
  {
    try
    {
      MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);
      isNonSingular_ = true;
      isDiagonalizable_ = true;
      for (unsigned int i = 0; i < size_ && isDiagonalizable_; i++)
      {
        if (abs(iEigenValues_[i]) > NumConstants::TINY())
          isDiagonalizable_ = false;
      }
    }
    catch (ZeroDivisionException& e)
    {
      ApplicationTools::displayMessage("Singularity during  diagonalization. Taylor series used instead.");
      
      isNonSingular_ = false;
      isDiagonalizable_ = false;
      MatrixTools::Taylor(generator_, 30, vPowGen_);
    }
    
    // and the exchangeability_
    for (size_t i = 0; i < size_; i++)
      for (size_t j = 0; j < size_; j++)
        exchangeability_(i,j) = generator_(i,j) / freq_[j];
  }
}

/******************************************************************************/
double RN95s::Pij_t(size_t i, size_t j, double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-l_);
  exp3_ = exp(-c3_ * l_);

  switch (i)
  {
    {
    // A
    case 0: {
      switch (j)
      {
      case 0: return freq_[0] - 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (alpha_ * (c3_ - 1.) - 0.5 * c4_) / (c3_ * (c3_ - 1.)) * exp3_; // A
      case 1: return freq_[1] + 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
      case 2: return freq_[2] - 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (alpha_ * (c3_ - 1.) - 0.5 * c4_) / (c3_ * (c3_ - 1.)) * exp3_;  // G
      case 3: return freq_[3] + 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;           // T, U
      }
    }
    // C
    case 1: {
      switch (j)
      {
      case 0: return freq_[0] + 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // A
      case 1: return freq_[1] - 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * (c3_ - 1.) - 0.5 * c8_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
      case 2: return freq_[2] + 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // G
      case 3: return freq_[3] - 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * (c3_ - 1.) - 0.5 * c8_) / (c3_ * (c3_ - 1.)) * exp3_;                  // T
      }
    }
    // G
    case 2: {
      switch (j)
      {
      case 0: return freq_[0] - 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (0.5 * c8_ - beta_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_; // A
      case 1: return freq_[1] + 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
      case 2: return freq_[2] - 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (0.5 * c8_ - beta_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;   // G
      case 3: return freq_[3] + 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;           // T, U
      }
    }
    }
  // T, U
  case 3: {
    switch (j)
    {
    case 0: return freq_[0] + 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // A
    case 1: return freq_[1] - 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (0.5 * c4_ - alpha_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
    case 2: return freq_[2] + 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // G
    case 3: return freq_[3] - 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (0.5 * c4_ - alpha_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;                  // T
    }
  }
  }
  return 0;
}

/******************************************************************************/
double RN95s::dPij_dt(size_t i, size_t j, double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = -1.* rate_* r_* exp(-1. * l_);
  exp3_ = -c3_* rate_* r_* exp(-c3_ * l_);

  switch (i)
  {
    {
    // A
    case 0: {
      switch (j)
      {
      case 0: return -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (alpha_ * (c3_ - 1.) - 0.5 * c4_) / (c3_ * (c3_ - 1.)) * exp3_; // A
      case 1: return 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
      case 2: return -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (alpha_ * (c3_ - 1.) - 0.5 * c4_) / (c3_ * (c3_ - 1.)) * exp3_;  // G
      case 3: return 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;           // T, U
      }
    }
    // C
    case 1: {
      switch (j)
      {
      case 0: return 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // A
      case 1: return -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * (c3_ - 1.) - 0.5 * c8_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
      case 2: return 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // G
      case 3: return -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * (c3_ - 1.) - 0.5 * c8_) / (c3_ * (c3_ - 1.)) * exp3_;                  // T
      }
    }
    // G
    case 2: {
      switch (j)
      {
      case 0: return -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (0.5 * c8_ - beta_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_; // A
      case 1: return 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
      case 2: return -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (0.5 * c8_ - beta_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;   // G
      case 3: return 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;           // T, U
      }
    }
    }
  // T, U
  case 3: {
    switch (j)
    {
    case 0: return 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // A
    case 1: return -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (0.5 * c4_ - alpha_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
    case 2: return 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // G
    case 3: return -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (0.5 * c4_ - alpha_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;                  // T
    }
  }
  }
  return 0;
}

/******************************************************************************/
double RN95s::d2Pij_dt2(size_t i, size_t j, double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = 1. * rate_ * r_ * 1.* rate_* r_* exp(-1. * l_);
  exp3_ = c3_ * rate_ * r_ * c3_ * rate_ * r_ * exp(-c3_ * l_);

  switch (i)
  {
    {
    // A
    case 0: {
      switch (j)
      {
      case 0: return -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (alpha_ * (c3_ - 1.) - 0.5 * c4_) / (c3_ * (c3_ - 1.)) * exp3_; // A
      case 1: return 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
      case 2: return -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (alpha_ * (c3_ - 1.) - 0.5 * c4_) / (c3_ * (c3_ - 1.)) * exp3_;  // G
      case 3: return 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;           // T, U
      }
    }
    // C
    case 1: {
      switch (j)
      {
      case 0: return 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // A
      case 1: return -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * (c3_ - 1.) - 0.5 * c8_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
      case 2: return 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // G
      case 3: return -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * (c3_ - 1.) - 0.5 * c8_) / (c3_ * (c3_ - 1.)) * exp3_;                  // T
      }
    }
    // G
    case 2: {
      switch (j)
      {
      case 0: return -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (0.5 * c8_ - beta_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_; // A
      case 1: return 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
      case 2: return -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (0.5 * c8_ - beta_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;   // G
      case 3: return 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;           // T, U
      }
    }
    }
  // T, U
  case 3: {
    switch (j)
    {
    case 0: return 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // A
    case 1: return -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (0.5 * c4_ - alpha_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
    case 2: return 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // G
    case 3: return -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (0.5 * c4_ - alpha_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;                  // T
    }
  }
  }
  return 0;
}

/******************************************************************************/

const Matrix<double>& RN95s::getPij_t(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-1. * l_);
  exp3_ = exp(-c3_ * l_);

  // A
  p_(0, 0) = freq_[0] - 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (alpha_ * (c3_ - 1.) - 0.5 * c4_) / (c3_ * (c3_ - 1.)) * exp3_; // A
  p_(0, 1) = freq_[1] + 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
  p_(0, 2) = freq_[2] - 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (alpha_ * (c3_ - 1.) - 0.5 * c4_) / (c3_ * (c3_ - 1.)) * exp3_;  // G
  p_(0, 3) = freq_[3] + 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;           // T, U
  // C
  p_(1, 0) = freq_[0] + 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // A
  p_(1, 1) = freq_[1] - 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * (c3_ - 1.) - 0.5 * c8_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
  p_(1, 2) = freq_[2] + 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // G
  p_(1, 3) = freq_[3] - 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * (c3_ - 1.) - 0.5 * c8_) / (c3_ * (c3_ - 1.)) * exp3_;                  // T
  // G
  p_(2, 0) = freq_[0] - 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (0.5 * c8_ - beta_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_; // A
  p_(2, 1) = freq_[1] + 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
  p_(2, 2) = freq_[2] - 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (0.5 * c8_ - beta_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;   // G
  p_(2, 3) = freq_[3] + 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;           // T, U
  // T, U
  p_(3, 0) = freq_[0] + 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // A
  p_(3, 1) = freq_[1] - 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (0.5 * c4_ - alpha_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
  p_(3, 2) = freq_[2] + 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // G
  p_(3, 3) = freq_[3] - 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (0.5 * c4_ - alpha_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;                  // T

  return p_;
}

/******************************************************************************/

const Matrix<double>&  RN95s::getdPij_dt(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = -1.* rate_* r_* exp(-1. * l_);
  exp3_ = -c3_* rate_* r_* exp(-c3_ * l_);

  // A
  p_(0, 0) = -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (alpha_ * (c3_ - 1.) - 0.5 * c4_) / (c3_ * (c3_ - 1.)) * exp3_; // A
  p_(0, 1) =  0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
  p_(0, 2) =  -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (alpha_ * (c3_ - 1.) - 0.5 * c4_) / (c3_ * (c3_ - 1.)) * exp3_;  // G
  p_(0, 3) = 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;           // T, U
  // C
  p_(1, 0) = 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // A
  p_(1, 1) = -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * (c3_ - 1.) - 0.5 * c8_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
  p_(1, 2) = 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // G
  p_(1, 3) = -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * (c3_ - 1.) - 0.5 * c8_) / (c3_ * (c3_ - 1.)) * exp3_;                  // T
  // G
  p_(2, 0) = -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (0.5 * c8_ - beta_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_; // A
  p_(2, 1) = 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
  p_(2, 2) = -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (0.5 * c8_ - beta_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;   // G
  p_(2, 3) = 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;           // T, U
  // T, U
  p_(3, 0) = 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // A
  p_(3, 1) = -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (0.5 * c4_ - alpha_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
  p_(3, 2) = 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // G
  p_(3, 3) = -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (0.5 * c4_ - alpha_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;                  // T
  return p_;
}

/******************************************************************************/

const Matrix<double>&  RN95s::getd2Pij_dt2(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = 1. * rate_ * r_ * 1.* rate_* r_* exp(-1. * l_);
  exp3_ = c3_ * rate_ * r_ * c3_ * rate_ * r_ * exp(-c3_ * l_);

  // A
  p_(0, 0) = -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (alpha_ * (c3_ - 1.) - 0.5 * c4_) / (c3_ * (c3_ - 1.)) * exp3_; // A
  p_(0, 1) =  0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
  p_(0, 2) =  -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (alpha_ * (c3_ - 1.) - 0.5 * c4_) / (c3_ * (c3_ - 1.)) * exp3_;  // G
  p_(0, 3) = 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;           // T, U
  // C
  p_(1, 0) = 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // A
  p_(1, 1) = -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * (c3_ - 1.) - 0.5 * c8_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
  p_(1, 2) = 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // G
  p_(1, 3) = -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * (c3_ - 1.) - 0.5 * c8_) / (c3_ * (c3_ - 1.)) * exp3_;                  // T
  // G
  p_(2, 0) = -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (0.5 * c8_ - beta_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_; // A
  p_(2, 1) = 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
  p_(2, 2) = -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (0.5 * c8_ - beta_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;   // G
  p_(2, 3) = 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (delta_ * alpha_ - gamma_ * beta_) / (c3_ * (c3_ - 1.)) * exp3_;           // T, U
  // T, U
  p_(3, 0) = 0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ + (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // A
  p_(3, 1) = -0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ + (0.5 * c4_ - alpha_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;                  // C
  p_(3, 2) = 0.5 * c4_ / (1. * (c3_ - 1.)) * exp1_ - (beta_ * gamma_ - delta_ * alpha_) / (c3_ * (c3_ - 1.)) * exp3_; // G
  p_(3, 3) = -0.5 * c8_ / (1. * (c3_ - 1.)) * exp1_ - (0.5 * c4_ - alpha_ * (c3_ - 1.)) / (c3_ * (c3_ - 1.)) * exp3_;

  return p_;
}

/******************************************************************************/
void RN95s::setFreq(map<int, double>& freqs)
{
  setParameterValue("thetaA", (freqs[0] + freqs[3]) / 2);

  updateMatrices();
}

/******************************************************************************/

