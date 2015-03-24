//
// File: T92.cpp
// Created by:  Julien Dutheil
// Created on: Mon May 26 14:41:24 2003
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

#include "T92.h"
#include "../FrequenciesSet/NucleotideFrequenciesSet.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/******************************************************************************/

T92::T92(const NucleicAlphabet* alpha, double kappa, double theta) :
  AbstractParameterAliasable("T92."),
  AbstractReversibleNucleotideSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "T92."),
  kappa_(kappa),
  theta_(theta),
  k_(),
  r_(),
  piA_((1. - theta_) / 2.),
  piC_(theta_ / 2.),
  piG_(theta_ / 2.),
  piT_((1. - theta_) / 2.),
  exp1_(),
  exp2_(),
  l_(),
  p_(size_, size_)
{
  addParameter_(new Parameter("T92.kappa", kappa, &Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("T92.theta", theta, &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  p_.resize(size_, size_);
  updateMatrices();
}

/******************************************************************************/

void T92::updateMatrices()
{
  kappa_ = getParameterValue("kappa");
  theta_ = getParameterValue("theta");
  piA_ = (1 - theta_) / 2;
  piC_ = theta_ / 2;
  piG_ = theta_ / 2;
  piT_ = (1 - theta_) / 2;
  k_ = (kappa_ + 1.) / 2.;
  r_ = 2. / (1. + 2. * theta_ * kappa_ - 2. * theta_ * theta_ * kappa_);

  freq_[0] = piA_;
  freq_[1] = piC_;
  freq_[2] = piG_;
  freq_[3] = piT_;

  generator_(0, 0) = -(1. +        theta_ * kappa_) / 2;
  generator_(1, 1) = -(1. + (1. - theta_) * kappa_) / 2;
  generator_(2, 2) = -(1. + (1. - theta_) * kappa_) / 2;
  generator_(3, 3) = -(1. +        theta_ * kappa_) / 2;

  generator_(1, 0) = (1. - theta_) / 2;
  generator_(3, 0) = (1. - theta_) / 2;
  generator_(0, 1) = theta_ / 2;
  generator_(2, 1) = theta_ / 2;
  generator_(1, 2) = theta_ / 2;
  generator_(3, 2) = theta_ / 2;
  generator_(0, 3) = (1. - theta_) / 2;
  generator_(2, 3) = (1. - theta_) / 2;

  generator_(2, 0) = kappa_ * (1. - theta_) / 2;
  generator_(3, 1) = kappa_ * theta_ / 2;
  generator_(0, 2) = kappa_ * theta_ / 2;
  generator_(1, 3) = kappa_ * (1. - theta_) / 2;

  // Normalization:
  MatrixTools::scale(generator_, r_);

  // Exchangeability:
  exchangeability_(0, 0) = generator_(0, 0) * 2. / (1. - theta_);
  exchangeability_(0, 1) = generator_(0, 1) * 2. / theta_;
  exchangeability_(0, 2) = generator_(0, 2) * 2. / theta_;
  exchangeability_(0, 3) = generator_(0, 3) * 2. / (1. - theta_);

  exchangeability_(1, 0) = generator_(1, 0) * 2. / (1. - theta_);
  exchangeability_(1, 1) = generator_(1, 1) * 2 / theta_;
  exchangeability_(1, 2) = generator_(1, 2) * 2 / theta_;
  exchangeability_(1, 3) = generator_(1, 3) * 2 / (1. - theta_);

  exchangeability_(2, 0) = generator_(2, 0) * 2. / (1. - theta_);
  exchangeability_(2, 1) = generator_(2, 1) * 2 / theta_;
  exchangeability_(2, 2) = generator_(2, 2) * 2 / theta_;
  exchangeability_(2, 3) = generator_(2, 3) * 2 / (1. - theta_);

  exchangeability_(3, 0) = generator_(3, 0) * 2. / (1. - theta_);
  exchangeability_(3, 1) = generator_(3, 1) * 2. / theta_;
  exchangeability_(3, 2) = generator_(3, 2) * 2. / theta_;
  exchangeability_(3, 3) = generator_(3, 3) * 2. / (1. - theta_);

  // Eigen values:
  eigenValues_[0] = 0;
  eigenValues_[1] = eigenValues_[2] = -r_ * (1. + kappa_) / 2;
  eigenValues_[3] = -r_;

  // Eigen vectors:
  leftEigenVectors_(0, 0) = -(theta_ - 1.) / 2.;
  leftEigenVectors_(0, 1) = theta_ / 2.;
  leftEigenVectors_(0, 2) = theta_ / 2.;
  leftEigenVectors_(0, 3) = -(theta_ - 1.) / 2.;

  leftEigenVectors_(1, 0) = 0.;
  leftEigenVectors_(1, 1) = -(theta_ - 1.);
  leftEigenVectors_(1, 2) = 0.;
  leftEigenVectors_(1, 3) = theta_ - 1.;

  leftEigenVectors_(2, 0) = theta_;
  leftEigenVectors_(2, 1) = 0.;
  leftEigenVectors_(2, 2) = -theta_;
  leftEigenVectors_(2, 3) = 0.;

  leftEigenVectors_(3, 0) = -(theta_ - 1.) / 2.;
  leftEigenVectors_(3, 1) = -theta_ / 2.;
  leftEigenVectors_(3, 2) = theta_ / 2.;
  leftEigenVectors_(3, 3) = (theta_ - 1.) / 2.;


  rightEigenVectors_(0, 0) = 1.;
  rightEigenVectors_(0, 1) = 0.;
  rightEigenVectors_(0, 2) = 1.;
  rightEigenVectors_(0, 3) = 1.;

  rightEigenVectors_(1, 0) = 1.;
  rightEigenVectors_(1, 1) = 1.;
  rightEigenVectors_(1, 2) = 0.;
  rightEigenVectors_(1, 3) = -1.;

  rightEigenVectors_(2, 0) = 1.;
  rightEigenVectors_(2, 1) = 0.;
  rightEigenVectors_(2, 2) = (theta_ - 1.) / theta_;
  rightEigenVectors_(2, 3) = 1.;

  rightEigenVectors_(3, 0) = 1.;
  rightEigenVectors_(3, 1) = theta_ / (theta_ - 1.);
  rightEigenVectors_(3, 2) = 0;
  rightEigenVectors_(3, 3) = -1.;
}

/******************************************************************************/

double T92::Pij_t(size_t i, size_t j, double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-l_);
  exp2_ = exp(-k_ * l_);

  switch (i)
  {
  // A
  case 0: {
    switch (j)
    {
    case 0: return piA_ * (1. + exp1_) + theta_ * exp2_; // A
    case 1: return piC_ * (1. - exp1_);                 // C
    case 2: return piG_ * (1. + exp1_) - theta_ * exp2_; // G
    case 3: return piT_ * (1. - exp1_);                 // T, U
    }
  }
  // C
  case 1: {
    switch (j)
    {
    case 0: return piA_ * (1. - exp1_);                        // A
    case 1: return piC_ * (1. + exp1_) + (1. - theta_) * exp2_; // C
    case 2: return piG_ * (1. - exp1_);                        // G
    case 3: return piT_ * (1. + exp1_) - (1. - theta_) * exp2_; // T, U
    }
  }
  // G
  case 2: {
    switch (j)
    {
    case 0: return piA_ * (1. + exp1_) - (1. - theta_) * exp2_; // A
    case 1: return piC_ * (1. - exp1_);                        // C
    case 2: return piG_ * (1. + exp1_) + (1. - theta_) * exp2_; // G
    case 3: return piT_ * (1. - exp1_);                        // T, U
    }
  }
  // T, U
  case 3: {
    switch (j)
    {
    case 0: return piA_ * (1. - exp1_);                 // A
    case 1: return piC_ * (1. + exp1_) - theta_ * exp2_; // C
    case 2: return piG_ * (1. - exp1_);                 // G
    case 3: return piT_ * (1. + exp1_) + theta_ * exp2_; // T, U
    }
  }
  }
  return 0;
}

/******************************************************************************/

double T92::dPij_dt(size_t i, size_t j, double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-l_);
  exp2_ = exp(-k_ * l_);

  switch (i)
  {
  // A
  case 0: {
    switch (j)
    {
    case 0: return rate_ * r_ * (piA_ * -exp1_ + theta_ * -k_ * exp2_); // A
    case 1: return rate_ * r_ * (piC_ *   exp1_);                       // C
    case 2: return rate_ * r_ * (piG_ * -exp1_ - theta_ * -k_ * exp2_); // G
    case 3: return rate_ * r_ * (piT_ *   exp1_);                       // T, U
    }
  }
  // C
  case 1: {
    switch (j)
    {
    case 0: return rate_ * r_ * (piA_ *   exp1_);                              // A
    case 1: return rate_ * r_ * (piC_ * -exp1_ + (1. - theta_) * -k_ * exp2_); // C
    case 2: return rate_ * r_ * (piG_ *   exp1_);                              // G
    case 3: return rate_ * r_ * (piT_ * -exp1_ - (1. - theta_) * -k_ * exp2_); // T, U
    }
  }
  // G
  case 2: {
    switch (j)
    {
    case 0: return rate_ * r_ * (piA_ * -exp1_ - (1. - theta_) * -k_ * exp2_); // A
    case 1: return rate_ * r_ * (piC_ *   exp1_);                              // C
    case 2: return rate_ * r_ * (piG_ * -exp1_ + (1. - theta_) * -k_ * exp2_); // G
    case 3: return rate_ * r_ * (piT_ *   exp1_);                              // T, U
    }
  }
  // T, U
  case 3: {
    switch (j)
    {
    case 0: return rate_ * r_ * (piA_ *   exp1_);                       // A
    case 1: return rate_ * r_ * (piC_ * -exp1_ - theta_ * -k_ * exp2_); // C
    case 2: return rate_ * r_ * (piG_ *   exp1_);                       // G
    case 3: return rate_ * r_ * (piT_ * -exp1_ + theta_ * -k_ * exp2_); // T, U
    }
  }
  }
  return 0;
}

/******************************************************************************/

double T92::d2Pij_dt2(size_t i, size_t j, double d) const
{
  double k2_ = k_ * k_;
  l_ = rate_ * r_ * d;
  double r2 = rate_ * rate_ * r_ * r_;
  exp1_ = exp(-l_);
  exp2_ = exp(-k_ * l_);

  switch (i)
  {
  // A
  case 0: {
    switch (j)
    {
    case 0: return r2 * (piA_ *   exp1_ + theta_ * k2_ * exp2_); // A
    case 1: return r2 * (piC_ * -exp1_);                       // C
    case 2: return r2 * (piG_ *   exp1_ - theta_ * k2_ * exp2_); // G
    case 3: return r2 * (piT_ * -exp1_);                       // T, U
    }
  }
  // C
  case 1: {
    switch (j)
    {
    case 0: return r2 * (piA_ * -exp1_);                              // A
    case 1: return r2 * (piC_ *   exp1_ + (1. - theta_) * k2_ * exp2_); // C
    case 2: return r2 * (piG_ * -exp1_);                              // G
    case 3: return r2 * (piT_ *   exp1_ - (1. - theta_) * k2_ * exp2_); // T, U
    }
  }
  // G
  case 2: {
    switch (j)
    {
    case 0: return r2 * (piA_ *   exp1_ - (1. - theta_) * k2_ * exp2_); // A
    case 1: return r2 * (piC_ * -exp1_);                              // C
    case 2: return r2 * (piG_ *   exp1_ + (1. - theta_) * k2_ * exp2_); // G
    case 3: return r2 * (piT_ * -exp1_);                              // T, U
    }
  }
  // T, U
  case 3: {
    switch (j)
    {
    case 0: return r2 * (piA_ * -exp1_);                       // A
    case 1: return r2 * (piC_ *   exp1_ - theta_ * k2_ * exp2_); // C
    case 2: return r2 * (piG_ * -exp1_);                       // G
    case 3: return r2 * (piT_ *   exp1_ + theta_ * k2_ * exp2_); // T, U
    }
  }
  }
  return 0;
}

/******************************************************************************/

const Matrix<double>& T92::getPij_t(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-l_);
  exp2_ = exp(-k_ * l_);

  // A
  p_(0, 0) = piA_ * (1. + exp1_) + theta_ * exp2_; // A
  p_(0, 1) = piC_ * (1. - exp1_);                  // C
  p_(0, 2) = piG_ * (1. + exp1_) - theta_ * exp2_; // G
  p_(0, 3) = piT_ * (1. - exp1_);                  // T, U

  // C
  p_(1, 0) = piA_ * (1. - exp1_);                         // A
  p_(1, 1) = piC_ * (1. + exp1_) + (1. - theta_) * exp2_; // C
  p_(1, 2) = piG_ * (1. - exp1_);                         // G
  p_(1, 3) = piT_ * (1. + exp1_) - (1. - theta_) * exp2_; // T, U

  // G
  p_(2, 0) = piA_ * (1. + exp1_) - (1. - theta_) * exp2_; // A
  p_(2, 1) = piC_ * (1. - exp1_);                         // C
  p_(2, 2) = piG_ * (1. + exp1_) + (1. - theta_) * exp2_; // G
  p_(2, 3) = piT_ * (1. - exp1_);                         // T, U

  // T, U
  p_(3, 0) = piA_ * (1. - exp1_);                  // A
  p_(3, 1) = piC_ * (1. + exp1_) - theta_ * exp2_; // C
  p_(3, 2) = piG_ * (1. - exp1_);                  // G
  p_(3, 3) = piT_ * (1. + exp1_) + theta_ * exp2_; // T, U

  return p_;
}

const Matrix<double>& T92::getdPij_dt(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-l_);
  exp2_ = exp(-k_ * l_);

  // A
  p_(0, 0) = rate_ * r_ * (piA_ * -exp1_ + theta_ * -k_ * exp2_); // A
  p_(0, 1) = rate_ * r_ * (piC_ *   exp1_);                        // C
  p_(0, 2) = rate_ * r_ * (piG_ * -exp1_ - theta_ * -k_ * exp2_); // G
  p_(0, 3) = rate_ * r_ * (piT_ *   exp1_);                        // T, U

  // C
  p_(1, 0) = rate_ * r_ * (piA_ *   exp1_);                               // A
  p_(1, 1) = rate_ * r_ * (piC_ * -exp1_ + (1. - theta_) * -k_ * exp2_); // C
  p_(1, 2) = rate_ * r_ * (piG_ *   exp1_);                               // G
  p_(1, 3) = rate_ * r_ * (piT_ * -exp1_ - (1. - theta_) * -k_ * exp2_); // T, U

  // G
  p_(2, 0) = rate_ * r_ * (piA_ * -exp1_ - (1. - theta_) * -k_ * exp2_); // A
  p_(2, 1) = rate_ * r_ * (piC_ *   exp1_);                               // C
  p_(2, 2) = rate_ * r_ * (piG_ * -exp1_ + (1. - theta_) * -k_ * exp2_); // G
  p_(2, 3) = rate_ * r_ * (piT_ *   exp1_);                               // T, U

  // T, U
  p_(3, 0) = rate_ * r_ * (piA_ *   exp1_);                        // A
  p_(3, 1) = rate_ * r_ * (piC_ * -exp1_ - theta_ * -k_ * exp2_); // C
  p_(3, 2) = rate_ * r_ * (piG_ *   exp1_);                        // G
  p_(3, 3) = rate_ * r_ * (piT_ * -exp1_ + theta_ * -k_ * exp2_); // T, U

  return p_;
}

const Matrix<double>& T92::getd2Pij_dt2(double d) const
{
  double k2 = k_ * k_;
  l_ = rate_ * r_ * d;
  double r2 = rate_ * rate_ * r_ * r_;
  exp1_ = exp(-l_);
  exp2_ = exp(-k_ * l_);

  // A
  p_(0, 0) = r2 * (piA_ *   exp1_ + theta_ * k2 * exp2_); // A
  p_(0, 1) = r2 * (piC_ * -exp1_);                      // C
  p_(0, 2) = r2 * (piG_ *   exp1_ - theta_ * k2 * exp2_); // G
  p_(0, 3) = r2 * (piT_ * -exp1_);                      // T, U

  // C
  p_(1, 0) = r2 * (piA_ * -exp1_);                             // A
  p_(1, 1) = r2 * (piC_ *   exp1_ + (1. - theta_) * k2 * exp2_); // C
  p_(1, 2) = r2 * (piG_ * -exp1_);                             // G
  p_(1, 3) = r2 * (piT_ *   exp1_ - (1. - theta_) * k2 * exp2_); // T, U

  // G
  p_(2, 0) = r2 * (piA_ *   exp1_ - (1. - theta_) * k2 * exp2_); // A
  p_(2, 1) = r2 * (piC_ * -exp1_);                             // C
  p_(2, 2) = r2 * (piG_ *   exp1_ + (1. - theta_) * k2 * exp2_); // G
  p_(2, 3) = r2 * (piT_ * -exp1_);                             // T, U

  // T, U
  p_(3, 0) = r2 * (piA_ * -exp1_);                      // A
  p_(3, 1) = r2 * (piC_ *   exp1_ - theta_ * k2 * exp2_); // C
  p_(3, 2) = r2 * (piG_ * -exp1_);                      // G
  p_(3, 3) = r2 * (piT_ *   exp1_ + theta_ * k2 * exp2_); // T, U

  return p_;
}

/******************************************************************************/

void T92::setFreq(std::map<int, double>& freqs)
{
  double f = (freqs[1] + freqs[2]) / (freqs[0] + freqs[1] + freqs[2] + freqs[3]);
  setParameterValue("theta", f);
  updateMatrices();
}

/******************************************************************************/

