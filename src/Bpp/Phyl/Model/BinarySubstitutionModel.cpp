//
// File: BinarySubstitutionModel.cpp
// Created by: Laurent Gueguen
// Created on: 2009
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

#include "BinarySubstitutionModel.h"

// From the STL:
#include <cmath>

using namespace bpp;

#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace std;

/******************************************************************************/

BinarySubstitutionModel::BinarySubstitutionModel(const BinaryAlphabet* alpha, double kappa) :
  AbstractParameterAliasable("Binary."),
  //AbstractSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "Binary."),
  AbstractReversibleSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "Binary."),
  kappa_(kappa),
  lambda_(0),
  exp_(0),
  p_(size_,size_)
{
  addParameter_(new Parameter(getNamespace() + "kappa", kappa_, &Parameter::R_PLUS_STAR));
  updateMatrices();
}

/******************************************************************************/

void BinarySubstitutionModel::updateMatrices()
{
  kappa_ = getParameterValue("kappa"); // alpha/beta
  lambda_ = (kappa_ + 1) * (kappa_ + 1) / (2 * kappa_);

  // Frequences:
  freq_[0] = 1 / (kappa_ + 1);
  freq_[1] = kappa_ / (kappa_ + 1);

  // Generator:
  generator_(0, 0) = rate_ * -(kappa_ + 1) / 2;
  generator_(0, 1) = rate_ * (kappa_ + 1) / 2;
  generator_(1, 0) = rate_ * (kappa_ + 1) / (2 * kappa_);
  generator_(1, 1) = rate_ * -(kappa_ + 1) / (2 * kappa_);

  // Eigen values:
  eigenValues_[0] = 0;
  eigenValues_[1] = -rate_ * lambda_;

  // Eigen vectors:
  leftEigenVectors_(0,0) = 1 / (kappa_ + 1);
  leftEigenVectors_(0,1) = kappa_ / (kappa_ + 1);
  if (kappa_ != 1.0)
  {
    leftEigenVectors_(1,0) = (kappa_ - 1) / (kappa_ + 1);
    leftEigenVectors_(1,1) = -(kappa_ - 1) / (kappa_ + 1);
  }
  else
  {
    leftEigenVectors_(1,0) = 1;
    leftEigenVectors_(1,1) = -1;
  }

  rightEigenVectors_(0,0) = 1;
  rightEigenVectors_(1,0) = 1;

  if (kappa_ != 1.0)
  {
    rightEigenVectors_(0,1) = kappa_ / (kappa_ - 1);
    rightEigenVectors_(1,1) = -1 / (kappa_ - 1);
  }
  else
  {
    rightEigenVectors_(0,1) = 1 / 2;
    rightEigenVectors_(1,1) = -1 / 2;
  }
}

/******************************************************************************/

double BinarySubstitutionModel::Pij_t(size_t i, size_t j, double d) const
{
  exp_ = exp(-lambda_ * rate_ * d);

  switch (i)
  {
  case 0: {
    switch (j)
    {
    case 0: return (1 + kappa_ * exp_) / (kappa_ + 1);
    case 1: return kappa_ / (kappa_ + 1) * (1 - exp_);
    }
  }
  case 1: {
    switch (j)
    {
    case 0: return (1 - exp_) / (kappa_ + 1);
    case 1: return (kappa_ + exp_) / (kappa_ + 1);
    }
  }
  }
  return 0;
}

/******************************************************************************/

double BinarySubstitutionModel::dPij_dt(size_t i, size_t j, double d) const
{
  exp_ = rate_ * exp(-lambda_ * rate_ * d);

  switch (i)
  {
  case 0: {
    switch (j)
    {
    case 0: return -(kappa_ + 1) / 2 * exp_;
    case 1: return (kappa_ + 1) / 2 * exp_;
    }
  }
  case 1: {
    switch (j)
    {
    case 0: return (kappa_ + 1) / (2 * kappa_) * exp_;
    case 1: return -(kappa_ + 1) / (2 * kappa_) * exp_;
    }
  }
  }
  return 0;
}

/******************************************************************************/

double BinarySubstitutionModel::d2Pij_dt2(size_t i, size_t j, double d) const
{
  exp_ = rate_ * rate_ * exp(-lambda_ * rate_ * d);

  switch (i)
  {
  case 0: {
    switch (j)
    {
    case 0: return lambda_ * (kappa_ + 1) / 2 * exp_;
    case 1: return -lambda_ * (kappa_ + 1) / 2 * exp_;
    }
  }
  case 1: {
    switch (j)
    {
    case 0: return -lambda_ * (kappa_ + 1) / (2 * kappa_) * exp_;
    case 1: return lambda_ * (kappa_ + 1) / (2 * kappa_) * exp_;
    }
  }
  }
  return 0;
}

/******************************************************************************/

const Matrix<double>& BinarySubstitutionModel::getPij_t(double d) const
{
  exp_ = exp(-lambda_ * rate_ * d);

  p_(0,0) = (1 + kappa_ * exp_) / (kappa_ + 1);
  p_(0,1) = kappa_ / (kappa_ + 1) * (1 - exp_);

  p_(1,0) =  (1 - exp_) / (kappa_ + 1);
  p_(1,1) = (kappa_ + exp_) / (kappa_ + 1);

  return p_;
}

const Matrix<double>& BinarySubstitutionModel::getdPij_dt(double d) const
{
  exp_ = rate_ * exp(-lambda_ * rate_ * d);

  p_(0,0) = -(kappa_ + 1) / 2 * exp_;
  p_(0,1) = (kappa_ + 1) / 2 * exp_;

  p_(1,0) = (kappa_ + 1) / (2 * kappa_) * exp_;
  p_(1,1) = -(kappa_ + 1) / (2 * kappa_) * exp_;

  return p_;
}

const Matrix<double>& BinarySubstitutionModel::getd2Pij_dt2(double d) const
{
  exp_ = rate_ * rate_ * exp(-lambda_ * rate_ * d);

  p_(0,0) = lambda_ * (kappa_ + 1) / 2 * exp_;
  p_(0,1) = -lambda_ * (kappa_ + 1) / 2 * exp_;
  p_(1,0) = -lambda_ * (kappa_ + 1) / (2 * kappa_) * exp_;
  p_(1,1) = lambda_ * (kappa_ + 1) / (2 * kappa_) * exp_;

  return p_;
}

/******************************************************************************/

void BinarySubstitutionModel::setFreq(std::map<int, double>& freqs)
{
  kappa_ = freqs[1] / freqs[0];
  setParameterValue("kappa",kappa_);
  updateMatrices();
}
