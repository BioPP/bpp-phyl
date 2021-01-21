//
// File: TwoParameterBinarySubstitutionModel.cpp
// Created by: Laurent Gueguen
// Created on: 2009
//

/*
   Copyright or � or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "TwoParameterBinarySubstitutionModel.h"

// From the STL:
#include <cmath>

using namespace bpp;

#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace std;

/******************************************************************************/

TwoParameterBinarySubstitutionModel::TwoParameterBinarySubstitutionModel(const BinaryAlphabet* alpha, double mu, double pi0):
  AbstractParameterAliasable("TwoParameterBinary."),
  //AbstractReversibleSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "TwoParameterBinary."),
  AbstractReversibleSubstitutionModel(alpha, std::shared_ptr<const StateMap>(new CanonicalStateMap(alpha, false)), "TwoParameterBinary."),
  mu_(mu),
  pi0_(pi0),
  lambda_(0),
  exp_(0),
  p_(size_,size_)
{
  addParameter_(new Parameter(getNamespace() + "mu", mu_, std::make_shared<IntervalConstraint>(NumConstants::MILLI(), 100, false, false)));
  addParameter_(new Parameter(getNamespace() + "pi0", pi0_, std::make_shared<IntervalConstraint>(0.05, 0.95, true, true)));
  updateMatrices();
}

/******************************************************************************/

void TwoParameterBinarySubstitutionModel::updateMatrices()
{
  mu_ = getParameterValue("mu");
  rate_ = mu_;
  pi0_ = getParameterValue("pi0");
  lambda_ = 1;

  // Frequences:
  freq_[0] = pi0_;
  freq_[1] = 1-pi0_;

  // Generator:
  generator_(0, 0) = -1 * rate_ * freq_[1];
  generator_(0, 1) = rate_ * freq_[1];
  generator_(1, 0) = rate_ * freq_[0];
  generator_(1, 1) = -1 * rate_ * freq_[0];

  // Eigen values:
  eigenValues_[0] = 0;
  eigenValues_[1] = -1 * mu_;

  // Eigen vectors:
  leftEigenVectors_(0,0) = pi0_;
  leftEigenVectors_(0,1) = 1-pi0_;
  leftEigenVectors_(1,0) = 1;
  leftEigenVectors_(1,1) = -1;

  rightEigenVectors_(0,0) = 1;
  rightEigenVectors_(1,0) = 1;
  rightEigenVectors_(0,1) = 1-pi0_;
  rightEigenVectors_(1,1) = -1* pi0_;
}

/******************************************************************************/

double TwoParameterBinarySubstitutionModel::Pij_t(size_t i, size_t j, double d) const
{
  exp_ = exp(-lambda_ * rate_ * d);

  switch (i) {
    case 0:
      switch (j) {
        case 0: return ( (1-pi0_) + pi0_ * exp_);
        case 1: return pi0_ * (1 - exp_);
        default: return 0;
      }
    case 1:
      switch (j) {
        case 0: return (1-pi0_) * (1 - exp_);
        case 1: return pi0_ + (1-pi0_) * exp_;
        default: return 0;
      }
    default: return 0;
  }
}

/******************************************************************************/

double TwoParameterBinarySubstitutionModel::dPij_dt(size_t i, size_t j, double d) const
{
  exp_ = rate_ * exp(-lambda_ * rate_ * d);

  switch (i) {
    case 0:
      switch (j) {
        case 0: return -1 * pi0_ * exp_;
        case 1: return pi0_ * exp_;
        default: return 0;
      }
    case 1:
      switch (j) {
        case 0: return (1 - pi0_) * exp_;
        case 1: return -1 * (1 - pi0_) * exp_;
        default: return 0;
      }
    default: return 0;
  }
}

/******************************************************************************/

double TwoParameterBinarySubstitutionModel::d2Pij_dt2(size_t i, size_t j, double d) const
{
  exp_ = rate_ * rate_ * exp(-lambda_ * rate_ * d);

  switch (i) {
    case 0:
      switch (j) {
        case 0: return  pi0_ * exp_;
        case 1: return  -1 * pi0_ * exp_;
        default: return 0;
      }
    case 1:
      switch (j) {
        case 0: return -1 * (1 - pi0_) * exp_;
        case 1: return  (1 - pi0_) * exp_;
        default: return 0;
      }
    default: return 0;
  }
  return 0;
}

/******************************************************************************/

const Matrix<double>& TwoParameterBinarySubstitutionModel::getPij_t(double d) const
{
  exp_ = exp(-lambda_ * rate_ * d);

  p_(0,0) = (1 - pi0_) + pi0_ * exp_;
  p_(0,1) = pi0_ * (1 - exp_);

  p_(1,0) =  (1 - pi0_) * (1 - exp_);
  p_(1,1) = pi0_ + (1 - pi0_) * exp_;

  return p_;
}

/******************************************************************************/

const Matrix<double>& TwoParameterBinarySubstitutionModel::getdPij_dt(double d) const
{
  exp_ = rate_ * exp(-lambda_ * rate_ * d);

  p_(0,0) = -1 * pi0_ * exp_;
  p_(0,1) = pi0_ * exp_;

  p_(1,0) = (1 - pi0_) * exp_;
  p_(1,1) = -1 * (1 - pi0_) * exp_;

  return p_;
}

/******************************************************************************/

const Matrix<double>& TwoParameterBinarySubstitutionModel::getd2Pij_dt2(double d) const
{
  exp_ = rate_ * rate_ * exp(-lambda_ * rate_ * d);

  p_(0,0) = pi0_ * exp_;
  p_(0,1) = -1 * pi0_ * exp_;
  p_(1,0) = -1 * (1 - pi0_) * exp_;
  p_(1,1) = (1 - pi0_) * exp_;

  return p_;
}

/******************************************************************************/

void TwoParameterBinarySubstitutionModel::setMuBounds(double lb, double ub)
{
  std::shared_ptr<IntervalConstraint> bounds(new IntervalConstraint(lb, ub, true, true)); 
  getParameter_("mu").setConstraint(bounds);
}
