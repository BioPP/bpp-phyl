//
// File: JCnuc.cpp
// Created by: Julien Dutheil
// Created on: Tue May 27 16:04:36 2003
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

#include "JCnuc.h"

using namespace bpp;

#include <cmath>

using namespace std;

/******************************************************************************/

JCnuc::JCnuc(const NucleicAlphabet* alpha) :
  AbstractParameterAliasable("JC69."),
  AbstractReversibleNucleotideSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "JC69."),
  exp_(),
  p_(size_, size_)
{
  updateMatrices();
}

/******************************************************************************/

void JCnuc::updateMatrices()
{
  // Frequencies:
  freq_[0] = freq_[1] = freq_[2] = freq_[3] = 1. / 4.;

  // Generator and exchangeabilities:
  for (size_t i = 0; i < 4; ++i)
  {
    for (size_t j = 0; j < 4; ++j)
    {
      generator_(i, j) = (i == j) ? -1. : 1. / 3.;
      exchangeability_(i, j) = generator_(i, j) * 4.;
    }
  }

  // Eigen values:
  eigenValues_[0] = 0;
  eigenValues_[1] = eigenValues_[2] = eigenValues_[3] = -4. / 3.;

  // Eigen vectors:
  for (size_t i = 0; i < 4; i++)
  {
    leftEigenVectors_(0, i) = 1. / 4.;
  }
  for (size_t i = 1; i < 4; i++)
  {
    for (size_t j = 0; j < 4; j++)
    {
      leftEigenVectors_(i, j) = -1. / 4.;
    }
  }
  leftEigenVectors_(1, 2) = 3. / 4.;
  leftEigenVectors_(2, 1) = 3. / 4.;
  leftEigenVectors_(3, 0) = 3. / 4.;

  for (size_t i = 0; i < 4; i++)
  {
    rightEigenVectors_(i, 0) = 1.;
  }
  for (size_t i = 1; i < 4; i++)
  {
    rightEigenVectors_(3, i) = -1.;
  }
  for (size_t i = 0; i < 3; i++)
  {
    for (size_t j = 1; j < 4; j++)
    {
      rightEigenVectors_(i, j) = 0.;
    }
  }
  rightEigenVectors_(2, 1) = 1.;
  rightEigenVectors_(1, 2) = 1.;
  rightEigenVectors_(0, 3) = 1.;
}

/******************************************************************************/

double JCnuc::Pij_t(size_t i, size_t j, double d) const
{
  if (i == j)
    return 1. / 4. + 3. / 4. * exp(-rate_ * 4. / 3. * d);
  else
    return 1. / 4. - 1. / 4. * exp(-rate_ * 4. / 3. * d);
}

/******************************************************************************/

double JCnuc::dPij_dt(size_t i, size_t j, double d) const
{
  if (i == j)
    return -exp(-rate_ * 4. / 3. * d) * rate_;
  else
    return 1. / 3. * exp(-rate_ * 4. / 3. * d) * rate_;
}

/******************************************************************************/

double JCnuc::d2Pij_dt2(size_t i, size_t j, double d) const
{
  if (i == j)
    return 4. / 3. * exp(-rate_ * 4. / 3. * d) * rate_ * rate_;
  else
    return -4. / 9. * exp(-rate_ * 4. / 3. * d) * rate_ * rate_;
}

/******************************************************************************/

const Matrix<double>& JCnuc::getPij_t(double d) const
{
  exp_ = exp(-4. / 3. * d * rate_);
  for (size_t i = 0; i < size_; i++)
  {
    for (size_t j = 0; j < size_; j++)
    {
      p_(i, j) = (i == j) ? 1. / 4. + 3. / 4. * exp_ : 1. / 4. - 1. / 4. * exp_;
    }
  }
  return p_;
}

const Matrix<double>& JCnuc::getdPij_dt(double d) const
{
  exp_ = exp(-4. / 3. * d * rate_);
  for (size_t i = 0; i < size_; i++)
  {
    for (size_t j = 0; j < size_; j++)
    {
      p_(i, j) = rate_ * ((i == j) ? -exp_ : 1. / 3. * exp_);
    }
  }
  return p_;
}

const Matrix<double>& JCnuc::getd2Pij_dt2(double d) const
{
  exp_ = exp(-4. / 3. * d * rate_);
  for (size_t i = 0; i < size_; i++)
  {
    for (size_t j = 0; j < size_; j++)
    {
      p_(i, j) = rate_ * rate_ * ((i == j) ? 4. / 3. * exp_ : -4. / 9. * exp_);
    }
  }
  return p_;
}

/******************************************************************************/

