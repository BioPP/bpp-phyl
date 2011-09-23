//
// File: BinarySubstitutionModel.cpp
// Created by: Laurent Gueguen
// Created on: 2009
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

#include "BinarySubstitutionModel.h"

// From the STL:
#include <cmath>

using namespace bpp;

#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace std;

/******************************************************************************/

BinarySubstitutionModel::BinarySubstitutionModel(const BinaryAlphabet* alpha, double kappa) :
  AbstractParameterAliasable("BinarySubstitutionModel."),
  AbstractSubstitutionModel(alpha, "BinarySubstitutionModel."),
  _kappa(kappa),
  _lambda(0),
  _exp(0),
  _p(size_,size_)
{
  Parameter kappap(getNamespace() + "kappa", _kappa, &Parameter::R_PLUS_STAR);
  addParameter_(kappap);
  updateMatrices();
}

/******************************************************************************/

void BinarySubstitutionModel::updateMatrices()
{
  _kappa = getParameterValue("kappa"); // alpha/beta
  _lambda = (_kappa + 1) * (_kappa + 1) / (2 * _kappa);

  // Frequences:
  freq_[0] = 1 / (_kappa + 1);
  freq_[1] = _kappa / (_kappa + 1);

  // Generator:
  generator_(0, 0) = rate_ * -(_kappa + 1) / 2;
  generator_(0, 1) = rate_ * (_kappa + 1) / 2;
  generator_(1, 0) = rate_ * (_kappa + 1) / (2 * _kappa);
  generator_(1, 1) = rate_ * -(_kappa + 1) / (2 * _kappa);

  // Eigen values:
  eigenValues_[0] = 0;
  eigenValues_[1] = -rate_ * _lambda;

  // Eigen vectors:
  leftEigenVectors_(0,0) = 1 / (_kappa + 1);
  leftEigenVectors_(0,1) = _kappa / (_kappa + 1);
  if (_kappa != 1.0)
  {
    leftEigenVectors_(1,0) = (_kappa - 1) / (_kappa + 1);
    leftEigenVectors_(1,1) = -(_kappa - 1) / (_kappa + 1);
  }
  else
  {
    leftEigenVectors_(1,0) = 1;
    leftEigenVectors_(1,1) = -1;
  }

  rightEigenVectors_(0,0) = 1;
  rightEigenVectors_(1,0) = 1;

  if (_kappa != 1.0)
  {
    rightEigenVectors_(0,1) = _kappa / (_kappa - 1);
    rightEigenVectors_(1,1) = -1 / (_kappa - 1);
  }
  else
  {
    rightEigenVectors_(0,1) = 1 / 2;
    rightEigenVectors_(1,1) = -1 / 2;
  }
}

/******************************************************************************/

double BinarySubstitutionModel::Pij_t(unsigned int i, unsigned int j, double d) const
{
  _exp = exp(-_lambda * d);

  switch (i)
  {
  case 0: {
    switch (j)
    {
    case 0: return (1 + _kappa * _exp) / (_kappa + 1);
    case 1: return _kappa / (_kappa + 1) * (1 - _exp);
    }
  }
  case 1: {
    switch (j)
    {
    case 0: return (1 - _exp) / (_kappa + 1);
    case 1: return (_kappa + _exp) / (_kappa + 1);
    }
  }
  }
  return 0;
}

/******************************************************************************/

double BinarySubstitutionModel::dPij_dt(unsigned int i, unsigned int j, double d) const
{
  _exp = exp(-_lambda * d);

  switch (i)
  {
  case 0: {
    switch (j)
    {
    case 0: return -(_kappa + 1) / 2 * _exp;
    case 1: return (_kappa + 1) / 2 * _exp;
    }
  }
  case 1: {
    switch (j)
    {
    case 0: return (_kappa + 1) / (2 * _kappa) * _exp;
    case 1: return -(_kappa + 1) / (2 * _kappa) * _exp;
    }
  }
  }
  return 0;
}

/******************************************************************************/

double BinarySubstitutionModel::d2Pij_dt2(unsigned int i, unsigned int j, double d) const
{
  _exp = exp(-_lambda * d);

  switch (i)
  {
  case 0: {
    switch (j)
    {
    case 0: return _lambda * (_kappa + 1) / 2 * _exp;
    case 1: return -_lambda * (_kappa + 1) / 2 * _exp;
    }
  }
  case 1: {
    switch (j)
    {
    case 0: return -_lambda * (_kappa + 1) / (2 * _kappa) * _exp;
    case 1: return _lambda * (_kappa + 1) / (2 * _kappa) * _exp;
    }
  }
  }
  return 0;
}

/******************************************************************************/

const Matrix<double>& BinarySubstitutionModel::getPij_t(double d) const
{
  _exp = exp(-_lambda * d);

  _p(0,0) = (1 + _kappa * _exp) / (_kappa + 1);
  _p(0,1) = _kappa / (_kappa + 1) * (1 - _exp);

  _p(1,0) =  (1 - _exp) / (_kappa + 1);
  _p(1,1) = (_kappa + _exp) / (_kappa + 1);

  return _p;
}

const Matrix<double>& BinarySubstitutionModel::getdPij_dt(double d) const
{
  _exp = exp(-_lambda * d);

  _p(0,0) = -(_kappa + 1) / 2 * _exp;
  _p(0,1) = (_kappa + 1) / 2 * _exp;

  _p(1,0) = (_kappa + 1) / (2 * _kappa) * _exp;
  _p(1,1) = -(_kappa + 1) / (2 * _kappa) * _exp;

  return _p;
}

const Matrix<double>& BinarySubstitutionModel::getd2Pij_dt2(double d) const
{
  _exp = exp(-_lambda * d);

  _p(0,0) = _lambda * (_kappa + 1) / 2 * _exp;
  _p(0,1) = -_lambda * (_kappa + 1) / 2 * _exp;
  _p(1,0) = -_lambda * (_kappa + 1) / (2 * _kappa) * _exp;
  _p(1,1) = _lambda * (_kappa + 1) / (2 * _kappa) * _exp;

  return _p;
}

/******************************************************************************/

void BinarySubstitutionModel::setFreq(std::map<int, double>& freqs)
{
  _kappa = freqs[1] / freqs[0];
  setParameterValue("kappa",_kappa);
  updateMatrices();
}
