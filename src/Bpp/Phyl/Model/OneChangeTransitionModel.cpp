//
// File: OneChangeTransitionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: samedi 24 octobre 2015, à 18h 50
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "OneChangeTransitionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

double OneChangeTransitionModel::Pij_t    (size_t i, size_t j, double t) const
{
  double qii = getSubstitutionModel().Qij(i, i);
  if (qii == 0)
  {
    return i == j ? 1 : 0;
  }
  else
  {
    if (t != 0)
    {
      double rate = getModel().getRate();
      double v = exp(qii * t * rate);
      return (getTransitionModel().Pij_t(i, j, t) - (i == j ? v : 0)) / (1 - v);
    }
    else
      return i == j ? 0 : -getSubstitutionModel().Qij(i, j) / qii;
  }
}

double OneChangeTransitionModel::dPij_dt  (size_t i, size_t j, double t) const
{
  const RowMatrix<double>& Q = getSubstitutionModel().getGenerator();
  double qii = Q(i, i);

  if (qii == 0)
    return 0;

  double rate = getTransitionModel().getRate();
  if (t != 0)
  {
    double v = exp(qii * t * rate);
    return (getTransitionModel().dPij_dt(i, j, t) * (1 - v) + (getTransitionModel().Pij_t(i, j, t) - (i == j ? 1 : 0)) * qii * rate * v) / ((1 - v) * (1 - v));
  }
  else
  {
    double q2ij = 0;
    for (size_t k = 0; k < size_; ++k)
    {
      q2ij += Q(i, k) * Q(k, j);
    }

    return rate * (-q2ij + qii * Q(i, j)) / (2 * qii);
  }
}

double OneChangeTransitionModel::d2Pij_dt2(size_t i, size_t j, double t) const
{
  double rate = getTransitionModel().getRate();
  if (t != 0)
  {
    double qii = getSubstitutionModel().Qij(i, i);

    if (qii == 0)
    {
      return 0;
    }
    else
    {
      double v = exp(qii * t * rate);
      double mv = 1 - v;

      return -((mv * getTransitionModel().d2Pij_dt2(i, j, t) + 2 * qii * rate * v * getTransitionModel().dPij_dt(i, j, t)) * mv + (getTransitionModel().Pij_t(i, j, t) - (i == j ? 1 : 0)) * qii * qii * rate * rate * (v + v * v)) / (mv * mv * mv);
    }
  }
  else
  {
    const RowMatrix<double>& Q = getSubstitutionModel().getGenerator();
    double qii = Q(i, i);

    if (qii == 0)
      return 0;

    double q2ik, q2ij = 0, q3ij = 0;

    for (size_t k = 0; k < size_; ++k)
    {
      q2ik = 0;
      for (size_t l = 0; l < size_; ++l)
      {
        q2ik += Q(i, l) * Q(l, k);
      }
      if (k == j)
        q2ij = q2ik;

      q3ij += q2ik * Q(k, j);
    }

    return -rate * rate * (2 * q3ij - 3 * qii * q2ij + qii * qii * Q(i, j)) / (12 * qii);
  }
}

const Matrix<double>& OneChangeTransitionModel::getPij_t(double t) const
{
  const RowMatrix<double>& origPij = getTransitionModel().getPij_t(t);
  const RowMatrix<double>& Q = getSubstitutionModel().getGenerator();
  double rate = getTransitionModel().getRate();

  for (unsigned int i = 0; i < size_; ++i)
  {
    vector<double>& pi_t = pij_t.getRow(i);

    double qii = Q(i, i);
    if (qii == 0)
    {
      for (unsigned int j = 0; j < size_; ++j)
      {
        pi_t[j] = (i == j ? 1 : 0);
      }
    }
    else
    {
      if (t != 0)
      {
        const vector<double>& origPij_i = origPij.getRow(i);

        double v = exp(qii * t * rate);
        for (unsigned int j = 0; j < size_; ++j)
        {
          pi_t[j] = (origPij_i[j] - (i == j ? v : 0)) / (1 - v);
        }
      }
      else
      {
        const vector<double>& Q_i = Q.getRow(i);
        for (unsigned int j = 0; j < size_; ++j)
        {
          pi_t[j] = (i == j ? 0 : -Q_i[j] / qii);
        }
      }
    }
  }
  return pij_t;
}


const Matrix<double>& OneChangeTransitionModel::getdPij_dt(double t) const
{
  double rate = getTransitionModel().getRate();
  if (t != 0)
  {
    const RowMatrix<double>& origPij = getTransitionModel().getPij_t(t);
    const RowMatrix<double>& origdPij = getTransitionModel().getdPij_dt(t);

    for (unsigned int i = 0; i < size_; ++i)
    {
      vector<double>& dpi_t = dpij_t.getRow(i);
      double qii = getSubstitutionModel().Qij(i, i);

      if (qii == 0)
      {
        for (unsigned int j = 0; j < size_; ++j)
        {
          dpi_t[j] = 0;
        }
      }
      else
      {
        {
          double v = exp(qii * t * rate);
          double mv2 = (1 - v) * (1 - v);
          const vector<double>& origPij_i = origPij.getRow(i);
          const vector<double>& origdPij_i = origdPij.getRow(i);

          for (unsigned int j = 0; j < size_; ++j)
          {
            dpi_t[j] = (origdPij_i[j] * (1 - v) + (origPij_i[j] - (i == j ? 1 : 0)) * qii * rate * v) / mv2;
          }
        }
      }
    }
  }
  else
  {
    const RowMatrix<double>& Q = getSubstitutionModel().getGenerator();
    RowMatrix<double> Q2;
    MatrixTools::mult<double>(Q, Q, Q2);

    for (unsigned int i = 0; i < size_; ++i)
    {
      const vector<double>& Q_i = Q.getRow(i);
      vector<double>& dpi_t = dpij_t.getRow(i);
      double qii = Q_i[i];

      if (qii == 0)
      {
        for (unsigned int j = 0; j < size_; ++j)
        {
          dpi_t[j] = 0;
        }
      }
      else
      {
        const vector<double>& Q2_i = Q2.getRow(i);

        for (unsigned int j = 0; j < size_; ++j)
        {
          dpi_t[j] = rate * (-Q2_i[j] + qii * Q_i[j]) / (2 * qii);
        }
      }
    }
  }

  return dpij_t;
}


const Matrix<double>& OneChangeTransitionModel::getd2Pij_dt2(double t) const
{
  double rate = getTransitionModel().getRate();
  double r2 = rate * rate;

  if (t != 0)
  {
    const RowMatrix<double>& origPij = getTransitionModel().getPij_t(t);
    const RowMatrix<double>& origdPij = getTransitionModel().getdPij_dt(t);
    const RowMatrix<double>& origd2Pij = getTransitionModel().getd2Pij_dt2(t);

    for (unsigned int i = 0; i < size_; ++i)
    {
      double qii = getSubstitutionModel().Qij(i, i);
      vector<double>& d2pi_t = d2pij_t.getRow(i);

      if (qii == 0)
      {
        for (unsigned int j = 0; j < size_; ++j)
        {
          d2pi_t[j] = 0;
        }
      }
      else
      {
        {
          double v = exp(rate * qii * t);
          double mv = 1 - v;
          double q2 = qii * qii;

          double mv3 = mv * mv * mv;
          double vpv2 = v + v * v;

          const vector<double>& origPij_i = origPij.getRow(i);
          const vector<double>& origdPij_i = origdPij.getRow(i);
          const vector<double>& origd2Pij_i = origd2Pij.getRow(i);

          for (unsigned int j = 0; j < size_; ++j)
          {
            d2pi_t[j] = -((mv * origd2Pij_i[j] + 2 * qii * rate * v * origdPij_i[j]) * mv + (origPij_i[j] - (i == j ? 1 : 0)) * q2 * r2 * vpv2) / mv3;
          }
        }
      }
    }
  }
  else
  {
    const RowMatrix<double>& Q = getSubstitutionModel().getGenerator();
    RowMatrix<double> Q2;
    MatrixTools::mult<double>(Q, Q, Q2);
    RowMatrix<double> Q3;
    MatrixTools::mult<double>(Q, Q2, Q3);

    for (unsigned int i = 0; i < size_; ++i)
    {
      const vector<double>& Q_i = Q.getRow(i);
      double qii = Q_i[i];
      vector<double>& d2pi_t = d2pij_t.getRow(i);

      if (qii == 0)
      {
        for (unsigned int j = 0; j < size_; ++j)
        {
          d2pi_t[j] = 0;
        }
      }
      else
      {
        const vector<double>& Q2_i = Q2.getRow(i);
        const vector<double>& Q3_i = Q3.getRow(i);

        for (unsigned int j = 0; j < size_; ++j)
        {
          d2pi_t[j] = -r2 * (2 * Q3_i[j] - 3 * qii * Q2_i[j] + qii * qii * Q_i[j]) / (12 * qii);
        }
      }
    }
  }
  return d2pij_t;
}
