// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "MultinomialFromTransitionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

void MultinomialFromTransitionModel::compute_Multinomial_(const Eigen::VectorXd& counts) const
{
  Pi_.array() = 0;

  for (auto i = 0; i < Eigen::Index(size_); i++)
  {
    double* Pi_i = &Pi_(i);
    for (auto j = 0; j < Eigen::Index(size_); j++)
    {
      if (counts(j) != 0)
      {
        *Pi_i += counts(j) * std::log((*Pij_t)(size_t(i), size_t(j)));
      }
    }
  }

  Pi_.array() = (Pi_.array() + getFact_(counts)).exp();
}

void MultinomialFromTransitionModel::compute_dMultinomial_dt_(const Eigen::VectorXd& counts) const
{
  Pi_.array() = 0; dPi_.array() = 0;

  for (size_t i = 0; i < size_; i++)
  {
    double* Pi_i(&Pi_(Eigen::Index(i))), *dPi_i(&dPi_(Eigen::Index(i)));

    for (size_t j = 0; j < size_; j++)
    {
      if (counts(Eigen::Index(j)) != 0)
      {
        *Pi_i += counts(Eigen::Index(j)) * std::log((*Pij_t)(i, j));
        *dPi_i += counts(Eigen::Index(j)) * (*dPij_dt)(i, j) / (*Pij_t)(i, j);
      }
    }
  }

  dPi_.array() *= (Pi_.array() + getFact_(counts)).exp();
}

void MultinomialFromTransitionModel::compute_d2Multinomial_dt2_(const Eigen::VectorXd& counts) const
{
  Pi_.array() = 0; d2Pi_.array() = 0;

  for (size_t i = 0; i < size_; i++)
  {
    double* Pi_i(&Pi_(Eigen::Index(i))), *d2Pi_i(&d2Pi_(Eigen::Index(i)));

    for (size_t j = 0; j < size_; j++)
    {
      if (counts(Eigen::Index(j)) != 0)
      {
        *Pi_i += counts(Eigen::Index(j)) * std::log((*Pij_t)(i, j));
        *d2Pi_i += counts(Eigen::Index(j)) * (*d2Pij_dt2)(i, j) / (*Pij_t)(i, j);
      }
    }
  }

  d2Pi_.array() *= (Pi_.array() + getFact_(counts)).exp();
}


const Eigen::VectorXd& MultinomialFromTransitionModel::Lik_t(const Eigen::VectorXd& to, double t) const
{
  if (t != tref_)
  {
    tref_ = t;
    Pij_t = &transitionModel().getPij_t(t);
    dPij_dt = 0; d2Pij_dt2 = 0;
  }

  compute_Multinomial_(to);
  return Pi_;
}

const Eigen::VectorXd& MultinomialFromTransitionModel::dLik_dt(const Eigen::VectorXd& to, double t) const
{
  if (t != tref_)
  {
    tref_ = t;
    Pij_t = &transitionModel().getPij_t(t);
    dPij_dt = &transitionModel().getdPij_dt(t);
    d2Pij_dt2 = 0;
  }
  else if (dPij_dt == 0)
    dPij_dt = &transitionModel().getdPij_dt(t);

  compute_dMultinomial_dt_(to);

  return dPi_;
}

const Eigen::VectorXd& MultinomialFromTransitionModel::d2Lik_dt2(const Eigen::VectorXd& to, double t) const
{
  if (t != tref_)
  {
    tref_ = t;
    Pij_t = &transitionModel().getPij_t(t);
    dPij_dt = 0;
    d2Pij_dt2 = &transitionModel().getd2Pij_dt2(t);
  }
  else if (d2Pij_dt2 == 0)
    d2Pij_dt2 = &transitionModel().getd2Pij_dt2(t);

  compute_d2Multinomial_dt2_(to);

  return d2Pi_;
}
