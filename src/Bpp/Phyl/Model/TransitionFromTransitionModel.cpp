// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "TransitionFromTransitionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

const Eigen::VectorXd& TransitionFromTransitionModel::Lik_t(const Eigen::VectorXd& to, double t) const
{
  if (t != tref_)
  {
    tref_ = t;
    dPij_dt = 0; d2Pij_dt2 = 0;
    Pij_t = &transitionModel().getPij_t(t);
  }
  else if (Pij_t == 0)
    Pij_t = &transitionModel().getPij_t(t);

  for (Eigen::Index i = 0; i < Eigen::Index(size_); ++i)
  {
    Pi_(i) = 0;
    for (Eigen::Index j = 0; j < Eigen::Index(size_); ++j)
    {
      Pi_(i) += (*Pij_t)(size_t(i), size_t(j)) * to[j];
    }
  }

  return Pi_;
}

const Eigen::VectorXd& TransitionFromTransitionModel::dLik_dt(const Eigen::VectorXd& to, double t) const
{
  if (t != tref_)
  {
    tref_ = t;
    Pij_t = 0; d2Pij_dt2 = 0;
    dPij_dt = &transitionModel().getdPij_dt(t);
  }
  else if (dPij_dt == 0)
    dPij_dt = &transitionModel().getdPij_dt(t);

  for (Eigen::Index i = 0; i < Eigen::Index(size_); ++i)
  {
    dPi_(i) = 0;
    for (auto j = 0; j < Eigen::Index(size_); ++j)
    {
      dPi_(i) += (*dPij_dt)(size_t(i), size_t(j)) * to[j];
    }
  }
  return dPi_;
}

const Eigen::VectorXd& TransitionFromTransitionModel::d2Lik_dt2(const Eigen::VectorXd& to, double t) const
{
  if (t != tref_)
  {
    tref_ = t;
    Pij_t = 0; dPij_dt = 0;
    d2Pij_dt2 = &transitionModel().getd2Pij_dt2(t);
  }
  else if (d2Pij_dt2 == 0)
    d2Pij_dt2 = &transitionModel().getd2Pij_dt2(t);

  for (Eigen::Index i = 0; i < Eigen::Index(size_); ++i)
  {
    d2Pi_(i) = 0;
    for (auto j = 0; j < Eigen::Index(size_); ++j)
    {
      d2Pi_(i) += (*d2Pij_dt2)(size_t(i), size_t(j)) * to[j];
    }
  }
  return d2Pi_;
}
