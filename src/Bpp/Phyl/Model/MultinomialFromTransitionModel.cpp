//
// File: MultinomialFromTransitionModel.cpp
// Created by: Laurent Gueguen
// Created on: dimanche 26 janvier 2020, à 07h 55
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

#include "MultinomialFromTransitionModel.h"
#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

void MultinomialFromTransitionModel::compute_Multinomial_(const Eigen::VectorXd& counts) const
{
  Pi_.array()=0;
  
  for (auto i=0; i< Eigen::Index(size_); i++)
  {
    double* Pi_i = &Pi_(i);
    for (auto j=0;j<Eigen::Index(size_);j++)
    {
      if (counts(j)!=0)
      {
        *Pi_i += counts(j)*std::log((*Pij_t)(size_t(i),size_t(j)));
      }
    }
  }
  
  Pi_.array() = Pi_.array().exp() * getFact_(counts);
}

void MultinomialFromTransitionModel::compute_dMultinomial_dt_(const Eigen::VectorXd& counts) const
{
  Pi_.array()=0; dPi_.array()=0;
  
  for (size_t i=0; i< size_; i++)
  {
    double *Pi_i(&Pi_(Eigen::Index(i))), *dPi_i(&dPi_(Eigen::Index(i)));
    
    for (size_t j=0;j<size_;j++)
    {
      if (counts(Eigen::Index(j))!=0)
      {
        *Pi_i += counts(Eigen::Index(j))*std::log((*Pij_t)(i,j));
        *dPi_i += counts(Eigen::Index(j))*(*dPij_dt)(i,j)/(*Pij_t)(i,j);
      }
    }
  }
  
  dPi_.array() *= Pi_.array().exp()  * getFact_(counts);
}

void MultinomialFromTransitionModel::compute_d2Multinomial_dt2_(const Eigen::VectorXd& counts) const
{
  Pi_.array()=0; d2Pi_.array()=0;
  
  for (size_t i=0; i< size_; i++)
  {
    double *Pi_i(&Pi_(Eigen::Index(i))), *d2Pi_i(&d2Pi_(Eigen::Index(i)));
    
    for (size_t j=0;j<size_;j++)
    {
      if (counts(Eigen::Index(j))!=0)
      {
        *Pi_i += counts(Eigen::Index(j))*std::log((*Pij_t)(i,j));
        *d2Pi_i += counts(Eigen::Index(j))*(*d2Pij_dt2)(i,j)/(*Pij_t)(i,j);
      }
    }
  }
  
  d2Pi_.array() *= Pi_.array().exp() * getFact_(counts);
}


const Eigen::VectorXd& MultinomialFromTransitionModel::Lik_t(const Eigen::VectorXd& to, double t) const
{
  if (t!=tref_)
  {
    tref_=t;
    Pij_t=&getTransitionModel().getPij_t(t);
    dPij_dt=0; d2Pij_dt2=0;
  }

  compute_Multinomial_(to);
  return Pi_;
}

const Eigen::VectorXd& MultinomialFromTransitionModel::dLik_dt(const Eigen::VectorXd& to, double t) const
{
  if (t!=tref_)
  {
    tref_=t;
    Pij_t=&getTransitionModel().getPij_t(t);
    dPij_dt=&getTransitionModel().getdPij_dt(t);
    d2Pij_dt2=0;
  }
  else if (dPij_dt==0)
    dPij_dt=&getTransitionModel().getdPij_dt(t);

  compute_dMultinomial_dt_(to);

  return dPi_;
}

const Eigen::VectorXd& MultinomialFromTransitionModel::d2Lik_dt2(const Eigen::VectorXd& to, double t) const
{
  if (t!=tref_)
  {
    tref_=t;
    Pij_t=&getTransitionModel().getPij_t(t);
    dPij_dt=0;
    d2Pij_dt2=&getTransitionModel().getd2Pij_dt2(t);
  }
  else if (d2Pij_dt2==0)
    d2Pij_dt2=&getTransitionModel().getd2Pij_dt2(t);

  compute_d2Multinomial_dt2_(to);

  return d2Pi_;
}

