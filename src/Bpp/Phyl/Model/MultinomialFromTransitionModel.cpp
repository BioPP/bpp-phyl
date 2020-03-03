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

const Eigen::VectorXd& MultinomialFromTransitionModel::Lik_t(const Eigen::VectorXd& to, double t) const
{
  if (t!=tref_)
  {
    tref_=t;
    dPij_dt=0; d2Pij_dt2=0;
    Pij_t=&getTransitionModel().getPij_t(t);
  }
  else if (Pij_t==0)
    Pij_t=&getTransitionModel().getPij_t(t);
    
  for (size_t i=0; i< size_; i++)
  {
    Pi_(i)=0;
    for (size_t j=0;j<size_;j++)
      Pi_(i)+=(*Pij_t)(i,j)*to[j];
  }
  
  return Pi_;
}

const Eigen::VectorXd& MultinomialFromTransitionModel::dLik_dt(const Eigen::VectorXd& to, double t) const
{
  if (t!=tref_)
  {
    tref_=t;
    Pij_t=0; d2Pij_dt2=0;
    dPij_dt=&getTransitionModel().getdPij_dt(t);
  }
  else if (dPij_dt==0)
    dPij_dt=&getTransitionModel().getdPij_dt(t);

  for (size_t i=0; i< size_; i++)
  {
    dPi_(i)=0;
    for (size_t j=0;j<size_;j++)
      dPi_(i)+=(*dPij_dt)(i,j)*to[j];
  }
  return dPi_;
}

const Eigen::VectorXd& MultinomialFromTransitionModel::d2Lik_dt2(const Eigen::VectorXd& to, double t) const
{
  if (t!=tref_)
  {
    tref_=t;
    Pij_t=0; dPij_dt=0;
    d2Pij_dt2=&getTransitionModel().getd2Pij_dt2(t);
  }
  else if (d2Pij_dt2==0)
    d2Pij_dt2=&getTransitionModel().getd2Pij_dt2(t);

  for (size_t i=0; i< size_; i++)
  {
    d2Pi_(i)=0;
    for (size_t j=0;j<size_;j++)
      d2Pi_(i)+=(*d2Pij_dt2)(i,j)*to[j];
  }
  return d2Pi_;
}

