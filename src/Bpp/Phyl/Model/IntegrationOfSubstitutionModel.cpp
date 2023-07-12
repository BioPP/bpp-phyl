//
// File: IntegrationOfSubstitutionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: samedi 24 octobre 2015, Ã  18h 50
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#include "IntegrationOfSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

IntegrationOfSubstitutionModel::IntegrationOfSubstitutionModel(std::unique_ptr<SubstitutionModelInterface>& originalModel, uint k, uint n) :
  AbstractParameterAliasable("Integrate."),
  AbstractWrappedModel("Integrate."),
  AbstractWrappedTransitionModel("Integrate."),
  AbstractFromSubstitutionModelTransitionModel(move(originalModel), "Integrate."),
  svar_(n,1,true),
  vZeta_(n),
  k_(k),
  vMatInv_(),
  time_(-1),
  tmp_(model().getNumberOfStates(), model().getNumberOfStates()),
  tmp2_(model().getNumberOfStates(), model().getNumberOfStates()),
  tmp3_(model().getNumberOfStates(), model().getNumberOfStates())
{
  for (size_t i=0; i<n;i++)
    vMatInv_.push_back(RowMatrix<double>(model().getNumberOfStates(), model().getNumberOfStates()));

  svar_.setNamespace("Integrate.");

  auto lpar = svar_.getParameters();    
  addParameters_(lpar);
  vZeta_ = svar_.getFrequencies();
}

IntegrationOfSubstitutionModel::IntegrationOfSubstitutionModel(std::unique_ptr<SubstitutionModelInterface>& originalModel, uint k, const Vdouble& vZetas) :
  AbstractParameterAliasable("Integrate."),
  AbstractWrappedModel("Integrate."),
  AbstractWrappedTransitionModel("Integrate."),
  AbstractFromSubstitutionModelTransitionModel(move(originalModel), "Integrate."),
  svar_(vZetas.size(),1,true),
  vZeta_(vZetas),
  k_(k),
  vMatInv_(),
  time_(-1),
  tmp_(model().getNumberOfStates(), model().getNumberOfStates()),
  tmp2_(model().getNumberOfStates(), model().getNumberOfStates()),
  tmp3_(model().getNumberOfStates(), model().getNumberOfStates())
{
  for (size_t i=0; i<vZeta_.size();i++)
    vMatInv_.push_back(RowMatrix<double>(model().getNumberOfStates(), model().getNumberOfStates()));

  // set constraints to have decreasing values in the Simplex
  svar_.setNamespace("Integrate.");

  svar_.setFrequencies(vZeta_);
  
  auto lpar = svar_.getParameters();    
  addParameters_(lpar);
}

void IntegrationOfSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  AbstractFromSubstitutionModelTransitionModel::fireParameterChanged(parameters);

  svar_.matchParametersValues(parameters);
  vZeta_ = svar_.getFrequencies();

  model_().matchParametersValues(parameters);
}

void IntegrationOfSubstitutionModel::computeInv_(double t) const
{
  if (time_==t) // no need to recompute
    return;

  time_=t;
  
  const auto& Q = substitutionModel().getGenerator();
  auto rate = substitutionModel().getRate();
  
  size_t s = vZeta_.size();
  
  for (size_t n=0;n<s;n++)
  {
    double x = rate * (double) s * vZeta_[n]/k_*time_;
   
    for (size_t i=0;i<(size_t)getNumberOfStates();i++)
      for (size_t j=0;j<(size_t)getNumberOfStates();j++)
        vMatInv_[n](i,j)=((i==j)?1:0) - x * Q(i,j);

    MatrixTools::inv(vMatInv_[n],tmp_);
    MatrixTools::pow(tmp_,k_,vMatInv_[n]);
  }
}

double IntegrationOfSubstitutionModel::Pij_t(size_t i, size_t j, double t) const
{
  throw Exception("I do not want to be here");
  // computeInv_(t);
  
  // double x=0;
  // for (const auto& MatInv:vMatInv_)
  //   x += MatInv(i,j);

  // return x;
}
  

double IntegrationOfSubstitutionModel::dPij_dt  (size_t i, size_t j, double t) const
{
  
  throw Exception("I do not want to be here");
  // computeInv_(t);
  
  // double sum=0;
  // for (size_t n=0;n<vZeta_.size();n++)
  // {
  //   double x = vZeta_[n]/k_*time_;
   
  //     for (size_t i=0;i<(size_t)getNumberOfStates();i++)
  //       for (size_t j=0;j<(size_t)getNumberOfStates();j++)
  //         vMatInv_[n](i,j)=((i==j)?1:0) - x * Q(i,j);

  //     MatrixTools::inv(vMatInv_[n],tmp_);
  //     MatrixTools::pow(tmp_,k_,vMatInv_[n]);
  // }

  // sum += MatInv_[k](i,j);
  // return sum;
}

double IntegrationOfSubstitutionModel::d2Pij_dt2(size_t i, size_t j, double t) const
{
  throw Exception("I do not want to be here");
  // double x = 0;
  // for (size_t k=0;k<getNumberOfStates();k++)
  //   x += ImzkQpmk_(i,k) + substitutionModel().d2Pij_dt2(k,j,t);

  // return x;
}

const Matrix<double>& IntegrationOfSubstitutionModel::getPij_t(double t) const
{
  computeInv_(t);
  
  MatrixTools::fill(pijt_,0);

  for (const auto& MatInv:vMatInv_)
  {    
    MatrixTools::add(pijt_,MatInv);
  }

  MatrixTools::scale(pijt_,1./(double)vZeta_.size());
  return pijt_;
}


const Matrix<double>& IntegrationOfSubstitutionModel::getdPij_dt(double t) const
{
  computeInv_(t);
  
  const auto& Q = substitutionModel().getGenerator();
  auto rate = substitutionModel().getRate();

  MatrixTools::fill(tmp3_,0);
  
  size_t s = vZeta_.size();
  
  for (size_t n=0;n<s;n++)
  {
    double x = rate * (double) s * vZeta_[n]/k_*time_;
    
    for (size_t i=0;i<(size_t)getNumberOfStates();i++)
      for (size_t j=0;j<(size_t)getNumberOfStates();j++)
        tmp_(i,j)=((i==j)?1:0) - x * Q(i,j);

    MatrixTools::inv(tmp_, tmp2_);

    MatrixTools::mult(vMatInv_[n],tmp2_,tmp_);
    MatrixTools::scale(tmp_,vZeta_[n]);
    
    MatrixTools::add(tmp3_,tmp_);
  }

  MatrixTools::mult(Q,tmp3_,dpijt_);
  MatrixTools::scale(dpijt_,rate);
  
  return dpijt_;
}


const Matrix<double>& IntegrationOfSubstitutionModel::getd2Pij_dt2(double t) const
{
  computeInv_(t);
  
  const auto& Q = substitutionModel().getGenerator();
  auto rate = substitutionModel().getRate();

  MatrixTools::fill(tmp3_,0);

  size_t s = vZeta_.size();
  
  for (size_t n=0;n<s;n++)
  {
    double x = rate * (double) s * vZeta_[n]/k_*time_;
    
    for (size_t i=0;i<(size_t)getNumberOfStates();i++)
      for (size_t j=0;j<(size_t)getNumberOfStates();j++)
        tmp_(i,j)=((i==j)?1:0) - x * Q(i,j);

    MatrixTools::pow(tmp_, 2, tmp2_);
    MatrixTools::inv(tmp2_, tmp_);
    
    MatrixTools::mult(vMatInv_[n],tmp2_,tmp_);
    MatrixTools::scale(tmp_,vZeta_[n]*vZeta_[n]);
    
    MatrixTools::add(tmp3_,tmp_);
  }

  MatrixTools::pow(Q,2,tmp_);
  
  MatrixTools::mult(tmp_,tmp3_,d2pijt_);

  MatrixTools::scale(d2pijt_, rate * rate * (double) s * (k_+1)/k_);
  
  return d2pijt_;
}
