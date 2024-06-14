// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
  
  const auto& Q = substitutionModel().generator();
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
  
  const auto& Q = substitutionModel().generator();
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
  
  const auto& Q = substitutionModel().generator();
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
