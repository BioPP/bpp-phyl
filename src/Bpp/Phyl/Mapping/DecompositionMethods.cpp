//
// File: DecompositionMethods.cpp
// Created by: Laurent Guéguen
// Created on: mardi 18 juillet 2017, à 22h 44
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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

#include "DecompositionMethods.h"

#include <vector>
#include <typeinfo>

using namespace bpp;
using namespace std;



/******************************************************************************/

DecompositionMethods::DecompositionMethods(const SubstitutionModel* model, SubstitutionRegister* reg) :
  model_(model),
  nbStates_(model->getNumberOfStates()),
  nbTypes_(reg->getNumberOfSubstitutionTypes()),
  jMat_(nbStates_, nbStates_),
  vjMat_(0), // will be initialized only if needed
  bMatrices_(reg->getNumberOfSubstitutionTypes()),
  insideProducts_(reg->getNumberOfSubstitutionTypes())
{
}				


DecompositionMethods::DecompositionMethods(const SubstitutionModel* model) :
  model_(model),
  nbStates_(model->getNumberOfStates()),
  nbTypes_(1),
  jMat_(nbStates_, nbStates_),
  vjMat_(0), // will be initialized only if needed
  bMatrices_(1),
  insideProducts_(1)
{
}				

void DecompositionMethods::computeProducts_()
{
  for (size_t i = 0; i < nbTypes_; ++i) {
    //vInv_ %*% bMatrices_[i] %*% v_;
    RowMatrix<double> tmp(nbStates_, nbStates_);
    MatrixTools::mult(model_->getRowLeftEigenVectors(), bMatrices_[i], tmp);
    MatrixTools::mult(tmp, model_->getColumnRightEigenVectors(), insideProducts_[i]);
  }
}

void DecompositionMethods::initStates_()
{
  jMat_.resize(nbStates_, nbStates_);
  vjMat_.clear();
}

void DecompositionMethods::initBMatrices_()
{
  //Re-initialize all B matrices according to substitution register.
  bMatrices_.resize(nbTypes_);
  insideProducts_.resize(nbTypes_);

  for (size_t i = 0; i < nbTypes_; ++i) {
    bMatrices_[i].resize(nbStates_, nbStates_);
    insideProducts_[i].resize(nbStates_, nbStates_);
  }
}


void DecompositionMethods::jFunction_(const std::vector<double>& lambda, double t, RowMatrix<double>& result) const
{
  vector<double> expLam = VectorTools::exp(lambda * t);
  for (unsigned int i = 0; i < nbStates_; ++i) {
    for (unsigned int j = 0; j < nbStates_; ++j) {
      double dd = lambda[i] - lambda[j];
      if (dd == 0) {
        result(i, j) = t * expLam[i];
      } else {
        result(i, j) = (expLam[i] - expLam[j]) / dd;
      }
    }
  }
}

void DecompositionMethods::computeExpectations(RowMatrix<double>& mapping, double length) const
{
  if (model_->isDiagonalizable())
  {
    jFunction_(model_->getEigenValues(), length, jMat_);
    
    RowMatrix<double> tmp1(nbStates_, nbStates_), tmp2(nbStates_, nbStates_);
    MatrixTools::hadamardMult(jMat_, insideProducts_[0], tmp1);
    MatrixTools::mult(model_->getColumnRightEigenVectors(), tmp1, tmp2);
    MatrixTools::mult(tmp2, model_->getRowLeftEigenVectors(), mapping);
  }
  else if (model_->isNonSingular())
  {
    vector<vector<RowMatrix<double> > > tmp(3);
    for (size_t i=0;i<3;i++)
    {
      tmp[i].resize(3);
      for (size_t j=0;j<3;j++)
        tmp[i][j].resize(nbStates_, nbStates_);
    }
    
    const RowMatrix<double>& Ui=model_->getRowLeftEigenVectors();
    const RowMatrix<double>& U=model_->getColumnRightEigenVectors();

    // Matrices "row translated" up and down
    
    RowMatrix<double> Uip(nbStates_, nbStates_), Uim(nbStates_, nbStates_);
    RowMatrix<double> insidep(nbStates_, nbStates_), insidem(nbStates_, nbStates_);

    MatrixTools::copyUp(Ui, Uip);
    MatrixTools::copyDown(Ui, Uim);
    
    jFunction_(model_->getEigenValues(), model_->getIEigenValues(), length, vjMat_);
    
    double un(-1.);
    
    const RowMatrix<double>& inside(insideProducts_[0]);      
    MatrixTools::copyUp(inside, insidep);
    MatrixTools::copyDown(inside, insidem);

    MatrixTools::hadamardMult(vjMat_[0][0], inside,  tmp[0][0]);
    MatrixTools::hadamardMult(vjMat_[0][1], insidem, tmp[0][1]);
    MatrixTools::hadamardMult(vjMat_[0][2], insidep, tmp[0][2]);
    
    MatrixTools::add(tmp[0][0],un,tmp[0][1]);
    MatrixTools::add(tmp[0][0],tmp[0][2]);
    
    MatrixTools::hadamardMult(vjMat_[1][0], inside,  tmp[1][0]);
    MatrixTools::hadamardMult(vjMat_[1][1], insidem, tmp[1][1]);
    MatrixTools::hadamardMult(vjMat_[1][2], insidep, tmp[1][2]);
    
    
    MatrixTools::scale(tmp[1][0],un);      
    MatrixTools::add(tmp[1][0],tmp[1][1]);
    MatrixTools::add(tmp[1][0],un,tmp[1][2]);
    
    MatrixTools::hadamardMult(vjMat_[2][0], inside,  tmp[2][0]);
    MatrixTools::hadamardMult(vjMat_[2][1], insidem, tmp[2][1]);
    MatrixTools::hadamardMult(vjMat_[2][2], insidep, tmp[2][2]);
    
    MatrixTools::add(tmp[2][0],un,tmp[2][1]);
    MatrixTools::add(tmp[2][0],tmp[2][2]);
    

    MatrixTools::mult(tmp[0][0], Ui, tmp[0][1]);
    MatrixTools::mult(tmp[1][0], Uim,tmp[1][1]);
    MatrixTools::mult(tmp[2][0], Uip,tmp[2][1]);
    
    
    MatrixTools::add(tmp[0][1],tmp[1][1]);
    MatrixTools::add(tmp[0][1],tmp[2][1]);

    MatrixTools::mult(model_->getColumnRightEigenVectors(), tmp[0][1], mapping);
  }
  else
    throw Exception("void DecompositionMethods::computeMappings : substitution mapping is not implemented for singular generators.");
}


void DecompositionMethods::computeExpectations(std::vector< RowMatrix<double> >& mappings, double length) const
{
  if (model_->isDiagonalizable())
  {
    jFunction_(model_->getEigenValues(), length, jMat_);
    
    for (size_t i = 0; i < nbTypes_; ++i) {
      RowMatrix<double> tmp1(nbStates_, nbStates_), tmp2(nbStates_, nbStates_);
      MatrixTools::hadamardMult(jMat_, insideProducts_[i], tmp1);
      MatrixTools::mult(model_->getColumnRightEigenVectors(), tmp1, tmp2);
      MatrixTools::mult(tmp2, model_->getRowLeftEigenVectors(), mappings[i]);
    }
  }
  else if (model_->isNonSingular())
  {
    vector<vector<RowMatrix<double> > > tmp(3);
    for (size_t i=0;i<3;i++)
    {
      tmp[i].resize(3);
      for (size_t j=0;j<3;j++)
        tmp[i][j].resize(nbStates_, nbStates_);
    }
    
    const RowMatrix<double>& Ui=model_->getRowLeftEigenVectors();
    const RowMatrix<double>& U=model_->getColumnRightEigenVectors();

    // Matrices "row translated" up and down
    
    RowMatrix<double> Uip(nbStates_, nbStates_), Uim(nbStates_, nbStates_);
    RowMatrix<double> insidep(nbStates_, nbStates_), insidem(nbStates_, nbStates_);

    MatrixTools::copyUp(Ui, Uip);
    MatrixTools::copyDown(Ui, Uim);
    
    jFunction_(model_->getEigenValues(), model_->getIEigenValues(), length, vjMat_);

    double un(-1.);
    
    for (size_t i = 0; i < nbTypes_; ++i) {

      const RowMatrix<double>& inside(insideProducts_[i]);      
      MatrixTools::copyUp(inside, insidep);
      MatrixTools::copyDown(inside, insidem);

      MatrixTools::hadamardMult(vjMat_[0][0], inside,  tmp[0][0]);
      MatrixTools::hadamardMult(vjMat_[0][1], insidem, tmp[0][1]);
      MatrixTools::hadamardMult(vjMat_[0][2], insidep, tmp[0][2]);

      MatrixTools::add(tmp[0][0],un,tmp[0][1]);
      MatrixTools::add(tmp[0][0],tmp[0][2]);

      MatrixTools::hadamardMult(vjMat_[1][0], inside,  tmp[1][0]);
      MatrixTools::hadamardMult(vjMat_[1][1], insidem, tmp[1][1]);
      MatrixTools::hadamardMult(vjMat_[1][2], insidep, tmp[1][2]);

      
      MatrixTools::scale(tmp[1][0],un);      
      MatrixTools::add(tmp[1][0],tmp[1][1]);
      MatrixTools::add(tmp[1][0],un,tmp[1][2]);

      MatrixTools::hadamardMult(vjMat_[2][0], inside,  tmp[2][0]);
      MatrixTools::hadamardMult(vjMat_[2][1], insidem, tmp[2][1]);
      MatrixTools::hadamardMult(vjMat_[2][2], insidep, tmp[2][2]);

      MatrixTools::add(tmp[2][0],un,tmp[2][1]);
      MatrixTools::add(tmp[2][0],tmp[2][2]);

      
      MatrixTools::mult(tmp[0][0], Ui, tmp[0][1]);
      MatrixTools::mult(tmp[1][0], Uim,tmp[1][1]);
      MatrixTools::mult(tmp[2][0], Uip,tmp[2][1]);
      
                        
      MatrixTools::add(tmp[0][1],tmp[1][1]);
      MatrixTools::add(tmp[0][1],tmp[2][1]);
      
      MatrixTools::mult(model_->getColumnRightEigenVectors(), tmp[0][1], mappings[i]);
    }
    
  }
  else
    throw Exception("void DecompositionMethods::computeMappings : substitution mapping is not implemented for singular generators.");
}


/******************************************************************************/

double DecompositionMethods::computeJcc_(const std::vector<double>& fact, const std::vector<double>& vfonc, double t) const
{
  double alpha(fact[0]),beta(fact[1]),gamma(fact[2]),delta(fact[3]);
  double epsilon(alpha-gamma);
  double f1(vfonc[0]*vfonc[1]-vfonc[3]*vfonc[4]),f2(-vfonc[3]*vfonc[5]), f3(vfonc[0]*vfonc[2]);
  
  if (epsilon!=0)
  {
    double beta2(beta*beta),delta2(delta*delta),eps2(epsilon*epsilon);

    return (epsilon*(eps2+delta2+beta2)*f1
            -delta*(eps2+delta2-beta2)*f2
            +beta*(eps2-delta2+beta2)*f3)/((eps2+(beta+delta)*(beta+delta))*(eps2+(beta-delta)*(beta-delta)));
  }
  else
  {
    if (abs(beta)!=abs(delta))
      return (-delta*f2-beta*f3)/(delta*delta-beta*beta);
    else
    {
      if (delta!=0)
        return vfonc[0]*(vfonc[5]/delta+t*vfonc[4])/2;
      else
        return t*vfonc[0];
    }
  }
}

double DecompositionMethods::computeJcs_(const std::vector<double>& fact, const std::vector<double>& vfonc, double t) const
{
  double alpha(fact[0]),beta(fact[1]),gamma(fact[2]),delta(fact[3]);
  double epsilon(alpha-gamma);
  double f1(vfonc[0]*vfonc[1]-vfonc[3]*vfonc[4]),f2(-vfonc[3]*vfonc[5]), f3(vfonc[0]*vfonc[2]);
  
  if (epsilon!=0)
  {
    double beta2(beta*beta),delta2(delta*delta),eps2(epsilon*epsilon);

    return (delta*(eps2+delta2-beta2)*f1
            +epsilon*(eps2+delta2+beta2)*f2
            +2*beta*delta*epsilon*f3)/((eps2+(beta+delta)*(beta+delta))*(eps2+(beta-delta)*(beta-delta)));

  }
  else
  {
    if (abs(beta)!=abs(delta))
      return delta*f1/(delta*delta-beta*beta);
    else
    {
      if (delta!=0)
        return vfonc[0]*t*vfonc[5]/2;
      else
        return 0;
    }
  }
}

double DecompositionMethods::computeJsc_(const std::vector<double>& fact, const std::vector<double>& vfonc, double t) const
{
  double alpha(fact[0]),beta(fact[1]),gamma(fact[2]),delta(fact[3]);
  double epsilon(alpha-gamma);
  double f1(vfonc[0]*vfonc[1]-vfonc[3]*vfonc[4]),f2(-vfonc[3]*vfonc[5]), f3(vfonc[0]*vfonc[2]);
  
  if (epsilon!=0)
  {
    double beta2(beta*beta),delta2(delta*delta),eps2(epsilon*epsilon);

    return (-beta*(eps2-delta2+beta2)*f1
            +2*beta*delta*epsilon*f2
            +epsilon*(eps2+delta2+beta2)*f3)/((eps2+(beta+delta)*(beta+delta))*(eps2+(beta-delta)*(beta-delta)));
  }
  else
  {
    if (abs(beta)!=abs(delta))
      return beta*f1/(delta*delta-beta*beta);
    else
    {
      if (delta!=0)
        return vfonc[0]*t*vfonc[5]/2;
      else
        return 0;
    }
  }
}

double DecompositionMethods::computeJss_(const std::vector<double>& fact, const std::vector<double>& vfonc, double t) const
{
  double alpha(fact[0]),beta(fact[1]),gamma(fact[2]),delta(fact[3]);
  double epsilon(alpha-gamma);
  double f1(vfonc[0]*vfonc[1]-vfonc[3]*vfonc[4]),f2(-vfonc[3]*vfonc[5]), f3(vfonc[0]*vfonc[2]);
  
  if (epsilon!=0)
  {
    double beta2(beta*beta),delta2(delta*delta),eps2(epsilon*epsilon);

    return (-2*beta*delta*epsilon*f1
            -beta*(eps2-delta2+beta2)*f2
            +delta*(eps2+delta2-beta2)*f3)/((eps2+(beta+delta)*(beta+delta))*(eps2+(beta-delta)*(beta-delta)));
  }
  else
  {
    if (abs(beta)!=abs(delta))
      return (-beta*f2+delta*f3)/(delta*delta-beta*beta);
    else
    {
      if (delta!=0)
         return vfonc[0]*(vfonc[5]/delta-t*vfonc[4])/2;
      else
        return 0;
    }
  }
}

double DecompositionMethods::computeJc_(const std::vector<double>& fact, const std::vector<double>& vfonc, double t) const
{
  double alpha(fact[0]),gamma(fact[1]),delta(fact[2]);
  double epsilon(alpha-gamma);
  double f1(vfonc[0]-vfonc[1]*vfonc[2]),f2(vfonc[1]*vfonc[3]);
  
  if (epsilon!=0)
    return (epsilon*f1+delta*f2)/(delta*delta+epsilon*epsilon);
  else
  {
    if (delta!=0)
      return f2/delta;
    else
      return t*vfonc[0];
  }
}

double DecompositionMethods::computeJs_(const std::vector<double>& fact, const std::vector<double>& vfonc, double t) const
{
  double alpha(fact[0]),gamma(fact[1]),delta(fact[2]);
  double epsilon(alpha-gamma);
  double f1(vfonc[0]-vfonc[1]*vfonc[2]),f2(vfonc[1]*vfonc[3]);
  
  if (epsilon!=0)
    return (delta*f1-epsilon*f2)/(delta*delta+epsilon*epsilon);
  else
  {
    if (delta!=0)
      return f1/delta;
    else
      return 0;
  }
}

double DecompositionMethods::computeKc_(const std::vector<double>& fact, const std::vector<double>& vfonc, double t) const
{
  double alpha(fact[0]),beta(fact[1]),gamma(fact[2]);
  double epsilon(alpha-gamma);
  double f3(vfonc[0]*vfonc[1]-vfonc[3]), f4(vfonc[0]*vfonc[2]);
  
  if (epsilon!=0)
    return (epsilon*f3+beta*f4)/(beta*beta+epsilon*epsilon);
  else
  {
    if (beta!=0)
      return f4/beta;
    else
      return t*vfonc[0];
  }
}

double DecompositionMethods::computeKs_(const std::vector<double>& fact, const std::vector<double>& vfonc, double t) const
{
  double alpha(fact[0]),beta(fact[1]),gamma(fact[2]);
  double epsilon(alpha-gamma);
  double f3(vfonc[0]*vfonc[1]-vfonc[3]), f4(vfonc[0]*vfonc[2]);
  
  if (epsilon!=0)
    return (-beta*f3+epsilon*f4)/(beta*beta+epsilon*epsilon);
  else
  {
    if (beta!=0)
      return -f3/beta;
    else
      return 0;
  }
}

double DecompositionMethods::computeD_(const std::vector<double>& fact, const std::vector<double>& vfonc, double t) const
{
  double dd = fact[1]-fact[0];
  if (dd == 0) 
    return t * vfonc[0];
  else 
    return (vfonc[1] - vfonc[0]) / dd;
}

void DecompositionMethods::jFunction_(const std::vector<double>& lambda, const std::vector<double>& ilambda, double t, vector< vector<RowMatrix<double> > >& vresult) const
{
  // Initialized here if needed
  if (vresult.size()!=3)
  {
    vresult.resize(3);
    for (size_t i=0;i<3;i++)
    {
      vresult[i].resize(3);
      for (size_t j=0;j<3;j++)
        vresult[i][j].resize(nbStates_, nbStates_);
    }
  }
  
  vector<double> expLam = VectorTools::exp(lambda * t);
  vector<double> cosLam(nbStates_), sinLam(nbStates_);

  vector<double> fact_(4), vfunc_(6);
  
  for (unsigned int i = 0; i < nbStates_; ++i) {
    if (ilambda[i]!=0)
    {
      cosLam[i]=cos(ilambda[i]*t);
      sinLam[i]=sin(ilambda[i]*t);
    }
  }

  // 1st block
  // 1st matrice

  for (unsigned int i = 0; i < nbStates_; ++i) {
    vfunc_[0]=expLam[i];
    fact_[0]=lambda[i];
    if (ilambda[i]==0)
    {
      for (unsigned int j = 0; j < nbStates_; ++j) {
        vfunc_[1]=expLam[j];
        fact_[1]=lambda[j];
        if (ilambda[j]==0)
          vresult[0][0](i, j) = computeD_(fact_, vfunc_, t);
        else
        {
          vfunc_[2]=cosLam[j];
          vfunc_[3]=sinLam[j];
          fact_[2]=ilambda[j];
          vresult[0][0](i,j) = computeJc_(fact_, vfunc_, t);
        }
      }
    }
    else
    {
      vfunc_[1]=cosLam[i];
      vfunc_[2]=sinLam[i];
      fact_[1]=ilambda[i];

      for (unsigned int j = 0; j < nbStates_; ++j) {
        vfunc_[3]=expLam[j];
        fact_[2]=lambda[j];
        if (ilambda[j]==0)
          vresult[0][0](i, j) = computeKc_(fact_, vfunc_, t);
        else
        {
          vfunc_[4]=cosLam[j];
          vfunc_[5]=sinLam[j];
          fact_[3]=ilambda[j];
          vresult[0][0](i,j) = computeJcc_(fact_, vfunc_, t);
        }
      }
    }
  }

  // 2nd matrice

  for (unsigned int i = 0; i < nbStates_; ++i) {
    if (ilambda[i]>=0)
    {
      for (unsigned int j = 0; j < nbStates_; ++j) 
        vresult[0][1](i,j)=0;
    }
    else
    {
      vfunc_[0]=expLam[i];
      vfunc_[1]=cosLam[i];
      vfunc_[2]=sinLam[i];
      fact_[0]=lambda[i];
      fact_[1]=ilambda[i];
      for (unsigned int j = 0; j < nbStates_; ++j) 
      {
        vfunc_[3]=expLam[j];
        fact_[2]=lambda[j];
        
        if (ilambda[j]==0)
          vresult[0][1](i, j) = computeKs_(fact_, vfunc_, t);
        else
        {
          vfunc_[4]=cosLam[j];
          vfunc_[5]=sinLam[j];
          fact_[3]=ilambda[j];
          vresult[0][1](i,j) = computeJsc_(fact_, vfunc_, t);
        }
      }
    }
  }
  

  // 3rd matrice

  for (unsigned int i = 0; i < nbStates_; ++i) {
    if (ilambda[i]<=0)
    {
      for (unsigned int j = 0; j < nbStates_; ++j) 
        vresult[0][2](i,j)=0;
    }
    else
    {
      vfunc_[0]=expLam[i];
      vfunc_[1]=cosLam[i];
      vfunc_[2]=sinLam[i];
      fact_[0]=lambda[i];
      fact_[1]=ilambda[i];
      for (unsigned int j = 0; j < nbStates_; ++j) 
      {
        vfunc_[3]=expLam[j];
        fact_[2]=lambda[j];
        
        if (ilambda[j]==0)
          vresult[0][2](i, j) = computeKs_(fact_, vfunc_, t);
        else
        {
          vfunc_[4]=cosLam[j];
          vfunc_[5]=sinLam[j];
          fact_[3]=ilambda[j];
          vresult[0][2](i,j) = computeJsc_(fact_, vfunc_, t);
        }
      }
    }
  }

  // 2nd block
  // 1st matrice

  for (unsigned int i = 0; i < nbStates_; ++i) {
    vfunc_[0]=expLam[i];
    fact_[0]=lambda[i];
    if (ilambda[i]==0)
    {
      for (unsigned int j = 0; j < nbStates_; ++j) {
        if (ilambda[j]>=0)
          vresult[1][0](i, j) = 0;
        else
        {
          fact_[1]=lambda[j];
          fact_[2]=ilambda[j];
          vfunc_[1]=expLam[j];
          vfunc_[2]=cosLam[j];
          vfunc_[3]=sinLam[j];
          vresult[1][0](i,j) = computeJs_(fact_, vfunc_, t);
        }
      }
    }
    else
    {
      vfunc_[1]=cosLam[i];
      vfunc_[2]=sinLam[i];
      fact_[1]=ilambda[i];

      for (unsigned int j = 0; j < nbStates_; ++j) {
        if (ilambda[j]>=0)
          vresult[1][0](i, j) = 0;
        else
        {
          fact_[2]=lambda[j];
          fact_[3]=ilambda[j];
          vfunc_[3]=expLam[j];
          vfunc_[4]=cosLam[j];
          vfunc_[5]=sinLam[j];
          vresult[1][0](i,j) = computeJcs_(fact_, vfunc_, t);
        }
      }
    }
  }

  // 2nd matrice

  for (unsigned int i = 0; i < nbStates_; ++i) {
    if (ilambda[i]>=0)
    {
      for (unsigned int j = 0; j < nbStates_; ++j) 
        vresult[1][1](i,j)=0;
    }
    else
    {
      fact_[0]=lambda[i];
      fact_[1]=ilambda[i];

      vfunc_[0]=expLam[i];
      vfunc_[1]=cosLam[i];
      vfunc_[2]=sinLam[i];
      for (unsigned int j = 0; j < nbStates_; ++j) 
      {
        if (ilambda[j]>=0)
          vresult[1][1](i, j) = 0;
        else
        {
          fact_[2]=lambda[j];
          fact_[3]=ilambda[j];
        
          vfunc_[3]=expLam[j];
          vfunc_[4]=cosLam[j];
          vfunc_[5]=sinLam[j];
          vresult[1][1](i,j) = computeJss_(fact_, vfunc_, t);
        }
      }
    }
  }
  

  // 3rd matrice

  for (unsigned int i = 0; i < nbStates_; ++i) {
    if (ilambda[i]<=0)
    {
      for (unsigned int j = 0; j < nbStates_; ++j) 
        vresult[1][2](i,j)=0;
    }
    else
    {
      vfunc_[0]=expLam[i];
      vfunc_[1]=cosLam[i];
      vfunc_[2]=sinLam[i];
      fact_[0]=lambda[i];
      fact_[1]=ilambda[i];
      for (unsigned int j = 0; j < nbStates_; ++j) 
      {
        if (ilambda[j]>=0)
          vresult[1][2](i, j) = 0;
        else
        {
          fact_[2]=lambda[j];
          fact_[3]=ilambda[j];
        
          vfunc_[3]=expLam[j];
          vfunc_[4]=cosLam[j];
          vfunc_[5]=sinLam[j];
          vresult[1][2](i,j) = computeJss_(fact_, vfunc_, t);
        }
      }
    }
  }
  
  // 3rd block
  // 1st matrice

  for (unsigned int i = 0; i < nbStates_; ++i) {
    vfunc_[0]=expLam[i];
    fact_[0]=lambda[i];
    if (ilambda[i]==0)
    {
      for (unsigned int j = 0; j < nbStates_; ++j) {
        if (ilambda[j]<=0)
          vresult[2][0](i, j) = 0;
        else
        {
          fact_[1]=lambda[j];
          fact_[2]=ilambda[j];
          vfunc_[1]=expLam[j];
          vfunc_[2]=cosLam[j];
          vfunc_[3]=sinLam[j];
          vresult[2][0](i,j) = computeJs_(fact_, vfunc_, t);
        }
      }
    }
    else
    {
      vfunc_[1]=cosLam[i];
      vfunc_[2]=sinLam[i];
      fact_[1]=ilambda[i];

      for (unsigned int j = 0; j < nbStates_; ++j) {
        if (ilambda[j]<=0)
          vresult[2][0](i, j) = 0;
        else
        {
          fact_[2]=lambda[j];
          fact_[3]=ilambda[j];
          vfunc_[3]=expLam[j];
          vfunc_[4]=cosLam[j];
          vfunc_[5]=sinLam[j];
          vresult[2][0](i,j) = computeJcs_(fact_, vfunc_, t);
        }
      }
    }
  }

  // 2nd matrice

  for (unsigned int i = 0; i < nbStates_; ++i) {
    if (ilambda[i]>=0)
    {
      for (unsigned int j = 0; j < nbStates_; ++j) 
        vresult[2][1](i,j)=0;
    }
    else
    {
      fact_[0]=lambda[i];
      fact_[1]=ilambda[i];

      vfunc_[0]=expLam[i];
      vfunc_[1]=cosLam[i];
      vfunc_[2]=sinLam[i];
      for (unsigned int j = 0; j < nbStates_; ++j) 
      {
        if (ilambda[j]<=0)
          vresult[2][1](i, j) = 0;
        else
        {
          fact_[2]=lambda[j];
          fact_[3]=ilambda[j];
        
          vfunc_[3]=expLam[j];
          vfunc_[4]=cosLam[j];
          vfunc_[5]=sinLam[j];
          vresult[2][1](i,j) = computeJss_(fact_, vfunc_, t);
        }
      }
    }
  }
  

  // 3rd matrice

  for (unsigned int i = 0; i < nbStates_; ++i) {
    if (ilambda[i]<=0)
    {
      for (unsigned int j = 0; j < nbStates_; ++j) 
        vresult[2][2](i,j)=0;
    }
    else
    {
      vfunc_[0]=expLam[i];
      vfunc_[1]=cosLam[i];
      vfunc_[2]=sinLam[i];
      fact_[0]=lambda[i];
      fact_[1]=ilambda[i];
      for (unsigned int j = 0; j < nbStates_; ++j) 
      {
        if (ilambda[j]<=0)
          vresult[2][2](i, j) = 0;
        else
        {
          fact_[2]=lambda[j];
          fact_[3]=ilambda[j];
        
          vfunc_[3]=expLam[j];
          vfunc_[4]=cosLam[j];
          vfunc_[5]=sinLam[j];
          vresult[2][2](i,j) = computeJss_(fact_, vfunc_, t);
        }
      }
    }
  }
}

/******************************************************************************/

void DecompositionMethods::setSubstitutionModel(const SubstitutionModel* model)
{
  model_ = model;
  size_t n = model->getNumberOfStates();
  if (n != nbStates_) {
    nbStates_ = n;
    initStates_();
    initBMatrices_();
  }
}


