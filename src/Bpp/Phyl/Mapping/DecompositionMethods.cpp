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

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include <vector>
#include <typeinfo>

using namespace std;

using namespace bpp;

/******************************************************************************/

DecompositionMethods::DecompositionMethods(const SubstitutionModel* model, SubstitutionRegister* reg) :
  model_(0), // not (model) so setSubstitutionModel is effective
  nbStates_(model->getNumberOfStates()),
  nbTypes_(reg->getNumberOfSubstitutionTypes()),
  jMat_(nbStates_, nbStates_),
  jIMat_(0,0),
  rightEigenVectors_(0,0),
  rightIEigenVectors_(0,0),
  leftEigenVectors_(0,0),
  leftIEigenVectors_(0,0),
  bMatrices_(reg->getNumberOfSubstitutionTypes()),
  insideProducts_(reg->getNumberOfSubstitutionTypes()),
  insideIProducts_(0)
{
  setSubstitutionModel(model);
}				


DecompositionMethods::DecompositionMethods(const SubstitutionModel* model) :
  model_(0), // not (model) so setSubstitutionModel is effective
  nbStates_(model->getNumberOfStates()),
  nbTypes_(1),
  jMat_(nbStates_, nbStates_),
  jIMat_(0,0),
  rightEigenVectors_(0,0),
  rightIEigenVectors_(0,0),
  leftEigenVectors_(0,0),
  leftIEigenVectors_(0,0),
  bMatrices_(1),
  insideProducts_(1),
  insideIProducts_(0)
{
  setSubstitutionModel(model);
}				

void DecompositionMethods::computeProducts_()
{
  //vInv_ %*% bMatrices_[i] %*% v_;
  if (model_->isDiagonalizable())
  {
    for (size_t i = 0; i < nbTypes_; ++i) {
      RowMatrix<double> tmp(nbStates_, nbStates_);
      MatrixTools::mult(model_->getRowLeftEigenVectors(), bMatrices_[i], tmp);
      MatrixTools::mult(tmp, model_->getColumnRightEigenVectors(), insideProducts_[i]);
    }
  }
  else
  {
    for (size_t i = 0; i < nbTypes_; ++i) {
      //vInv_ %*% bMatrices_[i] %*% v_;
      RowMatrix<double> tmp(nbStates_, nbStates_), itmp(nbStates_, nbStates_);
      
      MatrixTools::mult(leftEigenVectors_, bMatrices_[i], tmp);
      MatrixTools::mult(leftIEigenVectors_, bMatrices_[i], itmp);
      MatrixTools::mult(tmp, itmp, rightEigenVectors_, rightIEigenVectors_, insideProducts_[i], insideIProducts_[i]);
    }
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

void DecompositionMethods::jFunction_(const std::vector<double>& lambda, const std::vector<double>& ilambda, double t, RowMatrix<double>& result, RowMatrix<double>& iresult) const
{
  vector<double> expLam = VectorTools::exp(lambda * t);
  vector<double> cosLam = expLam * VectorTools::cos(ilambda * t);
  vector<double> sinLam = expLam * VectorTools::sin(ilambda * t);

  for (unsigned int i = 0; i < nbStates_; ++i) {
    for (unsigned int j = 0; j < nbStates_; ++j) {
      double dd = lambda[i] - lambda[j];
      double idd = ilambda[i] - ilambda[j];
      if (dd == 0 && idd == 0) {
        result(i, j) = t * cosLam[i];
        iresult(i, j) = t * sinLam[i];
      }
      else
      {
        double es = sinLam[i] - sinLam[j];
        double ec = cosLam[i] - cosLam[j];
        double num = dd * dd + idd * idd;
        
        result(i, j) = (dd * ec + idd * es) / num;
        iresult(i, j) = (dd * es - idd * ec) / num;
      }
    }
  }
}


void DecompositionMethods::computeExpectations(RowMatrix<double>& mapping, double length) const
{
  RowMatrix<double> tmp1(nbStates_, nbStates_), tmp2(nbStates_, nbStates_);
  if (model_->isDiagonalizable())
  {
    jFunction_(model_->getEigenValues(), length, jMat_);
    
    MatrixTools::hadamardMult(jMat_, insideProducts_[0], tmp1);
    MatrixTools::mult(model_->getColumnRightEigenVectors(), tmp1, tmp2);
    MatrixTools::mult(tmp2, model_->getRowLeftEigenVectors(), mapping);
  }
  else if (model_->isNonSingular())
  {
    jFunction_(model_->getEigenValues(), model_->getIEigenValues(), length, jMat_, jIMat_);
    
    RowMatrix<double> itmp1(nbStates_, nbStates_), itmp2(nbStates_, nbStates_);
    RowMatrix<double> imat(nbStates_, nbStates_);
    
    MatrixTools::hadamardMult(jMat_, jIMat_, insideProducts_[0], insideIProducts_[0], tmp1, itmp1);
    MatrixTools::mult(rightEigenVectors_, rightIEigenVectors_, tmp1, itmp1, tmp2, itmp2);
    MatrixTools::mult(tmp2, itmp2, leftEigenVectors_, leftIEigenVectors_, mapping, imat);
  }
  else
    throw Exception("void DecompositionMethods::computeMappings : substitution mapping is not implemented for singular generators.");
}


void DecompositionMethods::computeExpectations(std::vector< RowMatrix<double> >& mappings, double length) const
{
  RowMatrix<double> tmp1(nbStates_, nbStates_), tmp2(nbStates_, nbStates_);

  if (model_->isDiagonalizable())
  {
    jFunction_(model_->getEigenValues(), length, jMat_);

    for (size_t i = 0; i < nbTypes_; ++i) {
      MatrixTools::hadamardMult(jMat_, insideProducts_[i], tmp1);
      MatrixTools::mult(model_->getColumnRightEigenVectors(), tmp1, tmp2);
      MatrixTools::mult(tmp2, model_->getRowLeftEigenVectors(), mappings[i]);
    }
  }
  else if (model_->isNonSingular())
  {
    jFunction_(model_->getEigenValues(), model_->getIEigenValues(), length, jMat_, jIMat_);
    
    RowMatrix<double> itmp1(nbStates_, nbStates_), itmp2(nbStates_, nbStates_);
    RowMatrix<double> imat(nbStates_, nbStates_);

    for (size_t i = 0; i < nbTypes_; ++i) {
      MatrixTools::hadamardMult(jMat_, jIMat_, insideProducts_[i], insideIProducts_[i], tmp1, itmp1);
      MatrixTools::mult(rightEigenVectors_, rightIEigenVectors_, tmp1, itmp1, tmp2, itmp2);
      MatrixTools::mult(tmp2, itmp2, leftEigenVectors_, leftIEigenVectors_, mappings[i], imat);
    }
  }
  else
    throw Exception("void DecompositionMethods::computeMappings : substitution mapping is not implemented for singular generators.");

} 


void DecompositionMethods::initStates_()
{
  jMat_.resize(nbStates_, nbStates_);
}

/******************************************************************************/

void DecompositionMethods::setSubstitutionModel(const SubstitutionModel* model)
{
  if (model_==model)
    return;
  
  model_ = model;
  size_t n = model->getNumberOfStates();
  if (n != nbStates_)
  {
    nbStates_ = n;
    jMat_.resize(nbStates_, nbStates_);
    initBMatrices_();
  }
  
  if (!model_->isDiagonalizable())
  {
    jIMat_.resize(nbStates_, nbStates_);
    insideIProducts_.resize(nbTypes_);
    for (size_t i = 0; i < nbTypes_; ++i) 
      insideIProducts_[i].resize(nbStates_, nbStates_);

    // sets the imaginary parts of the eigenvectors
    if (rightEigenVectors_.getNumberOfRows()!=nbStates_)
    {
      rightEigenVectors_.resize(nbStates_,nbStates_);
      leftEigenVectors_.resize(nbStates_,nbStates_);
      rightIEigenVectors_.resize(nbStates_,nbStates_);
      leftIEigenVectors_.resize(nbStates_,nbStates_);
    }

    const ColMatrix<double>& rEV=model_->getColumnRightEigenVectors();
    const RowMatrix<double>& lEV=model_->getRowLeftEigenVectors();
    const vector<double>& vi=model_->getIEigenValues();

    for (size_t i=0; i<nbStates_; i++)
    {
      if (vi[i]==0)
      {
        rightEigenVectors_.getCol(i)=rEV.col(i);
        VectorTools::fill(rightIEigenVectors_.getCol(i),0.);
        leftEigenVectors_.getRow(i)=lEV.row(i);
        VectorTools::fill(leftIEigenVectors_.getRow(i),0.);
      }
      else if (vi[i]>0)
      {
        rightEigenVectors_.getCol(i)=rEV.col(i);
        rightIEigenVectors_.getCol(i)=rEV.col(i+1);
        leftEigenVectors_.getRow(i)=lEV.row(i)/2;
        leftIEigenVectors_.getRow(i)=lEV.row(i+1)*(-1./2);
      }
      else
      {
        rightEigenVectors_.getCol(i)=rEV.col(i-1);
        rightIEigenVectors_.getCol(i)=rEV.col(i)*(-1);
        leftEigenVectors_.getRow(i)=lEV.row(i-1)/2;
        leftIEigenVectors_.getRow(i)=lEV.row(i)/2;
      }
    }
  }
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

  if (!model_->isDiagonalizable())
  {
    insideIProducts_.resize(nbTypes_);
    for (size_t i = 0; i < nbTypes_; ++i) 
      insideIProducts_[i].resize(nbStates_, nbStates_);
  }
}



