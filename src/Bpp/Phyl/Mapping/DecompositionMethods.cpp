// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <typeinfo>
#include <vector>

#include "DecompositionMethods.h"

using namespace std;

using namespace bpp;

/******************************************************************************/

DecompositionMethods::DecompositionMethods(
    std::shared_ptr<const SubstitutionModelInterface> model,
    std::shared_ptr<const SubstitutionRegisterInterface> reg) :
  model_(model),
  nbStates_(model->getNumberOfStates()),
  nbTypes_(reg->getNumberOfSubstitutionTypes()),
  jMat_(nbStates_, nbStates_),
  jIMat_(0, 0),
  rightEigenVectors_(0, 0),
  rightIEigenVectors_(0, 0),
  leftEigenVectors_(0, 0),
  leftIEigenVectors_(0, 0),
  bMatrices_(reg->getNumberOfSubstitutionTypes()),
  insideProducts_(reg->getNumberOfSubstitutionTypes()),
  insideIProducts_(0)
{
  initBMatrices_();
  setSubstitutionModel(model);
}

DecompositionMethods::DecompositionMethods(
    std::shared_ptr<const SubstitutionRegisterInterface> reg) :
  model_(0),
  nbStates_(reg->stateMap().getNumberOfModelStates()),
  nbTypes_(reg->getNumberOfSubstitutionTypes()),
  jMat_(nbStates_, nbStates_),
  jIMat_(0, 0),
  rightEigenVectors_(0, 0),
  rightIEigenVectors_(0, 0),
  leftEigenVectors_(0, 0),
  leftIEigenVectors_(0, 0),
  bMatrices_(reg->getNumberOfSubstitutionTypes()),
  insideProducts_(reg->getNumberOfSubstitutionTypes()),
  insideIProducts_(0)
{
  initBMatrices_();
}

DecompositionMethods::DecompositionMethods(const StateMapInterface& stateMap) :
  model_(0),
  nbStates_(stateMap.getNumberOfModelStates()),
  nbTypes_(1),
  jMat_(nbStates_, nbStates_),
  jIMat_(0, 0),
  rightEigenVectors_(0, 0),
  rightIEigenVectors_(0, 0),
  leftEigenVectors_(0, 0),
  leftIEigenVectors_(0, 0),
  bMatrices_(1),
  insideProducts_(1),
  insideIProducts_(0)
{
  initBMatrices_();
}


DecompositionMethods::DecompositionMethods(
    std::shared_ptr<const SubstitutionModelInterface> model) :
  model_(model),
  nbStates_(model->getNumberOfStates()),
  nbTypes_(1),
  jMat_(nbStates_, nbStates_),
  jIMat_(0, 0),
  rightEigenVectors_(0, 0),
  rightIEigenVectors_(0, 0),
  leftEigenVectors_(0, 0),
  leftIEigenVectors_(0, 0),
  bMatrices_(1),
  insideProducts_(1),
  insideIProducts_(0)
{
  setSubstitutionModel(model);
}

void DecompositionMethods::computeProducts_()
{
  // vInv_ %*% bMatrices_[i] %*% v_;
  if (model_->isDiagonalizable())
  {
    for (size_t i = 0; i < nbTypes_; ++i)
    {
      RowMatrix<double> tmp(nbStates_, nbStates_);
      MatrixTools::mult(model_->getRowLeftEigenVectors(), bMatrices_[i], tmp);
      MatrixTools::mult(tmp, model_->getColumnRightEigenVectors(), insideProducts_[i]);
    }
  }
  else
  {
    for (size_t i = 0; i < nbTypes_; ++i)
    {
      // vInv_ %*% bMatrices_[i] %*% v_;
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
  for (unsigned int i = 0; i < nbStates_; ++i)
  {
    for (unsigned int j = 0; j < nbStates_; ++j)
    {
      double dd = lambda[i] - lambda[j];
      if (abs(dd) < NumConstants::TINY())
        result(i, j) = t * expLam[i];
      else
      {
        result(i, j) = (expLam[i] - expLam[j]) / dd;
      }
    }
  }
}

void DecompositionMethods::jFunction_(const std::vector<double>& lambda, const std::vector<double>& ilambda, double t, RowMatrix<double>& result, RowMatrix<double>& iresult) const
{
  vector<double> expLam = VectorTools::exp(lambda * t);
  auto x = ilambda * t;
  vector<double> cosLam = expLam * VectorTools::cos(x);
  vector<double> sinLam = expLam * VectorTools::sin(x);

  for (unsigned int i = 0; i < nbStates_; ++i)
  {
    for (unsigned int j = 0; j < nbStates_; ++j)
    {
      double dd = lambda[i] - lambda[j];
      double idd = ilambda[i] - ilambda[j];
      if ((abs(dd) < NumConstants::TINY()) && (abs(idd) < NumConstants::TINY()))
      {
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
  if (!model_)
    throw Exception("DecompositionMethods::computeExpectations: model not defined.");

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


void DecompositionMethods::computeExpectations(std::vector< RowMatrix<double>>& mappings, double length) const
{
  if (!model_)
    throw Exception("DecompositionMethods::computeExpectations: model not defined.");

  RowMatrix<double> tmp1(nbStates_, nbStates_), tmp2(nbStates_, nbStates_);

  if (model_->isDiagonalizable())
  {
    jFunction_(model_->getEigenValues(), length, jMat_);

    for (size_t i = 0; i < nbTypes_; ++i)
    {
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

    for (size_t i = 0; i < nbTypes_; ++i)
    {
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
  initBMatrices_();
}

/******************************************************************************/

void DecompositionMethods::setSubstitutionModel(
    std::shared_ptr<const SubstitutionModelInterface> model)
{
  model_ = model;
  if (!model)
    return;

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
    {
      insideIProducts_[i].resize(nbStates_, nbStates_);
    }

    // sets the imaginary parts of the eigenvectors
    if (rightEigenVectors_.getNumberOfRows() != nbStates_)
    {
      rightEigenVectors_.resize(nbStates_, nbStates_);
      leftEigenVectors_.resize(nbStates_, nbStates_);
      rightIEigenVectors_.resize(nbStates_, nbStates_);
      leftIEigenVectors_.resize(nbStates_, nbStates_);
    }

    const ColMatrix<double>& rEV = model_->getColumnRightEigenVectors();
    const RowMatrix<double>& lEV = model_->getRowLeftEigenVectors();
    const vector<double>& vi = model_->getIEigenValues();

    for (size_t i = 0; i < nbStates_; i++)
    {
      if (vi[i] == 0)
      {
        rightEigenVectors_.getCol(i) = rEV.col(i);
        VectorTools::fill(rightIEigenVectors_.getCol(i), 0.);
        leftEigenVectors_.getRow(i) = lEV.row(i);
        VectorTools::fill(leftIEigenVectors_.getRow(i), 0.);
      }
      else if (vi[i] > 0)
      {
        rightEigenVectors_.getCol(i) = rEV.col(i);
        rightIEigenVectors_.getCol(i) = rEV.col(i + 1);
        leftEigenVectors_.getRow(i) = lEV.row(i) / 2;
        leftIEigenVectors_.getRow(i) = lEV.row(i + 1) * (-1. / 2);
      }
      else
      {
        rightEigenVectors_.getCol(i) = rEV.col(i - 1);
        rightIEigenVectors_.getCol(i) = rEV.col(i) * (-1);
        leftEigenVectors_.getRow(i) = lEV.row(i - 1) / 2;
        leftIEigenVectors_.getRow(i) = lEV.row(i) / 2;
      }
    }
  }
}

void DecompositionMethods::initBMatrices_()
{
  // Re-initialize all B matrices according to substitution register.
  bMatrices_.resize(nbTypes_);
  insideProducts_.resize(nbTypes_);
  for (size_t i = 0; i < nbTypes_; ++i)
  {
    bMatrices_[i].resize(nbStates_, nbStates_);
    insideProducts_[i].resize(nbStates_, nbStates_);
  }
}
