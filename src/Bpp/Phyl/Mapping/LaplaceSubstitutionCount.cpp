// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "LaplaceSubstitutionCount.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

void LaplaceSubstitutionCount::computeCounts(double length) const
{
  RowMatrix<double> Q = model_->generator();
  double rate = model_->getRate();
  MatrixTools::scale(Q, rate);

  // L is the diagonal matrix with all substitution rates.
  size_t s = Q.getNumberOfRows();
  RowMatrix<double> QL(s, s);
  for (size_t i = 0; i < s; i++)
  {
    for (size_t j = 0; j < s; j++)
    {
      QL(i, j) = ((i == j) ? 0. : Q(i, j));
    }
  }

  MatrixTools::fill(m_, 0.);
  RowMatrix<double> M2(s, s);
  RowMatrix<double> M3(s, s);
  RowMatrix<double> M4(s, s);
  RowMatrix<double> M5(s, s);
  for (size_t n = 1; n < cutOff_; ++n)
  {
    MatrixTools::fill(M2, 0.);
    for (size_t p = 0; p < n; ++p)
    {
      MatrixTools::pow(Q, p, M3);         // Q^p -> M5
      MatrixTools::mult(M3, QL, M4);      // Q^p . QL -> M4
      MatrixTools::pow(Q, n - p - 1, M3); // Q^(n-p-1) -> M3
      MatrixTools::mult(M4, M3, M5);      // Q^p . QL . Q^(n-p-1) -> M5
      MatrixTools::add(M2, M5);
    }
    MatrixTools::scale(M2, pow(length, static_cast<double>(n)) / static_cast<double>(NumTools::fact(n)));
    MatrixTools::add(m_, M2);
  }

  // Now we must divide by pijt:
  RowMatrix<double> P = model_->getPij_t(length);
  for (size_t i = 0; i < s; i++)
  {
    for (size_t j = 0; j < s; j++)
    {
      m_(i, j) /= P(i, j);
    }
  }
}

/******************************************************************************/

double LaplaceSubstitutionCount::getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type) const
{
  if (!model_)
    throw Exception("LaplaceSubstitutionCount::getNumberOfSubstitutions: model not defined.");

  if (length == currentLength_)
    return m_(initialState, finalState);
  if (length < 0.000001)
    return initialState == finalState ? 0. : 1.; // Limit case!
  // Else we need to recompute M:
  computeCounts(length);

  currentLength_ = length;
  return m_(initialState, finalState);
}

/******************************************************************************/

unique_ptr<Matrix<double>> LaplaceSubstitutionCount::getAllNumbersOfSubstitutions(double length, size_t type) const
{
  if (!model_)
    throw Exception("LaplaceSubstitutionCount::getAllNumbersOfSubstitutions: model not defined.");

  if (length == currentLength_)
    return make_unique<RowMatrix<double>>(m_);
  if (length < 0.000001) // Limit case!
  {
    size_t s = model_->getAlphabet()->getSize();
    for (size_t i = 0; i < s; i++)
    {
      for (size_t j = 0; j < s; j++)
      {
        m_(i, j) = i == j ? 0. : 1.;
      }
    }
  }
  else
  {
    // Else we need to recompute M:
    computeCounts(length);
  }

  currentLength_ = length;

  return make_unique<RowMatrix<double>>(m_);
}

/******************************************************************************/

void LaplaceSubstitutionCount::storeAllNumbersOfSubstitutions(double length, size_t type, Eigen::MatrixXd& mat) const
{
  if (!model_)
    throw Exception("LaplaceSubstitutionCount::storeAllNumbersOfSubstitutions: model not defined.");

  auto s = Eigen::Index(model_->getAlphabet()->getSize());
  if (length == currentLength_)
    mat = Eigen::MatrixXd::Zero(s, s);

  if (length < 0.000001) // Limit case!
  {
    for (auto i = 0; i < s; i++)
    {
      for (auto j = 0; j < s; j++)
      {
        mat(i, j) = i == j ? 0. : 1.;
      }
    }
  }
  else
  {
    // Else we need to recompute M:
    computeCounts(length);
  }

  currentLength_ = length;

  mat.resize(s, s);

  for (auto i = 0; i < s; i++)
  {
    for (auto j = 0; j < s; j++)
    {
      mat(i, j) = std::isnan(m_(size_t(i), size_t(j))) ? 0 : m_(size_t(i), size_t(j));
    }
  }
}

/******************************************************************************/

void LaplaceSubstitutionCount::setSubstitutionModel(
    shared_ptr<const SubstitutionModelInterface> model)
{
  model_ = model;
  if (!model)
    return;

  size_t n = model->alphabet().getSize();
  m_.resize(n, n);
  // Recompute counts:
  if (currentLength_ > 0)
    computeCounts(currentLength_);
}

/******************************************************************************/
