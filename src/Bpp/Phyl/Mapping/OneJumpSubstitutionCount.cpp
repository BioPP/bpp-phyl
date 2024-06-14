// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "OneJumpSubstitutionCount.h"

using namespace bpp;
using namespace std;

unique_ptr<Matrix<double>> OneJumpSubstitutionCount::getAllNumbersOfSubstitutions(double length, size_t type) const
{
  if (!model_)
    throw Exception("OneJumpSubstitutionCount::getAllNumberOfSubstitutions: model not defined.");

  tmp_ = model_->getPij_t(length);
  size_t n = model_->getNumberOfStates();
  auto probs = make_unique<LinearMatrix<double>>(n, n);
  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = 0; j < n; ++j)
    {
      (*probs)(i, j) = (i == j ? 1. - tmp_(i, j) : 1.);
    }
  }
  return probs;
}

void OneJumpSubstitutionCount::storeAllNumbersOfSubstitutions(double length, size_t type, Eigen::MatrixXd& mat) const
{
  if (!model_)
    throw Exception("OneJumpSubstitutionCount::storeNumberOfSubstitutions: model not defined.");

  tmp_ = model_->getPij_t(length);
  auto n = Eigen::Index(model_->getNumberOfStates());

  mat.resize(n, n);
  for (auto i = 0; i < n; i++)
  {
    for (auto j = 0; j < n; j++)
    {
      mat(i, j) = (i == j ? 1. - tmp_(size_t(i), size_t(j)) : 1.);
    }
  }
}
