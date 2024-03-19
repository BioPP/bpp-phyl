// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "NaiveSubstitutionCount.h"

using namespace bpp;
using namespace std;

unique_ptr<Matrix<double>> NaiveSubstitutionCount::getAllNumbersOfSubstitutions(double length, size_t type) const
{
  size_t n = supportedChars_.size();
  auto mat = make_unique<RowMatrix<double>>(n, n);
  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = 0; j < n; ++j)
    {
      (*mat)(i, j) = (register_->getType(i, j) == type ? (weights_ ? weights_->getIndex(supportedChars_[i], supportedChars_[j]) : 1.) : 0.);
    }
  }
  return mat;
}

void NaiveSubstitutionCount::storeAllNumbersOfSubstitutions(double length, size_t type, Eigen::MatrixXd& mat) const
{
  auto n = Eigen::Index(supportedChars_.size());
  mat.resize(n, n);

  for (auto i = 0; i < n; ++i)
  {
    for (auto j = 0; j < n; ++j)
    {
      mat(i, j) = (register_->getType(size_t(i), size_t(j)) == type ? (weights_ ? weights_->getIndex(supportedChars_[size_t(i)], supportedChars_[size_t(j)]) : 1.) : 0.);
    }
  }
}


LabelSubstitutionCount::LabelSubstitutionCount(
    shared_ptr<const SubstitutionModelInterface> model) :
  AbstractSubstitutionCount(
    make_shared<TotalSubstitutionRegister>(model->getStateMap())),
  label_(model->getNumberOfStates(), model->getNumberOfStates()),
  supportedChars_(model->getAlphabetStates())
{
  size_t n = supportedChars_.size();
  double count = 0;
  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = 0; j < n; ++j)
    {
      if (i == j)
        label_(i, j) = 0;
      else
        label_(i, j) = ++count;
    }
  }
}

LabelSubstitutionCount::LabelSubstitutionCount(
    shared_ptr<const StateMapInterface> stateMap) :
  AbstractSubstitutionCount(
    make_shared<TotalSubstitutionRegister>(stateMap)),
  label_(stateMap->getNumberOfModelStates(), stateMap->getNumberOfModelStates()),
  supportedChars_(stateMap->getAlphabetStates())
{
  size_t n = supportedChars_.size();
  double count = 0;
  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = 0; j < n; ++j)
    {
      if (i == j)
        label_(i, j) = 0;
      else
        label_(i, j) = ++count;
    }
  }
}
