// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <typeinfo>
#include <vector>

#include "Bpp/Numeric/Matrix/MatrixTools.h"
#include "DecompositionSubstitutionCount.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

DecompositionSubstitutionCount::DecompositionSubstitutionCount(
    shared_ptr<const SubstitutionModelInterface> model,
    shared_ptr<const SubstitutionRegisterInterface> reg,
    shared_ptr<const AlphabetIndex2> weights,
    shared_ptr<const AlphabetIndex2> distances) :
  AbstractSubstitutionCount(reg),
  AbstractWeightedSubstitutionCount(weights),
  AbstractSubstitutionDistance(distances),
  DecompositionMethods(model, reg),
  counts_(reg->getNumberOfSubstitutionTypes()),
  currentLength_(0)
{
  // Check compatibility between model and substitution register:
  if (typeid(model->getAlphabet()) != typeid(reg->getAlphabet()))
    throw Exception("DecompositionSubstitutionCount (constructor): alphabets do not match between register and model.");

  initCounts_();

  initBMatrices_();
  fillBMatrices_();
  computeProducts_();
}

DecompositionSubstitutionCount::DecompositionSubstitutionCount(
    shared_ptr<const SubstitutionRegisterInterface> reg,
    shared_ptr<const AlphabetIndex2> weights,
    shared_ptr<const AlphabetIndex2> distances) :
  AbstractSubstitutionCount(reg),
  AbstractWeightedSubstitutionCount(weights),
  AbstractSubstitutionDistance(distances),
  DecompositionMethods(reg),
  counts_(reg->getNumberOfSubstitutionTypes()),
  currentLength_(0)
{
  initCounts_();
}


void DecompositionSubstitutionCount::initCounts_()
{
  counts_.resize(register_->getNumberOfSubstitutionTypes());
  // Re-initialize all count matrices according to substitution register.
  for (size_t i = 0; i < register_->getNumberOfSubstitutionTypes(); ++i)
  {
    counts_[i].resize(nbStates_, nbStates_);
  }
}

/*************************************************/

void DecompositionSubstitutionCount::fillBMatrices_()
{
  if (!model_)
    throw Exception("DecompositionSubstitutionCount::fillBMatrices_: model not defined.");

  vector<int> supportedStates = model_->getAlphabetStates();
  for (size_t j = 0; j < nbStates_; ++j)
  {
    for (size_t k = 0; k < nbStates_; ++k)
    {
      size_t i = register_->getType(j, k);
      if (i > 0 && k != j)
      {
        bMatrices_[i - 1](j, k) = model_->Qij(j, k);
      }
    }
  }

  if (distances_)
    setDistanceBMatrices_();
}

void DecompositionSubstitutionCount::setDistanceBMatrices_()
{
  if (!model_)
    throw Exception("DecompositionSubstitutionCount::setDistanceBMatrices_: model not defined.");

  vector<int> supportedStates = model_->getAlphabetStates();
  for (size_t j = 0; j < nbStates_; ++j)
  {
    for (size_t k = 0; k < nbStates_; ++k)
    {
      size_t i = register_->getType(j, k);
      if (i > 0 && k != j)
      {
        bMatrices_[i - 1](j, k) *= distances_->getIndex(supportedStates[j], supportedStates[k]);
      }
    }
  }
}

/******************************************************************************/

void DecompositionSubstitutionCount::computeCounts_(double length) const
{
  computeExpectations(counts_, length);

  // Now we must divide by pijt and account for putative weights:
  vector<int> supportedStates = model_->getAlphabetStates();
  RowMatrix<double> P = model_->getPij_t(length);
  for (size_t i = 0; i < nbTypes_; i++)
  {
    for (size_t j = 0; j < nbStates_; j++)
    {
      for (size_t k = 0; k < nbStates_; k++)
      {
        counts_[i](j, k) /= P(j, k);
        if (!std::isnormal(counts_[i](j, k)))
          counts_[i](j, k) = 0.;
        if (weights_)
          counts_[i](j, k) *= weights_->getIndex(supportedStates[j], supportedStates[k]);
      }
    }
  }
}

/******************************************************************************/

unique_ptr<Matrix<double>> DecompositionSubstitutionCount::getAllNumbersOfSubstitutions(
    double length, size_t type) const
{
  if (!model_)
    throw Exception("DecompositionSubstitutionCount::getAllNumbersOfSubstitutions: model not defined.");

  if (length < 0)
    throw Exception("DecompositionSubstitutionCount::getAllNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeCounts_(length);
    currentLength_ = length;
  }
  return make_unique<RowMatrix<double>>(counts_[type - 1]);
}

/******************************************************************************/

void DecompositionSubstitutionCount::storeAllNumbersOfSubstitutions(double length, size_t type, Eigen::MatrixXd& mat) const
{
  if (!model_)
    throw Exception("DecompositionSubstitutionCount::storeAllNumbersOfSubstitutions: model not defined.");

  if (length < 0)
    throw Exception("DecompositionSubstitutionCount::getAllNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeCounts_(length);
    currentLength_ = length;
  }

  mat.resize(Eigen::Index(nbStates_), Eigen::Index(nbStates_));

  const auto& ct = counts_[type - 1];
  for (size_t i = 0; i < nbStates_; i++)
  {
    for (size_t j = 0; j < nbStates_; j++)
    {
      mat(Eigen::Index(i), Eigen::Index(j)) = isnan(ct(i, j)) ? 0 : ct(i, j);
    }
  }
}

/******************************************************************************/

double DecompositionSubstitutionCount::getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type) const
{
  if (!model_)
    throw Exception("DecompositionSubstitutionCount::getNumberOfSubstitutions: model not defined.");

  if (length < 0)
    throw Exception("DecompositionSubstitutionCount::getNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeCounts_(length);
    currentLength_ = length;
  }
  return counts_[type - 1](initialState, finalState);
}

/******************************************************************************/

std::vector<double> DecompositionSubstitutionCount::getNumberOfSubstitutionsPerType(size_t initialState, size_t finalState, double length) const
{
  if (!model_)
    throw Exception("DecompositionSubstitutionCount::getNumberOfSubstitutionsPerType: model not defined.");

  if (length < 0)
    throw Exception("DecompositionSubstitutionCount::getNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");

  if (length != currentLength_)
  {
    computeCounts_(length);
    currentLength_ = length;
  }
  std::vector<double> v(getNumberOfSubstitutionTypes());
  for (size_t t = 0; t < getNumberOfSubstitutionTypes(); ++t)
  {
    v[t] = counts_[t](initialState, finalState);
  }
  return v;
}

/******************************************************************************/

void DecompositionSubstitutionCount::setSubstitutionModel(
    shared_ptr<const SubstitutionModelInterface> model)
{
  // Check compatibility between model and substitution register:
  if (typeid(model->getAlphabet()) != typeid(register_->getAlphabet()))
    throw Exception("DecompositionMethods::setSubstitutionModel: alphabets do not match between register and model.");

  DecompositionMethods::setSubstitutionModel(model);

  initCounts_();

  fillBMatrices_();
  computeProducts_();

  // Recompute counts:
  if (currentLength_ > 0)
    computeCounts_(currentLength_);
}

/******************************************************************************/

void DecompositionSubstitutionCount::substitutionRegisterHasChanged()
{
  if (!model_)
    return;

  // Check compatibility between model and substitution register:
  if (model_->getAlphabet()->getAlphabetType() != register_->getAlphabet()->getAlphabetType())
    throw Exception("DecompositionMethods::substitutionRegisterHasChanged: alphabets do not match between register and model.");

  initBMatrices_();
  initStates_();

  initCounts_();

  fillBMatrices_();
  computeProducts_();

  // Recompute counts:
  if (currentLength_ > 0)
    computeCounts_(currentLength_);
}

/******************************************************************************/

void DecompositionSubstitutionCount::weightsHaveChanged()
{
  if (typeid(weights_->getAlphabet()) != typeid(register_->getAlphabet()))
    throw Exception("DecompositionSubstitutionCount::weightsHaveChanged. Incorrect alphabet type.");

  if (!model_)
    return;

  // Recompute counts:
  if (currentLength_ > 0)
    computeCounts_(currentLength_);
}


void DecompositionSubstitutionCount::distancesHaveChanged()
{
  if (distances_->getAlphabet()->getAlphabetType() != register_->getAlphabet()->getAlphabetType())
    throw Exception("DecompositionSubstitutionCount::distancesHaveChanged. Incorrect alphabet type.");

  if (!model_)
    return;

  // Recompute counts:
  setDistanceBMatrices_();
  computeProducts_();

  if (currentLength_ > 0)
    computeCounts_(currentLength_);
}

/******************************************************************************/
