// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <typeinfo>
#include <vector>

#include "Bpp/Numeric/Matrix/MatrixTools.h"
#include "DecompositionReward.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

DecompositionReward::DecompositionReward(
    shared_ptr<const SubstitutionModelInterface> model,
    shared_ptr<const AlphabetIndex1> alphIndex) :
  AbstractReward(alphIndex),
  DecompositionMethods(model),
  rewards_(nbStates_, nbStates_),
  currentLength_(-1.)
{
  // Check compatiblity between model and alphabet Index:
  if (typeid(model->getAlphabet()) != typeid(alphIndex_->getAlphabet()))
    throw Exception("DecompositionReward (constructor): alphabets do not match between alphabet index and model.");

  // Initialize the B matrice. This is done once for all,
  // unless the number of states changes:

  initBMatrices_();
  initRewards_();

  fillBMatrice_();
  computeProducts_();
}

DecompositionReward::DecompositionReward(
    const StateMapInterface& stateMap,
    shared_ptr<const AlphabetIndex1> alphIndex) :
  AbstractReward(alphIndex),
  DecompositionMethods(stateMap),
  rewards_(nbStates_, nbStates_),
  currentLength_(-1.)
{}

/******************************************************************************/

void DecompositionReward::initRewards_()
{
  rewards_.resize(nbStates_, nbStates_);
}

/******************************************************************************/

void DecompositionReward::fillBMatrice_()
{
  vector<int> supportedStates = model_->getAlphabetStates();
  for (size_t j = 0; j < nbStates_; ++j)
  {
    bMatrices_[0](j, j) = getAlphabetIndex()->getIndex(supportedStates[j]);
  }
}

/******************************************************************************/

void DecompositionReward::computeRewards_(double length) const
{
  computeExpectations(rewards_, length);

  // Now we must divide by pijt:
  RowMatrix<double> P = model_->getPij_t(length);
  for (size_t j = 0; j < nbStates_; j++)
  {
    for (size_t k = 0; k < nbStates_; k++)
    {
      rewards_(j, k) /= P(j, k);
      if (std::isnan(rewards_(j, k)) || std::isnan(-rewards_(j, k)) || std::isinf(rewards_(j, k)))
        rewards_(j, k) = 0.;
    }
  }
}

/******************************************************************************/

Matrix<double>* DecompositionReward::getAllRewards(double length) const
{
  if (!model_)
    throw Exception("DecompositionReward::getAllRewards: model not defined.");

  if (length < 0)
    throw Exception("DecompositionReward::getAllRewards. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeRewards_(length);
    currentLength_ = length;
  }
  return new RowMatrix<double>(rewards_);
}

/******************************************************************************/

void DecompositionReward::storeAllRewards(double length, Eigen::MatrixXd& mat) const
{
  if (!model_)
    throw Exception("DecompositionReward::storeAllRewards: model not defined.");

  if (length < 0)
    throw Exception("DecompositionReward::storeAllRewards. Negative branch length: " + TextTools::toString(length) + ".");

  if (length != currentLength_)
  {
    computeRewards_(length);
    currentLength_ = length;
  }

  mat.resize(Eigen::Index(nbStates_), Eigen::Index(nbStates_));

  for (size_t j = 0; j < nbStates_; j++)
  {
    for (size_t k = 0; k < nbStates_; k++)
    {
      mat(Eigen::Index(j), Eigen::Index(k)) = rewards_(j, k);
    }
  }
}

/******************************************************************************/

double DecompositionReward::getReward(size_t initialState, size_t finalState, double length) const
{
  if (length < 0)
    throw Exception("DecompositionReward::getRewards. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeRewards_(length);
    currentLength_ = length;
  }
  return rewards_(initialState, finalState);
}

/******************************************************************************/

void DecompositionReward::setSubstitutionModel(
    shared_ptr<const SubstitutionModelInterface> model)
{
  DecompositionMethods::setSubstitutionModel(model);

  if (!model)
    return;

  // Check compatiblity between model and substitution register:
  if (typeid(model->getAlphabet()) != typeid(alphIndex_->getAlphabet()))
    throw Exception("DecompositionReward::setSubstitutionModel: alphabets do not match between alphabet index and model.");


  initRewards_();

  fillBMatrice_();
  computeProducts_();

  // Recompute rewards:

  computeRewards_(currentLength_);
}

/******************************************************************************/

void DecompositionReward::alphabetIndexHasChanged()
{
  if (!model_)
    return;

  // Check compatiblity between model and substitution register:
  if (typeid(model_->getAlphabet()) != typeid(alphIndex_->getAlphabet()))
    throw Exception("DecompositionReward::AlphabetIndexHasChanged: alphabets do not match between alphbaet index and model.");

  fillBMatrice_();
  computeProducts_();

  // Recompute rewards:
  if (currentLength_ > 0)
    computeRewards_(currentLength_);
}
