// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <vector>

#include "Bpp/Numeric/Matrix/MatrixTools.h"
#include "Bpp/Numeric/NumTools.h"
#include "UniformizationSubstitutionCount.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

UniformizationSubstitutionCount::UniformizationSubstitutionCount(
    shared_ptr<const SubstitutionModelInterface> model,
    shared_ptr<const SubstitutionRegisterInterface> reg,
    shared_ptr<const AlphabetIndex2> weights,
    shared_ptr<const AlphabetIndex2> distances) :
  AbstractSubstitutionCount(reg),
  AbstractWeightedSubstitutionCount(weights),
  AbstractSubstitutionDistance(distances),
  model_(model),
  nbStates_(model->getNumberOfStates()),
  bMatrices_(reg->getNumberOfSubstitutionTypes()),
  power_(),
  s_(reg->getNumberOfSubstitutionTypes()),
  miu_(0),
  counts_(reg->getNumberOfSubstitutionTypes()),
  currentLength_(0)
{
  // Check compatibility between model and substitution register:
  if (model->getAlphabet()->getAlphabetType() != reg->getAlphabet()->getAlphabetType())
    throw Exception("UniformizationSubstitutionCount (constructor): alphabets do not match between register and model.");

  // Initialize all B matrices according to substitution register. This is done once for all,
  // unless the number of states changes:
  initBMatrices_();
  fillBMatrices_();

  for (unsigned int i = 0; i < nbStates_; ++i)
  {
    double diagQ = abs(model_->Qij(i, i));
    if (diagQ > miu_)
      miu_ = diagQ;
  }

  if (miu_ > 10000)
    throw Exception("UniformizationSubstitutionCount::UniformizationSubstitutionCount The maximum diagonal values of generator is above 10000. Abort, chose another mapping method");
}

UniformizationSubstitutionCount::UniformizationSubstitutionCount(
    const StateMapInterface& stateMap,
    shared_ptr<const SubstitutionRegisterInterface> reg,
    shared_ptr<const AlphabetIndex2> weights,
    shared_ptr<const AlphabetIndex2> distances) :
  AbstractSubstitutionCount(reg),
  AbstractWeightedSubstitutionCount(weights),
  AbstractSubstitutionDistance(distances),
  model_(0),
  nbStates_(stateMap.getNumberOfModelStates()),
  bMatrices_(reg->getNumberOfSubstitutionTypes()),
  power_(),
  s_(reg->getNumberOfSubstitutionTypes()),
  miu_(0),
  counts_(reg->getNumberOfSubstitutionTypes()),
  currentLength_(0)
{}

/******************************************************************************/

void UniformizationSubstitutionCount::resetBMatrices_()
{
  size_t nbTypes = register_->getNumberOfSubstitutionTypes();
  bMatrices_.resize(nbTypes);
  counts_.resize(nbTypes);
  s_.resize(nbTypes);
}


void UniformizationSubstitutionCount::initBMatrices_()
{
  // Re-initialize all B matrices according to substitution register.
  for (size_t i = 0; i < register_->getNumberOfSubstitutionTypes(); ++i)
  {
    bMatrices_[i].resize(nbStates_, nbStates_);
    counts_[i].resize(nbStates_, nbStates_);
  }
}

void UniformizationSubstitutionCount::fillBMatrices_()
{
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

void UniformizationSubstitutionCount::setDistanceBMatrices_()
{
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

void UniformizationSubstitutionCount::computeCounts_(double length) const
{
  double lam = miu_ * length;
  RowMatrix<double> I;
  MatrixTools::getId(nbStates_, I);
  RowMatrix<double> R(model_->generator());
  MatrixTools::scale(R, model_->getRate());

  MatrixTools::scale(R, 1. / miu_);
  MatrixTools::add(R, I);

  // compute the stopping point
  // use the tail of Poisson distribution
  // can be approximated by 4 + 6 * sqrt(lam) + lam
  size_t nMax = static_cast<size_t>(ceil(4 + 6 * sqrt(lam) + lam));

  // compute the powers of R
  power_.resize(nMax + 1);
  power_[0] = I;
  for (size_t i = 1; i < nMax + 1; ++i)
  {
    MatrixTools::mult(power_[i - 1], R, power_[i]);
  }

  for (size_t i = 0; i < register_->getNumberOfSubstitutionTypes(); ++i)
  {
    s_[i].resize(nMax + 1);
    MatrixTools::mult(bMatrices_[i], power_[0], s_[i][0]);
    RowMatrix<double> tmp(nbStates_, nbStates_);
    for (size_t l = 1; l < nMax + 1; ++l)
    {
      MatrixTools::mult(R, s_[i][l - 1], s_[i][l]);
      MatrixTools::mult(bMatrices_[i], power_[l], tmp);
      MatrixTools::add(s_[i][l], tmp);
    }
    MatrixTools::fill(counts_[i], 0);
    for (size_t l = 0; l < nMax + 1; ++l)
    {
      tmp = s_[i][l];
      // double f = (pow(lam, static_cast<double>(l + 1)) * exp(-lam) / static_cast<double>(NumTools::fact(l + 1))) / miu_;
      double logF = static_cast<double>(l + 1) * log(lam) - lam - log(miu_) - NumTools::logFact(static_cast<double>(l + 1));
      MatrixTools::scale(tmp, exp(logF));
      MatrixTools::add(counts_[i], tmp);
    }
  }

  // Now we must divide by pijt and account for putative weights:
  vector<int> supportedStates = model_->getAlphabetStates();
  RowMatrix<double> P = model_->getPij_t(length);
  for (size_t i = 0; i < register_->getNumberOfSubstitutionTypes(); i++)
  {
    for (size_t j = 0; j < nbStates_; j++)
    {
      for (size_t k = 0; k < nbStates_; k++)
      {
        counts_[i](j, k) /= P(j, k);
        if (std::isinf(counts_[i](j, k)) || std::isnan(counts_[i](j, k)) || (!weights_ && counts_[i](j, k) < 0.))
          counts_[i](j, k) = 0;
        if (weights_)
          counts_[i](j, k) *= weights_->getIndex(supportedStates[j], supportedStates[k]);
      }
    }
  }
}

/******************************************************************************/

unique_ptr<Matrix<double>> UniformizationSubstitutionCount::getAllNumbersOfSubstitutions(double length, size_t type) const
{
  if (!model_)
    throw Exception("UniformizationSubstitutionCount::getAllNumbersOfSubstitutions: model not defined.");

  if (length < 0)
    throw Exception("UniformizationSubstitutionCount::getAllNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeCounts_(length);
    currentLength_ = length;
  }
  return make_unique<RowMatrix<double>>(counts_[type - 1]);
}

/******************************************************************************/

void UniformizationSubstitutionCount::storeAllNumbersOfSubstitutions(double length, size_t type, Eigen::MatrixXd& mat) const
{
  if (!model_)
    throw Exception("UniformizationSubstitutionCount::storeAllNumbersOfSubstitutions: model not defined.");

  if (length < 0)
    throw Exception("UniformizationSubstitutionCount::storeAllNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");
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

double UniformizationSubstitutionCount::getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type) const
{
  if (!model_)
    throw Exception("UniformizationSubstitutionCount::getNumberOfSubstitutions: model not defined.");

  if (length < 0)
    throw Exception("UniformizationSubstitutionCount::getNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeCounts_(length);
    currentLength_ = length;
  }
  return counts_[type - 1](initialState, finalState);
}

/******************************************************************************/

std::vector<double> UniformizationSubstitutionCount::getNumberOfSubstitutionsPerType(size_t initialState, size_t finalState, double length) const
{
  if (!model_)
    throw Exception("UniformizationSubstitutionCount::getNumberOfSubstitutionsPerTye: model not defined.");

  if (length < 0)
    throw Exception("UniformizationSubstitutionCount::getNumbersOfSubstitutions. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeCounts_(length);
    currentLength_ = length;
  }
  std::vector<double> v(getNumberOfSubstitutionTypes());
  for (unsigned int t = 0; t < getNumberOfSubstitutionTypes(); ++t)
  {
    v[t] = counts_[t](initialState, finalState);
  }
  return v;
}

/******************************************************************************/

void UniformizationSubstitutionCount::setSubstitutionModel(
    shared_ptr<const SubstitutionModelInterface> model)
{
  model_ = model;

  if (!model_)
    return;

  // Check compatibility between model and substitution register:
  if (model->alphabet().getAlphabetType() != register_->alphabet().getAlphabetType())
    throw Exception("UniformizationSubstitutionCount::setSubstitutionModel: alphabets do not match between register and model.");

  size_t n = model->getNumberOfStates();
  if (n != nbStates_)
  {
    nbStates_ = n;
    // Re-initialize all B matrices according to substitution register.
    initBMatrices_();
  }
  fillBMatrices_();

  miu_ = 0;
  for (size_t i = 0; i < nbStates_; ++i)
  {
    double diagQ = abs(model_->Qij(i, i));
    if (diagQ > miu_)
      miu_ = diagQ;
  }

  if (miu_ > 10000)
    throw Exception("UniformizationSubstitutionCount::setSubstitutionModel(). The maximum diagonal values of generator is above 10000. Abort, chose another mapping method.");

  // Recompute counts:
  if (currentLength_ > 0)
    computeCounts_(currentLength_);
}

/******************************************************************************/

void UniformizationSubstitutionCount::substitutionRegisterHasChanged()
{
  if (!model_)
    return;

  // Check compatibility between model and substitution register:
  if (model_->getAlphabet()->getAlphabetType() != register_->getAlphabet()->getAlphabetType())
    throw Exception("UniformizationSubstitutionCount::substitutionRegisterHasChanged: alphabets do not match between register and model.");

  resetBMatrices_();
  initBMatrices_();
  fillBMatrices_();

  // Recompute counts:
  if (currentLength_ > 0)
    computeCounts_(currentLength_);
}

/******************************************************************************/

void UniformizationSubstitutionCount::weightsHaveChanged()
{
  if (!model_)
    return;

  if (weights_->getAlphabet()->getAlphabetType() != register_->getAlphabet()->getAlphabetType())
    throw Exception("UniformizationSubstitutionCount::weightsHaveChanged. Incorrect alphabet type.");

  // Recompute counts:
  if (currentLength_ > 0)
    computeCounts_(currentLength_);
}

void UniformizationSubstitutionCount::distancesHaveChanged()
{
  if (!model_)
    return;

  if (distances_->getAlphabet()->getAlphabetType() != register_->getAlphabet()->getAlphabetType())
    throw Exception("UniformizationSubstitutionCount::distancesHaveChanged. Incorrect alphabet type.");

  // Recompute counts:
  setDistanceBMatrices_();

  if (currentLength_ > 0)
    computeCounts_(currentLength_);
}

/******************************************************************************/
