// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "RE08.h"

using namespace bpp;

#include <cmath>

using namespace std;

/******************************************************************************/

RE08::RE08(
    unique_ptr<ReversibleSubstitutionModelInterface> simpleModel,
    double lambda,
    double mu) :
  AbstractParameterAliasable("RE08."),
  AbstractReversibleSubstitutionModel(simpleModel->getAlphabet(), make_shared<CanonicalStateMap>(simpleModel->stateMap(), true), "RE08."),
  simpleModel_(std::move(simpleModel)),
  simpleGenerator_(),
  simpleExchangeabilities_(),
  exp_(), p_(), lambda_(lambda), mu_(mu),
  nestedPrefix_("model_" + simpleModel->getNamespace())
{
  addParameter_(new Parameter("RE08.lambda", lambda, Parameter::R_PLUS));
  addParameter_(new Parameter("RE08.mu", mu, Parameter::R_PLUS));
  simpleModel_->setNamespace("RE08." + nestedPrefix_);
  addParameters_(simpleModel_->getParameters());
  // We need to overrired this from the AbstractSubstitutionModel constructor,
  // since the number of states in the model is no longer equal to the size of the alphabet.
  size_ = simpleModel_->getNumberOfStates() + 1;
  generator_.resize(size_, size_);
  exchangeability_.resize(size_, size_);
  freq_.resize(size_);
  eigenValues_.resize(size_);
  leftEigenVectors_.resize(size_, size_);
  rightEigenVectors_.resize(size_, size_);
  p_.resize(size_, size_);
  updateMatrices_();
}

/******************************************************************************/

void RE08::updateMatrices_()
{
  double f = (lambda_ == 0 && mu_ == 0) ? 1 : lambda_ / (lambda_ + mu_);

  // Frequencies:
  for (size_t i = 0; i < size_ - 1; ++i)
  {
    freq_[i] = simpleModel_->freq(i) * f;
  }

  freq_[size_ - 1] = (1. - f);

  simpleGenerator_ = simpleModel_->generator();
  simpleExchangeabilities_ = simpleModel_->exchangeabilityMatrix();

  // Generator and exchangeabilities:
  for (size_t i = 0; i < size_ - 1; ++i)
  {
    for (size_t j = 0; j < size_ - 1; ++j)
    {
      generator_(i, j) = simpleGenerator_(i, j);
      exchangeability_(i, j) = simpleExchangeabilities_(i, j) / f;
      if (i == j)
      {
        generator_(i, j) -= mu_;
        exchangeability_(i, j) -= (mu_ / f) / simpleModel_->freq(i);
      }
    }
    generator_(i, size_ - 1) = mu_;
    generator_(size_ - 1, i) = lambda_ * simpleModel_->freq(i);
    exchangeability_(i, size_ - 1) = lambda_ + mu_;
    exchangeability_(size_ - 1, i) = lambda_ + mu_;
  }
  generator_(size_ - 1, size_ - 1) = -lambda_;
  exchangeability_(size_ - 1, size_ - 1) = -(lambda_ + mu_);

  // It is very likely that we are able to compute the eigen values and vector from the one of the simple model.
  // For now however, we will use a numerical diagonalization:
  AbstractSubstitutionModel::updateMatrices_();
  // We do not use the one from  AbstractReversibleSubstitutionModel, since we already computed the generator.
}

/******************************************************************************/

double RE08::Pij_t(size_t i, size_t j, double d) const
{
  double f = (lambda_ == 0 && mu_ == 0) ? 1. : lambda_ / (lambda_ + mu_);
  if (i < size_ - 1 && j < size_ - 1)
  {
    return (simpleModel_->Pij_t(i, j, d) - simpleModel_->freq(j)) * exp(-mu_ * d)
           + freq_[j] + (simpleModel_->freq(j) - freq_[j]) * exp(-(lambda_ + mu_) * d);
  }
  else
  {
    if (i == size_ - 1)
    {
      if (j < size_ - 1)
      {
        return freq_[j] * (1. - exp(-(lambda_ + mu_) * d));
      }
      else
      {
        return 1. - f * (1. - exp(-(lambda_ + mu_) * d));
      }
    }
    else
    {
      return freq_[j] * (1. - exp(-(lambda_ + mu_) * d));
    }
  }
}

/******************************************************************************/

double RE08::dPij_dt(size_t i, size_t j, double d) const
{
  double f = (lambda_ == 0 && mu_ == 0) ? 1. : lambda_ / (lambda_ + mu_);
  if (i < size_ - 1 && j < size_ - 1)
  {
    return simpleModel_->dPij_dt(i, j, d) * exp(-mu_ * d)
           - mu_ * (simpleModel_->Pij_t(i, j, d) - simpleModel_->freq(j)) * exp(-mu_ * d)
           - (lambda_ + mu_) * (simpleModel_->freq(j) - freq_[j]) * exp(-(lambda_ + mu_) * d);
  }
  else
  {
    if (i == size_ - 1)
    {
      if (j < size_ - 1)
      {
        return (lambda_ + mu_) * freq_[j] * exp(-(lambda_ + mu_) * d);
      }
      else
      {
        return -f * (lambda_ + mu_) * exp(-(lambda_ + mu_) * d);
      }
    }
    else
    {
      return (lambda_ + mu_) * freq_[j] * exp(-(lambda_ + mu_) * d);
    }
  }
}

/******************************************************************************/

double RE08::d2Pij_dt2(size_t i, size_t j, double d) const
{
  double f = (lambda_ == 0 && mu_ == 0) ? 1. : lambda_ / (lambda_ + mu_);
  if (i < size_ - 1 && j < size_ - 1)
  {
    return simpleModel_->d2Pij_dt2(i, j, d) * exp(-mu_ * d)
           - 2 * mu_ * simpleModel_->dPij_dt(i, j, d) * exp(-mu_ * d)
           + mu_ * mu_ * (simpleModel_->Pij_t(i, j, d) - simpleModel_->freq(j)) * exp(-mu_ * d)
           + (lambda_ + mu_) * (lambda_ + mu_) * (simpleModel_->freq(j) - freq_[j]) * exp(-(lambda_ + mu_) * d);
  }
  else
  {
    if (i == size_ - 1)
    {
      if (j < size_ - 1)
      {
        return -(lambda_ + mu_) * (lambda_ + mu_) * freq_[j] * exp(-(lambda_ + mu_) * d);
      }
      else
      {
        return f * (lambda_ + mu_) * (lambda_ + mu_) * exp(-(lambda_ + mu_) * d);
      }
    }
    else
    {
      return -(lambda_ + mu_) * (lambda_ + mu_) * freq_[j] * exp(-(lambda_ + mu_) * d);
    }
  }
}

/******************************************************************************/

const Matrix<double>& RE08::getPij_t(double d) const
{
  RowMatrix<double> simpleP = simpleModel_->getPij_t(d);
  double f = (lambda_ == 0 && mu_ == 0) ? 1. : lambda_ / (lambda_ + mu_);
  for (size_t i = 0; i < size_ - 1; i++)
  {
    for (size_t j = 0; j < size_ - 1; j++)
    {
      p_(i, j) = (simpleP(i, j) - simpleModel_->freq(j)) * exp(-mu_ * d)
          + freq_[j] + (simpleModel_->freq(j) - freq_[j]) * exp(-(lambda_ + mu_) * d);
    }
  }
  for (size_t j = 0; j < size_ - 1; j++)
  {
    p_(size_ - 1, j) = freq_[j] * (1. - exp(-(lambda_ + mu_) * d));
  }
  p_(size_ - 1, size_ - 1) = 1. - f * (1. - exp(-(lambda_ + mu_) * d));
  for (size_t i = 0; i < size_ - 1; i++)
  {
    p_(i, size_ - 1) = freq_[size_ - 1] * (1. - exp(-(lambda_ + mu_) * d));
  }
  return p_;
}

/******************************************************************************/

const Matrix<double>& RE08::getdPij_dt(double d) const
{
  RowMatrix<double> simpleP = simpleModel_->getPij_t(d);
  RowMatrix<double> simpleDP = simpleModel_->getdPij_dt(d);
  double f = (lambda_ == 0 && mu_ == 0) ? 1. : lambda_ / (lambda_ + mu_);
  for (size_t i = 0; i < size_ - 1; i++)
  {
    for (size_t j = 0; j < size_ - 1; j++)
    {
      p_(i, j) = simpleDP(i, j) * exp(-mu_ * d)
          - mu_ * (simpleP(i, j) - simpleModel_->freq(j)) * exp(-mu_ * d)
          - (lambda_ + mu_) * (simpleModel_->freq(j) - freq_[j]) * exp(-(lambda_ + mu_) * d);
    }
  }
  for (size_t j = 0; j < size_ - 1; j++)
  {
    p_(size_ - 1, j) = (lambda_ + mu_) * freq_[j] * exp(-(lambda_ + mu_) * d);
  }
  p_(size_ - 1, size_ - 1) = -f * (lambda_ + mu_) * exp(-(lambda_ + mu_) * d);
  for (size_t i = 0; i < size_ - 1; i++)
  {
    p_(i, size_ - 1) = (lambda_ + mu_) * freq_[size_ - 1] * exp(-(lambda_ + mu_) * d);
  }
  return p_;
}

/******************************************************************************/

const Matrix<double>& RE08::getd2Pij_dt2(double d) const
{
  RowMatrix<double> simpleP = simpleModel_->getPij_t(d);
  RowMatrix<double> simpleDP = simpleModel_->getdPij_dt(d);
  RowMatrix<double> simpleD2P = simpleModel_->getd2Pij_dt2(d);
  double f = (lambda_ == 0 && mu_ == 0) ? 1. : lambda_ / (lambda_ + mu_);
  for (size_t i = 0; i < size_ - 1; i++)
  {
    for (size_t j = 0; j < size_ - 1; j++)
    {
      p_(i, j) = simpleD2P(i, j) * exp(-mu_ * d)
          - 2 * mu_ * simpleDP(i, j) * exp(-mu_ * d)
          + mu_ * mu_ * (simpleP(i, j) - simpleModel_->freq(j)) * exp(-mu_ * d)
          + (lambda_ + mu_) * (lambda_ + mu_) * (simpleModel_->freq(j) - freq_[j]) * exp(-(lambda_ + mu_) * d);
    }
  }
  for (size_t j = 0; j < size_ - 1; j++)
  {
    p_(size_ - 1, j) = -(lambda_ + mu_) * (lambda_ + mu_) * freq_[j] * exp(-(lambda_ + mu_) * d);
  }
  p_(size_ - 1, size_ - 1) = f * (lambda_ + mu_) * (lambda_ + mu_) * exp(-(lambda_ + mu_) * d);
  for (size_t i = 0; i < size_ - 1; i++)
  {
    p_(i, size_ - 1) = -(lambda_ + mu_) * (lambda_ + mu_) * freq_[size_ - 1] * exp(-(lambda_ + mu_) * d);
  }
  return p_;
}

/******************************************************************************/

double RE08::getInitValue(size_t i, int state) const
{
  if (i >= size_)
    throw IndexOutOfBoundsException("RE08::getInitValue", i, 0, size_ - 1);
  if (state < -1 || !getAlphabet()->isIntInAlphabet(state))
    throw BadIntException(state, "RE08::getInitValue. Character " + alphabet_->intToChar(state) + " is not allowed in model.", alphabet_.get());
  if (i == size_ - 1 && state == -1)
    return 1.;
  vector<int> states = alphabet_->getAlias(state);
  for (size_t j = 0; j < states.size(); j++)
  {
    if ((int)i == states[j])
      return 1.;
  }
  return 0.;
}

/******************************************************************************/

void RE08::setNamespace(const string& prefix)
{
  AbstractSubstitutionModel::setNamespace(prefix);
  // We also need to update the namespace of the nested model:
  simpleModel_->setNamespace(prefix + nestedPrefix_);
}

/******************************************************************************/
