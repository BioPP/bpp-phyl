//
// File: RE08.cpp
// Created by: Julien Dutheil
// Created on: Mon Dec 29 10:15 2008
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "RE08.h"

using namespace bpp;

#include <cmath>

using namespace std;

/******************************************************************************/

RE08::RE08(ReversibleSubstitutionModel* simpleModel, double lambda, double mu) :
  AbstractParameterAliasable("RE08."),
  AbstractReversibleSubstitutionModel(simpleModel->getAlphabet(), new CanonicalStateMap(simpleModel->getStateMap(), true), "RE08."),
  simpleModel_(simpleModel),
  simpleGenerator_(),
  simpleExchangeabilities_(),
  exp_(), p_(), lambda_(lambda), mu_(mu),
  nestedPrefix_("model_" + simpleModel->getNamespace())
{
  addParameter_(new Parameter("RE08.lambda", lambda, &Parameter::R_PLUS));  
  addParameter_(new Parameter("RE08.mu", mu, &Parameter::R_PLUS));
  simpleModel_->setNamespace("RE08." + nestedPrefix_);
  addParameters_(simpleModel->getParameters());
  //We need to overrired this from the AbstractSubstitutionModel constructor,
  //since the number of states in the model is no longer equal to the size of the alphabet.
  size_ = simpleModel->getNumberOfStates() + 1;
  generator_.resize(size_, size_);
  exchangeability_.resize(size_, size_);
  freq_.resize(size_);
  eigenValues_.resize(size_);
  leftEigenVectors_.resize(size_, size_);
  rightEigenVectors_.resize(size_, size_);
  p_.resize(size_,size_);
  updateMatrices();
}

/******************************************************************************/
  
void RE08::updateMatrices()
{
  double f = (lambda_ == 0 && mu_ == 0) ? 1 : lambda_ / (lambda_ + mu_);
  
  // Frequencies:
  for(size_t i = 0; i < size_ - 1; i++)
    freq_[i] = simpleModel_->freq(i) * f;

  freq_[size_-1] = (1. - f);

  simpleGenerator_ = simpleModel_->getGenerator();
  simpleExchangeabilities_ = simpleModel_->getExchangeabilityMatrix();

  // Generator and exchangeabilities:
  for (size_t i = 0; i < size_ - 1; i++)
  {
    for (size_t j = 0; j < size_ - 1; j++)
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

  //It is very likely that we are able to compute the eigen values and vector from the one of the simple model.
  //For now however, we will use a numerical diagonalization:
  AbstractSubstitutionModel::updateMatrices();
  //We do not use the one from  AbstractReversibleSubstitutionModel, since we already computed the generator.
}
  
/******************************************************************************/

double RE08::Pij_t(size_t i, size_t j, double d) const
{
  double f = (lambda_ == 0 && mu_ == 0) ? 1. : lambda_ / (lambda_ + mu_);
  if(i < size_ - 1 && j < size_ - 1)
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
  if(i < size_ - 1 && j < size_ - 1)
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
        return - f * (lambda_ + mu_) * exp(-(lambda_ + mu_) * d);
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
  if(i < size_ - 1 && j < size_ - 1)
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
        return - (lambda_ + mu_) * (lambda_ + mu_) * freq_[j] * exp(-(lambda_ + mu_) * d);
      }
      else
      {
        return f * (lambda_ + mu_) * (lambda_ + mu_) * exp(-(lambda_ + mu_) * d);
      }
    }
    else
    {  
      return - (lambda_ + mu_) * (lambda_ + mu_) * freq_[j] * exp(-(lambda_ + mu_) * d);
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
  for(size_t j = 0; j < size_ - 1; j++)
  {
    p_(size_ - 1, j) = freq_[j] * (1. - exp(-(lambda_ + mu_) * d));
  }
  p_(size_ - 1, size_ - 1) = 1. - f * (1. - exp(-(lambda_ + mu_) * d));
  for(size_t i = 0; i < size_ - 1; i++)
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
  p_(size_ - 1, size_ - 1) = - f * (lambda_ + mu_) * exp(-(lambda_ + mu_) * d);
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
    p_(size_ - 1, j) = - (lambda_ + mu_) * (lambda_ + mu_) * freq_[j] * exp(-(lambda_ + mu_) * d);
  }
  p_(size_ - 1, size_ - 1) = f * (lambda_ + mu_) * (lambda_ + mu_) * exp(-(lambda_ + mu_) * d);
  for (size_t i = 0; i < size_ - 1; i++)
  {  
    p_(i, size_ - 1) = - (lambda_ + mu_) * (lambda_ + mu_) * freq_[size_ - 1] * exp(-(lambda_ + mu_) * d);
  }
  return p_;
}

/******************************************************************************/

double RE08::getInitValue(size_t i, int state) const throw (IndexOutOfBoundsException, BadIntException)
{
  if (i >= size_) throw IndexOutOfBoundsException("RE08::getInitValue", i, 0, size_ - 1);
  if (state < -1 || !getAlphabet()->isIntInAlphabet(state))
    throw BadIntException(state, "RE08::getInitValue. Character " + getAlphabet()->intToChar(state) + " is not allowed in model.");
  if (i == size_ - 1 && state == -1) return 1.;
  vector<int> states = getAlphabet()->getAlias(state);
  for (size_t j = 0; j < states.size(); j++)
    if ((int)i == states[j]) return 1.;
  return 0.;
}

/******************************************************************************/

void RE08::setNamespace(const string& prefix)
{
  AbstractSubstitutionModel::setNamespace(prefix);
  //We also need to update the namespace of the nested model:
  simpleModel_->setNamespace(prefix + nestedPrefix_);
}

/******************************************************************************/

