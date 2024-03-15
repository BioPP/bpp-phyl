// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <string>

#include "AbstractMixedTransitionModel.h"

using namespace bpp;
using namespace std;


AbstractMixedTransitionModel::AbstractMixedTransitionModel(
    shared_ptr<const Alphabet> alpha,
    shared_ptr<const StateMapInterface> stateMap,
    const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractTransitionModel(alpha, stateMap, prefix),
  modelsContainer_(),
  vProbas_(),
  vRates_()
{}

AbstractMixedTransitionModel::AbstractMixedTransitionModel(const AbstractMixedTransitionModel& msm) :
  AbstractParameterAliasable(msm),
  AbstractTransitionModel(msm),
  modelsContainer_(),
  vProbas_(),
  vRates_()
{
  for (size_t i = 0; i < msm.modelsContainer_.size(); ++i)
  {
    modelsContainer_.push_back(std::shared_ptr<TransitionModelInterface>(msm.modelsContainer_[i]->clone()));
    vProbas_.push_back(msm.vProbas_[i]);
    vRates_.push_back(msm.vRates_[i]);
  }
}

AbstractMixedTransitionModel& AbstractMixedTransitionModel::operator=(const AbstractMixedTransitionModel& model)
{
  AbstractTransitionModel::operator=(model);

  // Clear existing containers:
  modelsContainer_.clear();
  vProbas_.clear();
  vRates_.clear();

  for (size_t i = 0; i < model.modelsContainer_.size(); ++i)
  {
    modelsContainer_.push_back(std::shared_ptr<TransitionModelInterface>(model.modelsContainer_[i]->clone()));
    vProbas_.push_back(model.vProbas_[i]);
    vRates_.push_back(model.vRates_[i]);
  }

  return *this;
}

const Matrix<double>& AbstractMixedTransitionModel::getPij_t(double t) const
{
  vector<const Matrix<double>* > vM;
  double sP = 0;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
  {
    vM.push_back(&modelsContainer_[n]->getPij_t(t));
    sP += vProbas_[n];
  }

  for (unsigned int i = 0; i < getNumberOfStates(); i++)
  {
    for (unsigned int j = 0; j < getNumberOfStates(); j++)
    {
      double x = 0;
      for (unsigned int n = 0; n < modelsContainer_.size(); n++)
      {
        x += (*vM[n])(i, j) * vProbas_[n];
      }
      pijt_(i, j) = x / sP;
    }
  }
  return pijt_;
}


const Matrix<double>& AbstractMixedTransitionModel::getdPij_dt(double t) const
{
  vector<const Matrix<double>* > vM;
  double sP = 0;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
  {
    vM.push_back(&modelsContainer_[n]->getdPij_dt(t));
    sP += vProbas_[n];
  }

  for (unsigned int i = 0; i < getNumberOfStates(); i++)
  {
    for (unsigned int j = 0; j < getNumberOfStates(); j++)
    {
      double x = 0;
      for (unsigned int n = 0; n < modelsContainer_.size(); n++)
      {
        x += (*vM[n])(i, j) * vProbas_[n];
      }
      dpijt_(i, j) = x / sP;
    }
  }
  return dpijt_;
}


const Matrix<double>& AbstractMixedTransitionModel::getd2Pij_dt2(double t) const
{
  vector<const Matrix<double>* > vM;
  double sP = 0;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
  {
    vM.push_back(&modelsContainer_[n]->getd2Pij_dt2(t));
    sP += vProbas_[n];
  }

  for (unsigned int i = 0; i < getNumberOfStates(); i++)
  {
    for (unsigned int j = 0; j < getNumberOfStates(); j++)
    {
      double x = 0;
      for (unsigned int n = 0; n < modelsContainer_.size(); n++)
      {
        x += (*vM[n])(i, j) * vProbas_[n];
      }
      d2pijt_(i, j) = x / sP;
    }
  }
  return d2pijt_;
}


void AbstractMixedTransitionModel::setRate(double rate)
{
  AbstractTransitionModel::setRate(rate);

  double sum = 0;
  double sP = 0;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
  {
    sum += vRates_[n] * vProbas_[n];
    sP += vProbas_[n];
  }
  sum /= sP;

  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
  {
    vRates_[n] *= rate_ / sum;
    modelsContainer_[n]->setRate(vRates_[n]);
  }
}

void AbstractMixedTransitionModel::setVRates(const Vdouble& vd)
{
  if (vd.size() != modelsContainer_.size())
    throw Exception("AbstractMixedTransitionModel::setVRates  bad size of Vdouble argument.");

  for (unsigned int i = 0; i < vd.size(); i++)
  {
    vRates_[i] = vd[i];
  }

  normalizeVRates();
}

void AbstractMixedTransitionModel::normalizeVRates()
{
  double sum = 0;
  double sP = 0;
  for (unsigned int i = 0; i < vRates_.size(); i++)
  {
    sum += vRates_[i] * vProbas_[i];
    sP += vProbas_[i];
  }
  sum /= sP;

  for (unsigned int i = 0; i < vRates_.size(); i++)
  {
    vRates_[i] *= rate_ / sum;
    modelsContainer_[i]->setRate(vRates_[i]);
  }
}
