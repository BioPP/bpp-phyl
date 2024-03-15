// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Text/TextTools.h>
#include <string>

#include "MixtureOfTransitionModels.h"

using namespace bpp;
using namespace std;

MixtureOfTransitionModels::MixtureOfTransitionModels(
  shared_ptr<const Alphabet> alpha,
  vector<unique_ptr<TransitionModelInterface>>& vpModel) :
  AbstractParameterAliasable("Mixture."),
  AbstractTransitionModel(alpha, vpModel.size() ? vpModel[0]->getStateMap() : 0, "Mixture."),
  AbstractMixedTransitionModel(alpha, vpModel.size() ? vpModel[0]->getStateMap() : 0, "Mixture.")
{
  size_t i, nbmod = vpModel.size();
  if (nbmod == 0)
    throw Exception("MixtureOfTransitionModels::MixtureOfTransitionModels : empty vector of models.");
  
  for (i = 0; i < nbmod; ++i)
  {
    if (!vpModel[i])
      throw Exception("Empty model number " + TextTools::toString(i) + " in MixtureOfTransitionModels constructor");
    for (size_t j = i + 1; j < nbmod; ++j)
    {
      if (vpModel[i] == vpModel[j])
        throw Exception("Same model at positions " + TextTools::toString(i) + " and " +
                        TextTools::toString(j) + " in MixtureOfTransitionModels constructor");
    }
  }

  // Initialization of modelsContainer_.

  for (i = 0; i < nbmod; i++)
  {
    modelsContainer_.push_back(std::move(vpModel[i]));
    vProbas_.push_back(1.0 / static_cast<double>(nbmod));
    vRates_.push_back(1.0);
  }

  // Initialization of parameters_.

  // relative rates and probas
  for (i = 0; i < nbmod - 1; i++)
  {
    addParameter_(new Parameter("Mixture.relproba" + TextTools::toString(i + 1), 1.0 / static_cast<double>(nbmod - i), Parameter::PROP_CONSTRAINT_EX));
    addParameter_(new Parameter("Mixture.relrate" + TextTools::toString(i + 1), 1.0 / static_cast<double>(nbmod - i), Parameter::PROP_CONSTRAINT_EX));
  }

  // models parameters

  for (i = 0; i < nbmod; i++)
  {
    modelsContainer_[i]->setNamespace("Mixture." + TextTools::toString(i + 1) + "_" + modelsContainer_[i]->getNamespace());
    addParameters_(modelsContainer_[i]->getParameters());
  }

  updateMatrices_();
}

MixtureOfTransitionModels::MixtureOfTransitionModels(
  shared_ptr<const Alphabet> alpha,
  vector<unique_ptr<TransitionModelInterface>>& vpModel,
  Vdouble& vproba,
  Vdouble& vrate) :
  AbstractParameterAliasable("Mixture."),
  AbstractTransitionModel(alpha, vpModel.size() ? vpModel[0]->getStateMap() : 0, "Mixture."),
  AbstractMixedTransitionModel(alpha, vpModel.size() ? vpModel[0]->getStateMap() : 0, "Mixture.")
{
  size_t i, nbmod = vpModel.size();

  for (i = 0; i < nbmod; ++i)
  {
    if (!vpModel[i])
      throw Exception("Empty model number " + TextTools::toString(i) + " in MixtureOfTransitionModels constructor");
    for (size_t j = i + 1; j < nbmod; ++j)
    {
      if (vpModel[i] == vpModel[j])
        throw Exception("Same model at positions " + TextTools::toString(i) + " and " +
                        TextTools::toString(j) + " in MixtureOfTransitionModels constructor");
    }
  }

  double x = 0;
  double sumrates = 0;

  for (i = 0; i < nbmod; i++)
  {
    if (vrate[i] <= 0)
      throw Exception("Non positive rate: " + TextTools::toString(vrate[i]) + " in MixtureOfTransitionModels constructor.");
    if (vproba[i] <= 0)
      throw Exception("Non positive probability: " + TextTools::toString(vproba[i]) + " in MixtureOfTransitionModels constructor.");
    x += vproba[i];
    sumrates += vrate[i];
  }

  if (fabs(1. - x) > NumConstants::SMALL())
    throw Exception("Probabilities must equal 1 (sum = " + TextTools::toString(x) + ").");


  // Initialization of modelsContainer_.

  for (i = 0; i < nbmod; i++)
  {
    modelsContainer_.push_back(std::move(vpModel[i]));
  }

  // rates & probas

  vProbas_.resize(nbmod);
  vRates_.resize(nbmod);

  // Initialization of parameters_.

  // relative rates and probas
  x = 0;
  double y = 0;

  for (i = 0; i < nbmod - 1; ++i)
  {
    addParameter_(new Parameter("Mixture.relproba" + TextTools::toString(i + 1), vproba[i] / (1 - x), Parameter::PROP_CONSTRAINT_EX));
    x += vproba[i];
    addParameter_(new Parameter("Mixture.relrate" + TextTools::toString(i + 1), vrate[i] / (sumrates - y), Parameter::PROP_CONSTRAINT_EX));
    y += vrate[i];
  }

  // models parameters

  for (i = 0; i < nbmod; ++i)
  {
    modelsContainer_[i]->setNamespace("Mixture." + TextTools::toString(i + 1) + "_" + modelsContainer_[i]->getNamespace());
    addParameters_(modelsContainer_[i]->getParameters());
  }

  updateMatrices_();
}

MixtureOfTransitionModels::MixtureOfTransitionModels(const MixtureOfTransitionModels& msm) :
  AbstractParameterAliasable(msm),
  AbstractTransitionModel(msm),
  AbstractMixedTransitionModel(msm)
{}

MixtureOfTransitionModels& MixtureOfTransitionModels::operator=(const MixtureOfTransitionModels& msm)
{
  AbstractMixedTransitionModel::operator=(msm);

  return *this;
}


MixtureOfTransitionModels::~MixtureOfTransitionModels()
{}

const TransitionModelInterface& MixtureOfTransitionModels::model(const string& name) const
{
  size_t nbmod = getNumberOfModels();

  for (size_t i = 0; i < nbmod; ++i)
  {
    if (nModel(i).getName() == name)
      return nModel(i);
  }

  throw NullPointerException("MixtureOfTransitionModels::model. No model with name '" + name + "'.");
}

void MixtureOfTransitionModels::updateMatrices_()
{
  size_t i, j, nbmod = modelsContainer_.size();

  double x, y;
  x = 1.0;

  for (i = 0; i < nbmod - 1; i++)
  {
    y = getParameterValue("relproba" + TextTools::toString(i + 1));
    vProbas_[i] = x * y;
    x *= 1 - y;
  }
  vProbas_[nbmod - 1] = x;

  x = 1.0;
  double s = 0;
  for (i = 0; i < nbmod - 1; i++)
  {
    y = getParameterValue("relrate" + TextTools::toString(i + 1));
    vRates_[i] = x * y;
    x *= 1 - y;
    s += vProbas_[i] * vRates_[i];
  }

  vRates_[nbmod - 1] = x;
  s += vProbas_[nbmod - 1] * vRates_[nbmod - 1];

  for (i = 0; i < nbmod; i++)
  {
    vRates_[i] /= s;
  }

  // models

  for (i = 0; i < nbmod; i++)
  {
    modelsContainer_[i]->setRate(rate_ * vRates_[i]);
    modelsContainer_[i]->matchParametersValues(getParameters());
  }

  // / freq_

  for (i = 0; i < getNumberOfStates(); i++)
  {
    freq_[i] = 0;
    for (j = 0; j < modelsContainer_.size(); j++)
    {
      freq_[i] += vProbas_[j] * modelsContainer_[j]->freq(i);
    }
  }
}


void MixtureOfTransitionModels::setFreq(std::map<int, double>& m)
{
  ParameterList pl;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
  {
    modelsContainer_[n]->setFreq(m);
    pl.addParameters(modelsContainer_[n]->getParameters());
  }
  matchParametersValues(pl);
}

void MixtureOfTransitionModels::setVRates(const Vdouble& vd)
{
  AbstractMixedTransitionModel::setVRates(vd);

  size_t i, nbmod = modelsContainer_.size();
  double sP = 0;
  for (const auto& rate:vRates_)
  {
    sP += rate;
  }

  for (i = 0; i < nbmod - 1; i++)
  {
    setParameterValue("relrate" + TextTools::toString(i + 1), vRates_[i] / sP);
    sP -= vRates_[i];
  }
}

Vuint MixtureOfTransitionModels::getSubmodelNumbers(const string& desc) const
{
  size_t i;
  for (i = 0; i < getNumberOfModels(); i++)
  {
    if (getNModel(i)->getName() == desc)
      break;
  }
  if (i == getNumberOfModels())
    throw Exception("MixtureOfTransitionModels::getSubmodelNumbers model description do not match " + desc);

  Vuint submodnb(1, uint(i));

  return submodnb;
}
