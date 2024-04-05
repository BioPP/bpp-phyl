// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractBiblioSubstitutionModel.h"

using namespace bpp;
using namespace std;

AbstractBiblioTransitionModel::AbstractBiblioTransitionModel(const std::string& prefix) :
  // AbstractParameterAliasable(prefix),
  mapParNamesFromPmodel_(),
  lParPmodel_()
{}

/******************************************************************************/

AbstractBiblioTransitionModel::AbstractBiblioTransitionModel(const AbstractBiblioTransitionModel& model) :
  // AbstractParameterAliasable(model),
  mapParNamesFromPmodel_(model.mapParNamesFromPmodel_),
  lParPmodel_(model.lParPmodel_)
{}

/******************************************************************************/

AbstractBiblioTransitionModel& AbstractBiblioTransitionModel::operator=(const AbstractBiblioTransitionModel& model)
{
  AbstractParameterAliasable::operator=(model);
  mapParNamesFromPmodel_ = model.mapParNamesFromPmodel_;
  lParPmodel_            = model.lParPmodel_;
  return *this;
}


/******************************************************************************/

std::string AbstractBiblioTransitionModel::getParNameFromPmodel(const std::string& name) const
{
  auto it = mapParNamesFromPmodel_.find(name);
  if (it == mapParNamesFromPmodel_.end())
    throw Exception("AbstractBiblioTransitionModel::getParNameFromPmodel : unknown parameter " + name);
  return it->second;
}

/******************************************************************************/

std::string AbstractBiblioTransitionModel::getPmodelParName(const std::string& name) const
{
  for (auto it : mapParNamesFromPmodel_)
  {
    if (it.second == name)
      return it.first;
  }

  throw Exception("AbstractBiblioTransitionModel::getPmodelParName: unknown parameter name " + name);
}


/******************************************************************************/

void AbstractBiblioTransitionModel::updateMatrices_()
{
  for (size_t i = 0; i < lParPmodel_.size(); ++i)
  {
    if (mapParNamesFromPmodel_.find(lParPmodel_[i].getName()) != mapParNamesFromPmodel_.end())
    {
      lParPmodel_[i].setValue(parameter(getParameterNameWithoutNamespace(mapParNamesFromPmodel_[lParPmodel_[i].getName()])).getValue());
    }
  }

  model_().matchParametersValues(lParPmodel_);
}

/******************************************************************************/

void AbstractBiblioTransitionModel::addRateParameter()
{
  model_().addRateParameter();
  addParameter_(new Parameter(getNamespace() + "rate", model().getRate(), Parameter::R_PLUS_STAR));

  mapParNamesFromPmodel_[getNamespace() + "rate"] = "rate";
  lParPmodel_.reset();
  lParPmodel_.addParameters(model().getParameters());
}

/******************************************************************************/

void AbstractBiblioTransitionModel::setNamespace(const string& name)
{
  AbstractParameterAliasable::setNamespace(name);

  map<string, string> mapParNamesFromPmodel_new;

  for (const auto& it : mapParNamesFromPmodel_)
  {
    mapParNamesFromPmodel_new[name + model().getParameterNameWithoutNamespace(it.first)] = it.second;
  }

  mapParNamesFromPmodel_.clear();
  mapParNamesFromPmodel_ = mapParNamesFromPmodel_new;

  model_().setNamespace(name);

  lParPmodel_.reset();
  lParPmodel_.addParameters(model().getParameters());
}


/******************************************************************************/

void AbstractBiblioTransitionModel::setFreq(std::map<int, double>& frequ)
{
  AbstractTotallyWrappedTransitionModel::setFreq(frequ);

  ParameterList pl;
  for (const auto& it : mapParNamesFromPmodel_)
  {
    pl.addParameter(Parameter(getNamespace() + it.second, model().getParameterValue(model().getParameterNameWithoutNamespace(it.first))));
  }

  matchParametersValues(pl);
}


void AbstractBiblioTransitionModel::setFreqFromData(
    const SequenceDataInterface& data, double pseudoCount)
{
  map<int, double> freqs;
  SequenceContainerTools::getFrequencies(data, freqs, pseudoCount);
  setFreq(freqs);
}
