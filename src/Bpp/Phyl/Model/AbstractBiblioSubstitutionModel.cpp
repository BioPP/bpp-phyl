//
// File: AbstractBiblioSubstitutionModel.cpp
// Created by: Laurent Guéguen
// Created on: lundi 11 juillet 2011, à 21h 12
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

#include "AbstractBiblioSubstitutionModel.h"

using namespace bpp;
using namespace std;

AbstractBiblioTransitionModel::AbstractBiblioTransitionModel(const std::string& prefix) : AbstractParameterAliasable(prefix),
  mapParNamesFromPmodel_(),
  lParPmodel_()
{}

/******************************************************************************/

AbstractBiblioTransitionModel::AbstractBiblioTransitionModel(const AbstractBiblioTransitionModel& model) :
  AbstractParameterAliasable(model),
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
  if (it==mapParNamesFromPmodel_.end())
    throw Exception("AbstractBiblioTransitionModel::getParNameFromPmodel : unknown parameter " + name);
  return it->second;
}

/******************************************************************************/

std::string AbstractBiblioTransitionModel::getPmodelParName(const std::string& name) const
{
  for (auto it : mapParNamesFromPmodel_)
    if (it.second == name)
      return it.first;

  throw Exception("AbstractBiblioTransitionModel::getParNameFromPmodel: unknown parameter name " + name);
}


/******************************************************************************/

void AbstractBiblioTransitionModel::updateMatrices()
{
  for (size_t i = 0; i < lParPmodel_.size(); i++)
  {
    if (mapParNamesFromPmodel_.find(lParPmodel_[i].getName()) != mapParNamesFromPmodel_.end()) {
      lParPmodel_[i].setValue(getParameter(getParameterNameWithoutNamespace(mapParNamesFromPmodel_[lParPmodel_[i].getName()])).getValue());
    }
  }
  
  getModel().matchParametersValues(lParPmodel_);
}

/******************************************************************************/

void AbstractBiblioTransitionModel::addRateParameter()
{
  getModel().addRateParameter();
  addParameter_(new Parameter(getNamespace() + "rate", getModel().getRate(), Parameter::R_PLUS_STAR));
  
  mapParNamesFromPmodel_[getNamespace() + "rate"] = "rate";
  lParPmodel_.reset();
  lParPmodel_.addParameters(getModel().getParameters());
}

/******************************************************************************/

void AbstractBiblioTransitionModel::setNamespace(const std::string& name)
{
  AbstractParameterAliasable::setNamespace(name);

  std::map<std::string, std::string> mapParNamesFromPmodel_new;

  for (const auto& it : mapParNamesFromPmodel_)
    mapParNamesFromPmodel_new[name+getModel().getParameterNameWithoutNamespace(it.first)]=it.second;
  
  mapParNamesFromPmodel_.clear();
  mapParNamesFromPmodel_=mapParNamesFromPmodel_new;

  getModel().setNamespace(name);

  lParPmodel_.reset();
  lParPmodel_.addParameters(getModel().getParameters());
}


/******************************************************************************/

void AbstractBiblioTransitionModel::setFreq(std::map<int, double>& frequ)
{
  AbstractTotallyWrappedTransitionModel::setFreq(frequ);

  ParameterList pl;
  for (const auto& it : mapParNamesFromPmodel_)
  {
    pl.addParameter(Parameter(getNamespace() + it.second, getModel().getParameterValue(getModel().getParameterNameWithoutNamespace(it.first))));
  }

  matchParametersValues(pl);
}


void AbstractBiblioTransitionModel::setFreqFromData(const SequencedValuesContainer& data, double pseudoCount)
{
  map<int, double> freqs;
  SequenceContainerTools::getFrequencies(data, freqs, pseudoCount);
  setFreq(freqs);
}

