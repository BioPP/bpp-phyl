// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractCodonFrequenciesSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

AbstractCodonFrequenciesSubstitutionModel::AbstractCodonFrequenciesSubstitutionModel(
    unique_ptr<CodonFrequencySetInterface> pfreq,
    const string& prefix) :
  AbstractParameterAliasable(prefix),
  pfreqset_(std::move(pfreq)),
  freqName_("")
{
  freqName_ = pfreqset_->getNamespace();
  pfreqset_->setNamespace(prefix + freqName_);
  addParameters_(pfreqset_->getParameters());
}

/******************************************************************************/

AbstractCodonFrequenciesSubstitutionModel::~AbstractCodonFrequenciesSubstitutionModel() {}

/******************************************************************************/

void AbstractCodonFrequenciesSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  pfreqset_->matchParametersValues(parameters);
}

/******************************************************************************/

void AbstractCodonFrequenciesSubstitutionModel::setFreq(map<int, double>& frequencies)
{
  pfreqset_->setFrequenciesFromAlphabetStatesFrequencies(frequencies);
  matchParametersValues(pfreqset_->getParameters());
}

/******************************************************************************/

double AbstractCodonFrequenciesSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  return pfreqset_->getFrequencies()[j];
}

/******************************************************************************/

