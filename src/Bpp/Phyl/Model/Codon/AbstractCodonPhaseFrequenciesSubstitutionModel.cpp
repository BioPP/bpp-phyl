// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../FrequencySet/NucleotideFrequencySet.h"
#include "AbstractCodonPhaseFrequenciesSubstitutionModel.h"

#include <memory>

using namespace bpp;
using namespace std;

/******************************************************************************/

AbstractCodonPhaseFrequenciesSubstitutionModel::AbstractCodonPhaseFrequenciesSubstitutionModel(
    unique_ptr<CodonFrequencySetInterface> pfreq,
    const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  posFreqSet_(),
  freqName_("")
{
  if (dynamic_cast<CodonFromUniqueFrequencySet*>(pfreq.get())
      || dynamic_cast<CodonFromIndependentFrequencySet*>(pfreq.get()))
    posFreqSet_.reset(pfreq->clone());
  else
  {
    vector<unique_ptr<FrequencySetInterface>> vFS;
    if (dynamic_cast<FixedCodonFrequencySet*>(pfreq.get()))
    {
      for (unsigned int i = 0; i < 3; ++i)
      {
        vFS.push_back(make_unique<FixedNucleotideFrequencySet>(pfreq->getCodonAlphabet()->getNucleicAlphabet()));
      }
    }
    else
    {
      for (unsigned int i = 0; i < 3; ++i)
      {
        vFS.push_back(make_unique<FullNucleotideFrequencySet>(pfreq->getCodonAlphabet()->getNucleicAlphabet()));
      }
    }
    posFreqSet_.reset(new CodonFromIndependentFrequencySet(
                        pfreq->getGeneticCode(),
                        vFS, ""));

    posFreqSet_->setFrequencies(pfreq->getFrequencies());
  }

  freqName_ = pfreq->getNamespace();
  posFreqSet_->setNamespace(prefix + pfreq->getNamespace());
  addParameters_(posFreqSet_->getParameters());
  fireParameterChanged(posFreqSet_->getParameters());
}

/******************************************************************************/

AbstractCodonPhaseFrequenciesSubstitutionModel::~AbstractCodonPhaseFrequenciesSubstitutionModel() {}

/******************************************************************************/

void AbstractCodonPhaseFrequenciesSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  posFreqSet_->matchParametersValues(parameters);
}

/******************************************************************************/

void AbstractCodonPhaseFrequenciesSubstitutionModel::setFreq(map<int, double>& frequencies)
{
  posFreqSet_->setFrequenciesFromAlphabetStatesFrequencies(frequencies);
  matchParametersValues(posFreqSet_->getParameters());
}

/******************************************************************************/

double AbstractCodonPhaseFrequenciesSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  size_t i2(i), j2(j);

  double x = 1.;
  for (size_t k = 0; k < 3; k++)
  {
    if ((i2 % 4) != (j2 % 4))
      x *= dynamic_cast<WordFrequencySetInterface&>(*posFreqSet_).frequencySetForLetter(2 - k).getFrequencies()[j2 % 4];
    i2 /= 4;
    j2 /= 4;
  }
  return x;
}

/******************************************************************************/

