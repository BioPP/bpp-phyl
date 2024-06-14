// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "LG08.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

#include <map>

using namespace std;

/******************************************************************************/

LG08::LG08(
    shared_ptr<const ProteicAlphabet> alpha) :
  AbstractParameterAliasable("LG08."),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "LG08."),
  freqSet_(nullptr)
{
  #include "__LG08ExchangeabilityCode"
  #include "__LG08FrequenciesCode"
  freqSet_.reset(new FixedProteinFrequencySet(alpha, freq_));
  updateMatrices_();
}

LG08::LG08(
    shared_ptr<const ProteicAlphabet> alpha,
    unique_ptr<ProteinFrequencySetInterface> freqSet,
    bool initFreqs) :
  AbstractParameterAliasable("LG08+F."),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "LG08+F."),
  freqSet_(std::move(freqSet))
{
  #include "__LG08ExchangeabilityCode"
  #include "__LG08FrequenciesCode"
  freqSet_->setNamespace("LG08+F." + freqSet_->getNamespace());
  if (initFreqs)
    freqSet_->setFrequencies(freq_);
  else
    freq_ = freqSet_->getFrequencies();

  addParameters_(freqSet_->getParameters());

  updateMatrices_();
}

/******************************************************************************/

void LG08::setFreqFromData(const SequenceDataInterface& data, double pseudoCount)
{
  map<int, double> counts;
  SequenceContainerTools::getFrequencies(data, counts, pseudoCount);
  for (auto i : counts)
  {
    freq_[(size_t)i.first] = i.second;
  }

  freqSet_->setFrequencies(freq_);
  // Update parameters and re-compute generator and eigen values:
  matchParametersValues(freqSet_->getParameters());
}

/******************************************************************************/
