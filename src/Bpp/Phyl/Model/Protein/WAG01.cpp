// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "WAG01.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

WAG01::WAG01(shared_ptr<const ProteicAlphabet> alpha) :
  AbstractParameterAliasable("WAG01."),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "WAG01."),
  freqSet_(nullptr)
{
  #include "__WAG01ExchangeabilityCode"
  #include "__WAG01FrequenciesCode"
  freqSet_.reset(new FixedProteinFrequencySet(alpha, freq_));
  updateMatrices_();
}

WAG01::WAG01(
    shared_ptr<const ProteicAlphabet> alpha,
    unique_ptr<ProteinFrequencySetInterface> freqSet,
    bool initFreqs) :
  AbstractParameterAliasable("WAG01+F."),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "WAG01+F."),
  freqSet_(std::move(freqSet))
{
  #include "__WAG01ExchangeabilityCode"
  #include "__WAG01FrequenciesCode"
  freqSet_->setNamespace("WAG01+F." + freqSet_->getNamespace());
  if (initFreqs)
    freqSet_->setFrequencies(freq_);
  else
    freq_ = freqSet_->getFrequencies();
  addParameters_(freqSet_->getParameters());
  updateMatrices_();
}

/******************************************************************************/

void WAG01::setFreqFromData(const SequenceDataInterface& data, double pseudoCount)
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
