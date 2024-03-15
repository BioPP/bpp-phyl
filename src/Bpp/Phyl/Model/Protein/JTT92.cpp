// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "JTT92.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

#include <map>

using namespace std;

/******************************************************************************/

JTT92::JTT92(shared_ptr<const ProteicAlphabet> alpha) :
  AbstractParameterAliasable("JTT92."),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "JTT92."),
  freqSet_(nullptr)
{
  #include "__JTT92ExchangeabilityCode"
  #include "__JTT92FrequenciesCode"
  freqSet_.reset(new FixedProteinFrequencySet(alpha, freq_));
  updateMatrices_();
}

JTT92::JTT92(
    shared_ptr<const ProteicAlphabet> alpha,
    unique_ptr<ProteinFrequencySetInterface> freqSet,
    bool initFreqs) :
  AbstractParameterAliasable("JTT92+F."),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "JTT92+F."),
  freqSet_(std::move(freqSet))
{
  #include "__JTT92ExchangeabilityCode"
  #include "__JTT92FrequenciesCode"
  freqSet_->setNamespace("JTT92+F." + freqSet_->getNamespace());
  if (initFreqs)
    freqSet_->setFrequencies(freq_);
  else
    freq_ = freqSet_->getFrequencies();
  addParameters_(freqSet_->getParameters());
  updateMatrices_();
}

/******************************************************************************/

void JTT92::setFreqFromData(const SequenceDataInterface& data, double pseudoCount)
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
