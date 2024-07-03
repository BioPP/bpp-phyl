// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "DSO78.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

#include <map>

using namespace std;

/******************************************************************************/

DSO78::DSO78(std::shared_ptr<const ProteicAlphabet> alpha) :
  AbstractParameterAliasable("DSO78."),
  AbstractReversibleProteinSubstitutionModel(alpha, std::shared_ptr<const StateMapInterface>(new CanonicalStateMap(alpha, false)), "DSO78."),
  freqSet_(nullptr)
{
  #include "__DSO78ExchangeabilityCode"
  #include "__DSO78FrequenciesCode"
  freqSet_.reset(new FixedProteinFrequencySet(alpha, freq_));
  updateMatrices_();
}

DSO78::DSO78(
    std::shared_ptr<const ProteicAlphabet> alpha,
    std::unique_ptr<ProteinFrequencySetInterface> freqSet,
    bool initFreqs) :
  AbstractParameterAliasable("DSO78+F."),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "DSO78+F."),
  freqSet_(std::move(freqSet))
{
  #include "__DSO78ExchangeabilityCode"
  #include "__DSO78FrequenciesCode"
  freqSet_->setNamespace("DSO78+F." + freqSet_->getNamespace());
  if (initFreqs)
    freqSet_->setFrequencies(freq_);
  else
    freq_ = freqSet_->getFrequencies();
  addParameters_(freqSet_->getParameters());
  updateMatrices_();
}

/******************************************************************************/

void DSO78::setFreqFromData(const SequenceDataInterface& data, double pseudoCount)
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
