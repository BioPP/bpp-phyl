// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "D1WalkSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

#include <cmath>
#include <map>

using namespace std;

/******************************************************************************/

D1WalkSubstitutionModel::D1WalkSubstitutionModel(std::shared_ptr<const IntegerAlphabet> alpha, unsigned short method) :
  AbstractParameterAliasable("D1Walk."),
  AbstractReversibleSubstitutionModel(alpha, std::shared_ptr<StateMapInterface>(new CanonicalStateMap(alpha, false)), "D1Walk."),
  freqSet_(0)
{
  freqSet_ = std::make_shared<FullFrequencySet>(getStateMap(), true, method);
  freqSet_->setNamespace("D1Walk.");
  computeFrequencies(false);

  addParameters_(freqSet_->getParameters());
  // Exchangeability Matrix:
  for (unsigned int i = 0; i < size_-1; i++)
    exchangeability_(i, i+1) = 1;

  for (unsigned int i = 1; i < size_; i++)
    exchangeability_(i, i-1) = 1;
  
  updateMatrices_();
}

/******************************************************************************/

void D1WalkSubstitutionModel::updateMatrices_()
{
  // Frequencies:
  freq_ = freqSet_->getFrequencies();

  AbstractReversibleSubstitutionModel::updateMatrices_();
}


/******************************************************************************/

void D1WalkSubstitutionModel::setFreq(std::map<int, double>& freqs)
{
  for (auto i : freqs)
  {
    freq_[(size_t)i.first] = i.second;
  }

  freqSet_->setFrequencies(freq_);
  // Update parameters and re-compute generator and eigen values:
  matchParametersValues(freqSet_->getParameters());
}

/******************************************************************************/
