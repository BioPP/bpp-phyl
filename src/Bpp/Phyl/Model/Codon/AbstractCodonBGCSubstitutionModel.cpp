// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/NumConstants.h>

#include "AbstractCodonBGCSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractCodonBGCSubstitutionModel::AbstractCodonBGCSubstitutionModel(
  shared_ptr<const GeneticCode> pgencode,
  const string& prefix) :
  AbstractParameterAliasable(prefix),
  pgencode_(pgencode),
  B_(0),
  S_(0),
  stateMap_(new CanonicalStateMap(pgencode->getSourceAlphabet(), false))
{
  addParameter_(new Parameter(prefix + "S", 0));
  addParameter_(new Parameter(prefix + "B", 0));
}

void AbstractCodonBGCSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  B_ = getParameterValue("B");
  S_ = getParameterValue("S");
}

double AbstractCodonBGCSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  int si(stateMap_->getAlphabetStateAsInt(i)), sj(stateMap_->getAlphabetStateAsInt(j));

  int epsilon = pgencode_->codonAlphabet().getGCinCodon(sj)
                - pgencode_->codonAlphabet().getGCinCodon(si);

  switch (epsilon)
  {
  case 0:
    return pgencode_->areSynonymous(si, sj) ? 1.
            : (S_ == 0 ? 1. : S_ / (1 - exp(-S_)));
  case 1:
    return pgencode_->areSynonymous(si, sj)
            ? (B_ == 0 ? 1 : B_ / (1 - exp(-B_)))
            : (B_ + S_ == 0 ? 1. : (B_ + S_) / (1 - exp(-(B_ + S_))));
  case -1:
    return pgencode_->areSynonymous(si, sj)
            ? (B_ == 0 ? 1 : -B_ / (1 - exp(B_)))
            : (-B_ + S_ == 0 ? 1. : (-B_ + S_) / (1 - exp(B_ - S_)));
  }
  return 0;
}
