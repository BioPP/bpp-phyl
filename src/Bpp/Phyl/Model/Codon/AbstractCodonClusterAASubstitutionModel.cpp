// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/NumConstants.h>

#include "AbstractCodonClusterAASubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractCodonClusterAASubstitutionModel::AbstractCodonClusterAASubstitutionModel(
    std::shared_ptr<const GeneticCode> pgencode,
    const string& prefix,
    const vector<uint>& assign) :
  AbstractParameterAliasable(prefix),
  pgencode_(pgencode),
  omegaR_(1),
  omegaC_(1),
  assign_(assign),
  stateMap_(new CanonicalStateMap(pgencode->getSourceAlphabet(), false))
{
  if (assign_.size() != 20)
    throw BadSizeException("AbstractCodonClusterAASubstitutionModel::AbstractCodonClusterAASubstitutionModel: assign_", assign_.size(), 20);

  addParameter_(new Parameter(prefix + "omegaR", 1, std::make_shared<IntervalConstraint>(NumConstants::SMALL(), 999, true, true)));
  addParameter_(new Parameter(prefix + "omegaC", 1, std::make_shared<IntervalConstraint>(NumConstants::SMALL(), 999, true, true)));
}

void AbstractCodonClusterAASubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  omegaR_ = getParameterValue("omegaR");
  omegaC_ = getParameterValue("omegaC");
}

double AbstractCodonClusterAASubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  int si(stateMap_->getAlphabetStateAsInt(i)), sj(stateMap_->getAlphabetStateAsInt(j));

  return pgencode_->areSynonymous(si, sj) ? 1 :
         (assign_[static_cast<size_t>(pgencode_->translate(si))] == assign_[static_cast<size_t>(pgencode_->translate(sj))] ? omegaC_ :
         omegaR_);
}
