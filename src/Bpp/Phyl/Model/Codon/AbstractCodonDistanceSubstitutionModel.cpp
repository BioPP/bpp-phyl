// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/NumConstants.h>

#include "AbstractCodonDistanceSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

AbstractCodonDistanceSubstitutionModel::AbstractCodonDistanceSubstitutionModel(
    std::shared_ptr<const AlphabetIndex2> pdist,
    std::shared_ptr<const GeneticCode> pgencode,
    const string& prefix,
    bool paramSynRate) :
  AbstractParameterAliasable(prefix),
  pdistance_(pdist),
  pgencode_(pgencode),
  alpha_(10000),
  beta_(1),
  gamma_(1),
  stateMap_(new CanonicalStateMap(pgencode->getSourceAlphabet(), false))
{
  if (pdistance_)
    addParameter_(new Parameter(prefix + "alpha", 10000, Parameter::R_PLUS_STAR));

  if (paramSynRate)
    addParameter_(new Parameter(prefix + "gamma", 1, std::make_shared<IntervalConstraint>(NumConstants::SMALL(), 999, true, true)));

  addParameter_(new Parameter(prefix + "beta", 1, Parameter::R_PLUS_STAR));
}

void AbstractCodonDistanceSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  if (pdistance_)
    alpha_ = getParameterValue("alpha");

  if (hasParameter("gamma"))
    gamma_ = getParameterValue("gamma");
  beta_ = getParameterValue("beta");
}

double AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  int si(stateMap_->getAlphabetStateAsInt(i)), sj(stateMap_->getAlphabetStateAsInt(j));

  return pgencode_->areSynonymous(si, sj) ? gamma_ :
         beta_ * (pdistance_ ? exp(-pdistance_->getIndex(
        pgencode_->translate(si),
        pgencode_->translate(sj)) / alpha_) : 1);
}
