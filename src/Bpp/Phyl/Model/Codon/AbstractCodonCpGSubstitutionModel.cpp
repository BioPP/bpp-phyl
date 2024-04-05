// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/NumConstants.h>

#include "AbstractCodonCpGSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractCodonCpGSubstitutionModel::AbstractCodonCpGSubstitutionModel(
    std::shared_ptr<const CodonAlphabet> alphabet,
    const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  rho_(1),
  stateMap_(new CanonicalStateMap(alphabet, false))
{
  addParameter_(new Parameter(prefix + "rho", 1, std::make_shared<IntervalConstraint>(NumConstants::SMALL(), 999, true, true)));
}

void AbstractCodonCpGSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  rho_ = getParameterValue("rho");
}

double AbstractCodonCpGSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  int si(stateMap_->getAlphabetStateAsInt(i)), sj(stateMap_->getAlphabetStateAsInt(j));

  return ((si % 16 == 7  && (si - sj == 2 || sj - si == 8)) || ((si - 1) / 4 == 6 && (si - sj == 8 || (sj - si == 32)))) ? rho_ : 1;
}
