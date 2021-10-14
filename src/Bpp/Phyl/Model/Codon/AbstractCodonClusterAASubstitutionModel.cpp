//
// File: AbstractCodonClusterAASubstitutionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: samedi 25 aoÃÂ»t 2018, ÃÂ  09h 06
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#include <Bpp/Numeric/NumConstants.h>

#include "AbstractCodonClusterAASubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractCodonClusterAASubstitutionModel::AbstractCodonClusterAASubstitutionModel(
  const GeneticCode* pgencode,
  const std::string& prefix,
  const std::vector<uint>& assign) :
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
