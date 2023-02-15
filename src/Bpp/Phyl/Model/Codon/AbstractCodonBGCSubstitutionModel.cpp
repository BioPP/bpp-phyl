//
// File: AbstractCodonBGCSubstitutionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: 2009-02-08 00:00:00
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
