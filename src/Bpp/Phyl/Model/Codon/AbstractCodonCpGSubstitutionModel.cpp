//
// File: AbstractCodonCpGSubstitutionModel.cpp
// Created by:  Laurent Gueguen
// Created on: Feb 2009
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)
   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#include "AbstractCodonCpGSubstitutionModel.h"
#include <Bpp/Numeric/NumConstants.h>

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractCodonCpGSubstitutionModel::AbstractCodonCpGSubstitutionModel(
  const CodonAlphabet& alphabet,
  const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  rho_(1),
  stateMap_(new CanonicalStateMap(&alphabet, false))
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

  return ((si%16==7  && (si-sj==2 || sj-si==8)) || ((si-1)/4==6 && (si-sj==8 || (sj-si==32))))?rho_:1;
}

