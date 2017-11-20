//
// File: AbstractCodonBGCSubstitutionModel.cpp
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

#include "AbstractCodonBGCSubstitutionModel.h"
#include <Bpp/Numeric/NumConstants.h>

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractCodonBGCSubstitutionModel::AbstractCodonBGCSubstitutionModel(
  const GeneticCode* pgencode,
  const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  pgencode_(pgencode),
  B_(0),
  S_(0)
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
  int epsilon = pgencode_->getSourceAlphabet()->getGCinCodon(static_cast<int>(j))
    - pgencode_->getSourceAlphabet()->getGCinCodon(static_cast<int>(i));

  switch (epsilon)
  {
  case 0:
    return (pgencode_->areSynonymous(static_cast<int>(i), static_cast<int>(j))?1.
            :(S_==0?1.:S_/(1-exp(-S_))));
  case 1:
    return (pgencode_->areSynonymous(static_cast<int>(i), static_cast<int>(j))
            ?(B_==0?1:B_/(1-exp(-B_)))
            :(B_+S_==0?1.:(B_+S_)/(1-exp(-(B_+S_)))));
  case -1:
    return (pgencode_->areSynonymous(static_cast<int>(i), static_cast<int>(j))
            ?(B_==0?1:-B_/(1-exp(B_)))
            :(-B_+S_==0?1.:(-B_+S_)/(1-exp(B_-S_))));
  }
  return 0;
}

