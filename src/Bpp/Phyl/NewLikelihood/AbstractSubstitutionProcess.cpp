//
// File: AbstractSubstitutionProcess.cpp
// Created by: Julien Dutheil
// Created on: Tue Marc 22 21:17 2013
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

#include "AbstractSubstitutionProcess.h"

using namespace bpp;
using namespace std;

AbstractSubstitutionProcess::AbstractSubstitutionProcess(ParametrizableTree* tree, size_t nbClasses, const string& prefix) :
  AbstractParameterAliasable(prefix),
  pTree_(tree),
  nbClasses_(nbClasses),
  probabilities_(),
  probabilitiesD1_(),
  probabilitiesD2_(),
  computeProbabilities_(),
  computeProbabilitiesD1_(),
  computeProbabilitiesD2_()
{
  if (!tree)
    throw Exception("AbstractSubstitutionProcess. A tree instance must be provided.");
    
  // Set array sizes:
  size_t n = tree->getNumberOfBranches() * nbClasses;
  probabilities_.resize(n);
  probabilitiesD1_.resize(n);
  probabilitiesD2_.resize(n);
  computeProbabilities_.resize(n);
  computeProbabilitiesD1_.resize(n);
  computeProbabilitiesD2_.resize(n);
  for (size_t i = 0; i < n; ++i) {
    computeProbabilities_[i] = true;
    computeProbabilitiesD1_[i] = true;
    computeProbabilitiesD2_[i] = true;
  }
}
 
AbstractSubstitutionProcess::AbstractSubstitutionProcess(const AbstractSubstitutionProcess& asp) :
  AbstractParameterAliasable(asp),
  pTree_(asp.pTree_->clone()),
  nbClasses_(asp.nbClasses_),
  probabilities_(asp.probabilities_),
  probabilitiesD1_(asp.probabilitiesD1_),
  probabilitiesD2_(asp.probabilitiesD2_),
  computeProbabilities_(asp.computeProbabilities_),
  computeProbabilitiesD1_(asp.computeProbabilitiesD1_),
  computeProbabilitiesD2_(asp.computeProbabilitiesD2_)
{}

AbstractSubstitutionProcess& AbstractSubstitutionProcess::operator=(const AbstractSubstitutionProcess& asp)
{
  AbstractParameterAliasable::operator=(*this);
  
  pTree_.reset(pTree_->clone());
  nbClasses_ = asp.nbClasses_;
  probabilities_ = asp.probabilities_;
  probabilitiesD1_ = asp.probabilitiesD1_;
  probabilitiesD2_ = asp.probabilitiesD2_;
  computeProbabilities_ = asp.computeProbabilities_;
  computeProbabilitiesD1_ = asp.computeProbabilitiesD1_;
  computeProbabilitiesD2_ = asp.computeProbabilitiesD2_;
  return *this;
}

size_t AbstractSubstitutionProcess::getModelIndex_(int nodeId, size_t modelClass) const throw (NodeNotFoundException, IndexOutOfBoundsException)
{
  size_t i = pTree_->getNodeIndex(nodeId);
  if (modelClass >= nbClasses_)
    throw IndexOutOfBoundsException("AbstractSubstitutionProcess::getModelIndex_().", modelClass, 0, nbClasses_);
  return i * nbClasses_ + modelClass;
}

void AbstractSubstitutionProcess::fireParameterChanged(const ParameterList& pl)
{
  AbstractParameterAliasable::fireParameterChanged(pl);

  pTree_->matchParametersValues(pl);
  
  for (size_t i = 0; i < computeProbabilities_.size(); ++i) {
    if (modelChangesWithParameter_(i, pl)) {
      computeProbabilities_[i] = true;
      computeProbabilitiesD1_[i] = true;
      computeProbabilitiesD2_[i] = true;
    }
  }
}

