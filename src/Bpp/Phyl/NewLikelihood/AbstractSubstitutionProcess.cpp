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

AbstractSubstitutionProcess::AbstractSubstitutionProcess(ParametrizableTree* tree, size_t nbClasses) :
  pTree_(tree),
  nodeIndex_(),
  nbClasses_(nbClasses),
  probabilities_(),
  probabilitiesD1_(),
  probabilitiesD2_(),
  computeProbability_(),
  computeProbabilityD1_(),
  computeProbabilityD2_()
{
  if (!tree)
    throw Exception("AbstractSubstitutionProcess. A tree instance must be provided.");
    
  //Build node index (NB: if we allow to change the tree, this will have to be recomputed):
  std::vector<int> ids = pTree_->getBranchesId();
  for (size_t i = 0; i < ids.size(); ++i) {
    nodeIndex_[ids[i]] = i;
  }

  // Set array sizes:
  size_t n = tree->getNumberOfBranches() * nbClasses;
  probabilities_.resize(n);
  probabilitiesD1_.resize(n);
  probabilitiesD2_.resize(n);
  computeProbability_.resize(n);
  computeProbabilityD1_.resize(n);
  computeProbabilityD2_.resize(n);
  for (size_t i = 0; i < n; ++i) {
    computeProbability_[i] = false;
    computeProbabilityD1_[i] = false;
    computeProbabilityD2_[i] = false;
  }
}
 
AbstractSubstitutionProcess::AbstractSubstitutionProcess(const AbstractSubstitutionProcess& asp) :
  pTree_(pTree_->clone()),
  nodeIndex_(asp.nodeIndex_),
  nbClasses_(asp.nbClasses_),
  probabilities_(asp.probabilities_),
  probabilitiesD1_(asp.probabilitiesD1_),
  probabilitiesD2_(asp.probabilitiesD2_),
  computeProbability_(asp.computeProbability_),
  computeProbabilityD1_(asp.computeProbabilityD1_),
  computeProbabilityD2_(asp.computeProbabilityD2_)
{}

AbstractSubstitutionProcess& AbstractSubstitutionProcess::operator=(const AbstractSubstitutionProcess& asp)
{
  pTree_.reset(pTree_->clone());
  nodeIndex_ = asp.nodeIndex_;
  nbClasses_ = asp.nbClasses_;
  probabilities_ = asp.probabilities_;
  probabilitiesD1_ = asp.probabilitiesD1_;
  probabilitiesD2_ = asp.probabilitiesD2_;
  computeProbability_ = asp.computeProbability_;
  computeProbabilityD1_ = asp.computeProbabilityD1_;
  computeProbabilityD2_ = asp.computeProbabilityD2_;
  return *this;
}

size_t AbstractSubstitutionProcess::getNodeIndex_(int nodeId) const throw (NodeNotFoundException)
{
  std::map<int, size_t>::const_iterator it = nodeIndex_.find(nodeId);
  if (it != nodeIndex_.end())
    return it->second;
  else
    throw NodeNotFoundException("AbstractSubstitutionProcess::getNodeIndex(int).", nodeId);
}

size_t AbstractSubstitutionProcess::getModelIndex_(int nodeId, size_t modelClass) const throw (NodeNotFoundException, IndexOutOfBoundsException)
{
  size_t i = getNodeIndex_(nodeId);
  if (modelClass >= nbClasses_)
    throw IndexOutOfBoundsException("AbstractSubstitutionProcess::getModelIndex_().", modelClass, 0, nbClasses_);
  return i * nbClasses_ + modelClass;
}

void AbstractSubstitutionProcess::fireParameterChanged(const ParameterList& pl)
{
  for (size_t i = 0; i < computeProbability_.size(); ++i) {
    computeProbability_[i] = false;
    computeProbabilityD1_[i] = false;
    computeProbabilityD2_[i] = false;
  }
}

