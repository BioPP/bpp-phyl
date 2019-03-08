//
// File: SubstitutionProcessCollectionMember.cpp
// Created by: Laurent Guéguen
// Created on: lundi 1 juillet 2013, à 14h 51
//

/*
   Copyright or <A9> or Copr. CNRS, (November 16, 2004)

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

#include "SubstitutionProcessCollectionMember.h"

#include "SubstitutionProcessCollection.h"
#include "../Model/MixedTransitionModel.h"

#include <Bpp/Utils/MapTools.h>

using namespace bpp;
using namespace std;

SubstitutionProcessCollectionMember::SubstitutionProcessCollectionMember( SubstitutionProcessCollection* pSubProColl, size_t nProc, size_t nTree, size_t nDist) :
  AbstractParameterAliasable(""),
  pSubProColl_(pSubProColl),
  nProc_(nProc),
  nodeToModel_(),
  modelToNodes_(),
  nTree_(nTree),
  nDist_(nDist),
  stationarity_(true),
  nRoot_(0),
  computingTree_(new ComputingTree(pSubProColl_, nTree_, nDist_))
{
  updateParameters();
}


SubstitutionProcessCollectionMember::SubstitutionProcessCollectionMember(const SubstitutionProcessCollectionMember& set) :
  AbstractParameterAliasable(set),
  pSubProColl_(set.pSubProColl_),
  nProc_(set.nProc_),
  nodeToModel_(set.nodeToModel_),
  modelToNodes_(set.modelToNodes_),
  nTree_(set.nTree_),
  nDist_(set.nDist_),
  stationarity_(set.stationarity_),
  nRoot_(set.nRoot_),
  computingTree_(new ComputingTree(pSubProColl_, nTree_, nDist_))
{}

SubstitutionProcessCollectionMember& SubstitutionProcessCollectionMember::operator=(const SubstitutionProcessCollectionMember& set)
{
  AbstractParameterAliasable::operator=(set);
  pSubProColl_ = set.pSubProColl_;
  nProc_ = set.nProc_;

  nodeToModel_ = set.nodeToModel_;
  modelToNodes_ = set.modelToNodes_;
  nTree_ = set.nTree_;
  nDist_ = set.nDist_;
  stationarity_ = set.stationarity_;
  nRoot_ = set.nRoot_;

  computingTree_.release();

  computingTree_ = unique_ptr<ComputingTree>(new ComputingTree(pSubProColl_, nTree_, nDist_));

  return *this;
}

void SubstitutionProcessCollectionMember::clear()
{
  nodeToModel_.clear();
  modelToNodes_.clear();

  stationarity_ = true;
}


inline const Alphabet* SubstitutionProcessCollectionMember::getAlphabet() const
{
  return (getCollection()->getModel(modelToNodes_.begin()->first))->getAlphabet();
}

inline const TransitionModel* SubstitutionProcessCollectionMember::getModel(size_t i) const
{
  return getCollection()->getModel(i);
}

inline bool SubstitutionProcessCollectionMember::matchParametersValues(const ParameterList& parameters)
{
  return getCollection()->matchParametersValues(parameters);
}

std::vector<size_t> SubstitutionProcessCollectionMember::getModelNumbers() const
{
  vector<size_t> vMN;
  std::map<size_t, std::vector<unsigned int> >::const_iterator it;
  for (it = modelToNodes_.begin(); it != modelToNodes_.end(); it++)
  {
    vMN.push_back(it->first);
  }

  return vMN;
}

inline const DiscreteDistribution* SubstitutionProcessCollectionMember::getRateDistribution() const
{
  return &getCollection()->getRateDistribution(nDist_);
}

ParameterList SubstitutionProcessCollectionMember::getRateDistributionParameters(bool independent) const
{
  return getCollection()->getRateDistributionParameters(nDist_, independent);
}

ParameterList SubstitutionProcessCollectionMember::getBranchLengthParameters(bool independent) const
{
  return getCollection()->getBranchLengthParameters(nTree_, independent);
}

inline bool SubstitutionProcessCollectionMember::hasBranchLengthParameter(const std::string& name) const
{
  return getCollection()->hasBranchLengthParameter(name);
}

ParameterList SubstitutionProcessCollectionMember::getRootFrequenciesParameters(bool independent) const
{
  if (!stationarity_)
    return getCollection()->getRootFrequenciesParameters(nRoot_, independent);
  else
    return ParameterList();
}

void SubstitutionProcessCollectionMember::updateParameters()
{
  resetParameters_();
  addParameters_(getSubstitutionModelParameters(true));
  addParameters_(getRootFrequenciesParameters(true));
  addParameters_(getRateDistributionParameters(true));
  addParameters_(getBranchLengthParameters(true));
}

ParameterList SubstitutionProcessCollectionMember::getNonDerivableParameters() const
{
  ParameterList pl=getSubstitutionModelParameters(false);
  pl.includeParameters(getRootFrequenciesParameters(false));
  pl.includeParameters(getRateDistributionParameters(false));

  pl.includeParameters(getCollection()->getAliasedParameters(pl));
  pl.includeParameters(getCollection()->getFromParameters(pl));

  return pl;
}


ParameterList SubstitutionProcessCollectionMember::getSubstitutionModelParameters(bool independent) const
{
  ParameterList pl;

  // Then we update all models in the set:
  std::map<size_t, std::vector<unsigned int> >::const_iterator it;
  for (it = modelToNodes_.begin(); it != modelToNodes_.end(); it++)
  {
    pl.includeParameters(getCollection()->getSubstitutionModelParameters(it->first, independent));
  }

  return pl;
}


inline const FrequenciesSet* SubstitutionProcessCollectionMember::getRootFrequenciesSet() const
{
  if (stationarity_)
    return 0;
  else
    return &getCollection()->getFrequencies(nRoot_);
}

inline const std::vector<double>& SubstitutionProcessCollectionMember::getRootFrequencies() const
{
  if (stationarity_)
    return (getCollection()->getModel(modelToNodes_.begin()->first))->getFrequencies();
  else
    return (getCollection()->getFrequencies(nRoot_)).getFrequencies();
}

inline const ParametrizablePhyloTree& SubstitutionProcessCollectionMember::getParametrizablePhyloTree() const
{
  return getCollection()->getTree(nTree_);
}

inline size_t SubstitutionProcessCollectionMember::getNumberOfClasses() const
{
  return getCollection()->getRateDistribution(nDist_).getNumberOfCategories();
}


void SubstitutionProcessCollectionMember::addModel(size_t numModel, const std::vector<unsigned int>& nodesId)
{
  const TransitionModel& nmod = *getCollection()->getModel(numModel);

  if (modelToNodes_.size() > 0)
  {
    const TransitionModel& modi = *getCollection()->getModel(modelToNodes_.begin()->first);
    if (nmod.getAlphabet()->getAlphabetType() !=  modi.getAlphabet()->getAlphabetType())
      throw Exception("SubstitutionProcessCollectionMember::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same alphabet.");
    if (nmod.getNumberOfStates() != modi.getNumberOfStates())
      throw Exception("SubstitutionProcessCollectionMember::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same number of states.");
  }
  else if (!stationarity_)
  {
    const FrequenciesSet& freq = getCollection()->getFrequencies(nRoot_);
    if (freq.getAlphabet()->getAlphabetType() != nmod.getAlphabet()->getAlphabetType())
      throw Exception("SubstitutionProcessCollectionMember::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same alphabet as the root frequencies.");
    if (freq.getFrequencies().size() != nmod.getNumberOfStates())
      throw Exception("SubstitutionProcessCollectionMember::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same number of states as the root frequencies.");
  }

  // Associate this model to specified nodes:
  for (size_t i = 0; i < nodesId.size(); i++)
  {
    nodeToModel_[nodesId[i]] = numModel;
    modelToNodes_[numModel].push_back(nodesId[i]);
  }

  computingTree_->addModel(&nmod, nodesId);

  updateParameters();
}


void SubstitutionProcessCollectionMember::setRootFrequencies(size_t numFreq)
{
  const FrequenciesSet& freq = getCollection()->getFrequencies(numFreq);
  if (modelToNodes_.size() > 0)
  {
    const TransitionModel& modi = *getCollection()->getModel(modelToNodes_.begin()->first);

    if (freq.getAlphabet()->getAlphabetType() != modi.getAlphabet()->getAlphabetType())
      throw Exception("SubstitutionProcessCollectionMember::setRootFrequencies. A Frequencies Set cannot be added to a Model Set if it does not have the same alphabet as the models.");
    if (freq.getFrequencies().size() != modi.getNumberOfStates())
      throw Exception("SubstitutionProcessCollectionMember::setRootFrequencies. A Frequencies Set cannot be added to a Model Set if it does not have the same number of states as the models.");
  }

  stationarity_ = false;
  nRoot_ = numFreq;

  updateParameters();
}

bool SubstitutionProcessCollectionMember::checkOrphanNodes(bool throwEx) const
{
  vector<unsigned int> ids = getParametrizablePhyloTree().getAllNodesIndexes();
  unsigned int rootId = getParametrizablePhyloTree().getNodeIndex(getParametrizablePhyloTree().getRoot());
  for (size_t i = 0; i < ids.size(); i++)
  {
    if (ids[i] != rootId && nodeToModel_.find(ids[i]) == nodeToModel_.end())
    {
      if (throwEx)
        throw Exception("SubstitutionProcessCollectionMember::checkOrphanNodes(). Node '" + TextTools::toString(ids[i]) + "' in tree has no model associated.");
      return false;
    }
  }
  return true;
}

bool SubstitutionProcessCollectionMember::checkUnknownNodes(bool throwEx) const
{
  vector<unsigned int> ids = getParametrizablePhyloTree().getAllNodesIndexes();

  unsigned int id;
  unsigned int rootId = getParametrizablePhyloTree().getNodeIndex(getParametrizablePhyloTree().getRoot());
  
  std::map<size_t, std::vector<unsigned int> >::const_iterator it;

  for (it = modelToNodes_.begin(); it != modelToNodes_.end(); it++)
  {
    for (size_t j = 0; j < it->second.size(); j++)
    {
      id = it->second[j];
      if (id == rootId || !VectorTools::contains(ids, id))
      {
        if (throwEx)
          throw Exception("SubstitutionProcessCollectionMember::checkUnknownNodes(). Node '" + TextTools::toString(id) + "' is not found in tree or is the root node.");
        return false;
      }
    }
  }
  return true;
}

bool SubstitutionProcessCollectionMember::hasMixedTransitionModel() const
{
  std::map<size_t, std::vector<unsigned int> >::const_iterator it;
  for (it = modelToNodes_.begin(); it != modelToNodes_.end(); it++)
  {
    if (dynamic_cast<const MixedTransitionModel*>(getCollection()->getModel(it->first)) != NULL)
      return true;
  }
  return false;
}

/*
 * Inheriting from SubstitutionProcess
 */
bool SubstitutionProcessCollectionMember::isCompatibleWith(const AlignedValuesContainer& data) const
{
  if (modelToNodes_.size() > 0)
    return data.getAlphabet()->getAlphabetType() == getCollection()->getModel(modelToNodes_.begin()->first)->getAlphabet()->getAlphabetType();
  else
    return true;
}


inline const TransitionModel* SubstitutionProcessCollectionMember::getModelForNode(unsigned int nodeId) const 
{
  std::map<unsigned int, size_t>::const_iterator i = nodeToModel_.find(nodeId);
  if (i == nodeToModel_.end())
    throw Exception("SubstitutionProcessCollectionMember::getModelForNode(). No model associated to node with id " + TextTools::toString(nodeId));
  return getModel(i->second);
}

inline size_t SubstitutionProcessCollectionMember::getNumberOfStates() const
{
  if (modelToNodes_.size() == 0)
    return 0;
  else
    return getModel(modelToNodes_.begin()->first)->getNumberOfStates();
}

inline const TransitionModel* SubstitutionProcessCollectionMember::getModel(unsigned int nodeId, size_t classIndex) const
{
  return getModel(nodeToModel_.at(nodeId));
}

inline double SubstitutionProcessCollectionMember::getInitValue(size_t i, int state) const 
{
  if (modelToNodes_.size() == 0)
    throw Exception("SubstitutionProcessCollectionMember::getInitValue : no model associated");
  else
    return getModel(modelToNodes_.begin()->first)->getInitValue(i, state);
}

inline double SubstitutionProcessCollectionMember::getProbabilityForModel(size_t classIndex) const
{
  if (classIndex >= getRateDistribution()->getNumberOfCategories())
    throw IndexOutOfBoundsException("SubstitutionProcessCollectionMember::getProbabilityForModel.", classIndex, 0, getRateDistribution()->getNumberOfCategories());
  return getRateDistribution()->getProbability(classIndex);
}

inline Vdouble SubstitutionProcessCollectionMember::getClassProbabilities() const
{
  Vdouble vProb;

  for (size_t i = 0; i < getRateDistribution()->getNumberOfCategories(); i++)
  {
    vProb.push_back(getRateDistribution()->getProbability(i));
  }

  return vProb;
}

inline double SubstitutionProcessCollectionMember::getRateForModel(size_t classIndex) const
{
  if (classIndex >= getRateDistribution()->getNumberOfCategories())
    throw IndexOutOfBoundsException("SubstitutionProcessCollectionMember::getRateForModel.", classIndex, 0, getRateDistribution()->getNumberOfCategories());
  return getRateDistribution()->getCategory(classIndex);
}
