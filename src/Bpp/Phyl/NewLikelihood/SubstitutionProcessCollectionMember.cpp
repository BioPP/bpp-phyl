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
#include "../Model/MixedSubstitutionModel.h"

#include <Bpp/Utils/MapTools.h>

using namespace bpp;
using namespace std;

SubstitutionProcessCollectionMember::SubstitutionProcessCollectionMember(const SubstitutionProcessCollection* pSubProColl, size_t nTree, size_t nDist) :
  AbstractParameterAliasable(""),
  pSubProColl_(pSubProColl),
  nodeToModel_(),
  modelToNodes_(),
  nTree_(nTree),
  nDist_(nDist),
  stationarity_(true),
  nRoot_(0),
  computingTree_(pSubProColl_->getTree(nTree_)->getTree(),*pSubProColl_->getDistribution(nDist_))
{
  addParameters_(pSubProColl_->getTree(nTree_)->getParameters());
  addParameters_(pSubProColl_->getDistribution(nDist_)->getParameters());
}

SubstitutionProcessCollectionMember::SubstitutionProcessCollectionMember(const SubstitutionProcessCollectionMember& set) :
  AbstractParameterAliasable(set),
  pSubProColl_(set.pSubProColl_),
  nodeToModel_(set.nodeToModel_),
  modelToNodes_(set.modelToNodes_),
  nTree_(set.nTree_),
  nDist_(set.nDist_),
  stationarity_(set.stationarity_),
  nRoot_(set.nRoot_),
  computingTree_(set.computingTree_)
{
}

SubstitutionProcessCollectionMember& SubstitutionProcessCollectionMember::operator=(const SubstitutionProcessCollectionMember& set)
{
  AbstractParameterAliasable::operator=(set);
  pSubProColl_ = set.pSubProColl_;
  nodeToModel_ = set.nodeToModel_;
  modelToNodes_ = set.modelToNodes_;
  nTree_ = set.nTree_;
  nDist_ = set.nDist_;
  stationarity_ = set.stationarity_;
  nRoot_ = set.nRoot_;
  computingTree_ = set.computingTree_;

  return *this;
}

void SubstitutionProcessCollectionMember::clear()
{
  nodeToModel_.clear();
  modelToNodes_.clear();

  stationarity_=true;
}


const Alphabet* SubstitutionProcessCollectionMember::getAlphabet() const
{
  return pSubProColl_->getModel(modelToNodes_.begin()->first)->getAlphabet();
}

const SubstitutionModel* SubstitutionProcessCollectionMember::getModel(size_t i) const
{
  return pSubProColl_->getModel(i);
}


const DiscreteDistribution* SubstitutionProcessCollectionMember::getDistribution() const
{
  return pSubProColl_->getDistribution(nDist_);
}

const FrequenciesSet* SubstitutionProcessCollectionMember::getRootFrequenciesSet() const
{
  if (stationarity_)
    return 0;
  else
    return pSubProColl_->getFrequencies(nRoot_);
}

const std::vector<double>& SubstitutionProcessCollectionMember::getRootFrequencies() const
{
  if (stationarity_)
    return (pSubProColl_->getModel(modelToNodes_.begin()->first))->getFrequencies();
  else
    return (pSubProColl_->getFrequencies(nRoot_))->getFrequencies();
}

const TreeTemplate<Node>& SubstitutionProcessCollectionMember::getTree() const
{
  return pSubProColl_->getTree(nTree_)->getTree();
}

const ParametrizableTree& SubstitutionProcessCollectionMember::getParametrizableTree() const
{
  return *pSubProColl_->getTree(nTree_);
}

size_t SubstitutionProcessCollectionMember::getNumberOfClasses() const
{
  return pSubProColl_->getDistribution(nDist_)->getNumberOfCategories();
}


void SubstitutionProcessCollectionMember::addModel(size_t numModel, const std::vector<int>& nodesId)
{
  const SubstitutionModel* nmod=pSubProColl_->getModel(numModel);
  
  if (modelToNodes_.size()>0){
    const SubstitutionModel* pmodi=pSubProColl_->getModel(modelToNodes_.begin()->first);
    if (nmod->getAlphabet()->getAlphabetType() !=  pmodi->getAlphabet()->getAlphabetType())
      throw Exception("SubstitutionProcessCollectionMember::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same alphabet.");
    if (nmod->getNumberOfStates() != pmodi->getNumberOfStates())
      throw Exception("SubstitutionProcessCollectionMember::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same number of states.");
  }
  else if (!stationarity_) {
    const FrequenciesSet* pfreq=pSubProColl_->getFrequencies(nRoot_);
    if (pfreq->getAlphabet()->getAlphabetType() != nmod->getAlphabet()->getAlphabetType())
      throw Exception("SubstitutionProcessCollectionMember::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same alphabet as the root frequencies.");
    if (pfreq->getFrequencies().size() != nmod->getNumberOfStates())
      throw Exception("SubstitutionProcessCollectionMember::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same number of states as the root frequencies.");
  }
  
  // Associate this model to specified nodes:
  for (size_t i = 0; i < nodesId.size(); i++)
    {
      nodeToModel_[nodesId[i]] = numModel;
      modelToNodes_[numModel].push_back(nodesId[i]);
    }


  computingTree_.addModel(nmod, nodesId);
}

void SubstitutionProcessCollectionMember::setRootFrequencies(size_t numFreq)
{
  const FrequenciesSet* pfreq=pSubProColl_->getFrequencies(numFreq);
  if (modelToNodes_.size()>0){
    const SubstitutionModel* pmodi=pSubProColl_->getModel(modelToNodes_.begin()->first);

    if (pfreq->getAlphabet()->getAlphabetType() != pmodi->getAlphabet()->getAlphabetType())
      throw Exception("SubstitutionProcessCollectionMember::addRootFrequencies. A Frequencies Set cannot be added to a Model Set if it does not have the same alphabet as the models.");
    if (pfreq->getFrequencies().size() != pmodi->getNumberOfStates())
      throw Exception("SubstitutionProcessCollectionMember::addRootFrequencies. A Frequencies Set cannot be added to a Model Set if it does not have the same number of states as the models.");
  }
  
  stationarity_=false;
  nRoot_=numFreq;
}

bool SubstitutionProcessCollectionMember::checkOrphanNodes(bool throwEx) const
  throw (Exception)
{
  vector<int> ids = getTree().getNodesId();
  int rootId = getTree().getRootId();
  for (size_t i = 0; i < ids.size(); i++)
    {
      if (ids[i] != rootId && nodeToModel_.find(ids[i]) == nodeToModel_.end())
        {
          if (throwEx) throw Exception("SubstitutionProcessCollectionMember::checkOrphanNodes(). Node '" + TextTools::toString(ids[i]) + "' in tree has no model associated.");
          return false;
        }
    }
  return true;
}

bool SubstitutionProcessCollectionMember::checkUnknownNodes(bool throwEx) const
  throw (Exception)
{
  vector<int> ids = getTree().getNodesId();
  int id;
  int rootId = getTree().getRootId();
  for (size_t i = 0; i < modelToNodes_.size(); i++)
    {
      for (size_t j = 0; j < modelToNodes_.at(i).size(); j++)
        {
          id = modelToNodes_.at(i)[j];
          if (id == rootId || !VectorTools::contains(ids, id))
            {
              if (throwEx) throw Exception("SubstitutionProcessCollectionMember::checkUnknownNodes(). Node '" + TextTools::toString(id) + "' is not found in tree or is the root node.");
              return false;
            }
        }
    }
  return true;
}

bool SubstitutionProcessCollectionMember::hasMixedSubstitutionModel() const
{
  std::map<size_t, std::vector<int> >::const_iterator it;
  for (it= modelToNodes_.begin(); it != modelToNodes_.end() ; it++)
    {
      if (dynamic_cast<const MixedSubstitutionModel*>(pSubProColl_->getModel(it->first)) != NULL)
        return true;
    }
  return false;
}

/*
 * Inheriting from SubstitutionProcess
 */
  
bool SubstitutionProcessCollectionMember::isCompatibleWith(const SiteContainer& data) const
{
  if (modelToNodes_.size() > 0) 
    return data.getAlphabet()->getAlphabetType() == pSubProColl_->getModel(modelToNodes_.begin()->first)->getAlphabet()->getAlphabetType();
  else
    return true;
}


bool SubstitutionProcessCollectionMember::hasTransitionProbabilitiesParameter(const std::string& name) const
{
  return false;
}

const SubstitutionModel* SubstitutionProcessCollectionMember::getModelForNode(int nodeId) const throw (Exception)
{
  std::map<int, size_t>::const_iterator i = nodeToModel_.find(nodeId);
  if (i == nodeToModel_.end())
    throw Exception("SubstitutionProcessCollectionMember::getModelForNode(). No model associated to node with id " + TextTools::toString(nodeId));
  return getModel(i->second);
}

size_t SubstitutionProcessCollectionMember::getNumberOfStates() const
{
  if (modelToNodes_.size()==0)
    return 0;
  else 
    return getModel(modelToNodes_.begin()->first)->getNumberOfStates();
}

const SubstitutionModel& SubstitutionProcessCollectionMember::getSubstitutionModel(int nodeId, size_t classIndex) const
{
  return *getModel(nodeToModel_.at(nodeId));
}

double SubstitutionProcessCollectionMember::getInitValue(size_t i, int state) const throw (BadIntException)
{
  if (modelToNodes_.size()==0)
    throw Exception("SubstitutionProcessCollectionMember::getInitValue : no model associated");
  else
    return getModel(modelToNodes_.begin()->first)->getInitValue(i,state);
}

double SubstitutionProcessCollectionMember::getProbabilityForModel(size_t classIndex) const {
  if (classIndex >= getDistribution()->getNumberOfCategories())
    throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::getProbabilityForModel.", classIndex, 0, getDistribution()->getNumberOfCategories());
  return getDistribution()->getProbability(classIndex);
}

