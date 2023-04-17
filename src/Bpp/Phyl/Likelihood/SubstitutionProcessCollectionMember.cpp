//
// File: SubstitutionProcessCollectionMember.cpp
// Authors:
//   Laurent Guéguen
// Created: lundi 1 juillet 2013, ÃÂ  14h 51
//

/*
  Copyright or <A9> or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#include <Bpp/Utils/MapTools.h>

#include "../Model/MixedTransitionModel.h"
#include "SubstitutionProcessCollection.h"
#include "SubstitutionProcessCollectionMember.h"

using namespace bpp;
using namespace std;

std::vector<size_t> SubstitutionProcessCollectionMember::getModelNumbers() const
{
  vector<size_t> vMN;
  for (const auto& it : modelToNodes_)
  {
    vMN.push_back(it.first);
  }

  return vMN;
}

const BranchModelInterface& SubstitutionProcessCollectionMember::model(size_t n) const
{
  return collection().model(n);
}
  
std::shared_ptr<const BranchModelInterface> SubstitutionProcessCollectionMember::getModel(size_t n) const
{
  return collection().getModel(n);
}

std::shared_ptr<BranchModelInterface> SubstitutionProcessCollectionMember::getModel(size_t n)
{
  return collection().getModel(n);
}

std::shared_ptr<const ModelScenario> SubstitutionProcessCollectionMember::getModelScenario() const
{
  return collection().getModelScenario(nPath_);
}

std::shared_ptr<ModelScenario> SubstitutionProcessCollectionMember::getModelScenario()
{
  return collection().getModelScenario(nPath_);
}

std::shared_ptr<const DiscreteDistribution> SubstitutionProcessCollectionMember::getRateDistribution() const
{
  return collection().getRateDistribution(nDist_);
}

std::shared_ptr<DiscreteDistribution> SubstitutionProcessCollectionMember::getRateDistribution()
{
  return collection().getRateDistribution(nDist_);
}

const DiscreteDistribution& SubstitutionProcessCollectionMember::rateDistribution() const
{
  return collection().rateDistribution(nDist_);
}

DiscreteDistribution& SubstitutionProcessCollectionMember::rateDistribution()
{
  return collection().rateDistribution(nDist_);
}

ParameterList SubstitutionProcessCollectionMember::getRateDistributionParameters(bool independent) const
{
  return collection().getRateDistributionParameters(nDist_, independent);
}

ParameterList SubstitutionProcessCollectionMember::getBranchLengthParameters(bool independent) const
{
  if (nTree_!=0)
    return collection().getBranchLengthParameters(nTree_, independent);
  else
    return ParameterList();
}

ParameterList SubstitutionProcessCollectionMember::getRootFrequenciesParameters(bool independent) const
{
  if (!isStationary())
    return collection().getRootFrequenciesParameters(nRoot_, independent);
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
  ParameterList pl = getSubstitutionModelParameters(false);
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


std::shared_ptr<const FrequencySetInterface> SubstitutionProcessCollectionMember::getRootFrequencySet() const
{
  if (isStationary())
    return nullptr;
  else
    return collection().getFrequencySet(nRoot_);
}

std::shared_ptr<FrequencySetInterface> SubstitutionProcessCollectionMember::getRootFrequencySet()
{
  if (isStationary())
    return nullptr;
  else
    return collection().getFrequencySet(nRoot_);
}

const FrequencySetInterface& SubstitutionProcessCollectionMember::rootFrequencySet() const
{
  if (isStationary())
    throw NullPointerException("SubstitutionProcessCollectionMember::rootFrequencySet(). No root frequency parameters, this is a stationary model.");
  else
    return collection().frequencySet(nRoot_);
}

FrequencySetInterface& SubstitutionProcessCollectionMember::rootFrequencySet()
{
  if (isStationary())
    throw NullPointerException("SubstitutionProcessCollectionMember::rootFrequencySet(). No root frequency parameters, this is a stationary model.");
  else
    return collection().frequencySet(nRoot_);
}

const std::vector<double>& SubstitutionProcessCollectionMember::getRootFrequencies() const
{
  auto model = dynamic_pointer_cast<const TransitionModelInterface>(getCollection()->getModel(modelToNodes_.begin()->first));

  if (isStationary() && model)
    return model->getFrequencies();
  else
    return collection().frequencySet(nRoot_).getFrequencies();
}

void SubstitutionProcessCollectionMember::setModelScenario(size_t numPath)
{
  if (!getCollection()->hasModelScenario(numPath))
    throw Exception("SubstitutionProcessCollectionMember::setModelScenario: Collection does not have ModelScenario number");

  auto modelScenario = getCollection()->getModelScenario(numPath);

  // Now check all the models of the path are included in the process.

  auto models = modelScenario->getModels();

  auto modnum = getModelNumbers();

  for (const auto& model:models)
  {
    bool ok = false;
    for (auto num:modnum)
    {
      if (getModel(num) == model)
      {
        ok = true;
        break;
      }
    }

    if (!ok)
      throw Exception("SubstitutionProcessCollectionMember::setModelScenario: Unknown model " + model->getName());
  }

  nPath_ = numPath;
}


const ParametrizablePhyloTree& SubstitutionProcessCollectionMember::parametrizablePhyloTree() const
{
  if (collection().hasTreeNumber(nTree_))
    return collection().tree(nTree_);
  else
    throw Exception("SubstitutionProcessCollectionMember::parametrizablePhyloTree(). No associated tree.");
}

std::shared_ptr<const ParametrizablePhyloTree> SubstitutionProcessCollectionMember::getParametrizablePhyloTree() const
{
  return collection().hasTreeNumber(nTree_) ? collection().getTree(nTree_) : nullptr;
}

ParametrizablePhyloTree& SubstitutionProcessCollectionMember::parametrizablePhyloTree()
{
  if (collection().hasTreeNumber(nTree_))
    return collection().tree(nTree_);
  else
    throw Exception("SubstitutionProcessCollectionMember::parametrizablePhyloTree(). No associated tree.");
}

std::shared_ptr<ParametrizablePhyloTree> SubstitutionProcessCollectionMember::getParametrizablePhyloTree()
{
  return collection().hasTreeNumber(nTree_) ? collection().getTree(nTree_) : nullptr;
}

void SubstitutionProcessCollectionMember::addModel(size_t numModel, const std::vector<unsigned int>& nodesId)
{
  auto& nmod = getCollection()->model(numModel);

  if (modelToNodes_.size() > 0)
  {
    auto& modi = getCollection()->model(modelToNodes_.begin()->first);
    if (nmod.getAlphabet()->getAlphabetType() !=  modi.getAlphabet()->getAlphabetType())
      throw Exception("SubstitutionProcessCollectionMember::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same alphabet.");
    if (nmod.getNumberOfStates() != modi.getNumberOfStates())
      throw Exception("SubstitutionProcessCollectionMember::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same number of states.");
  }
  else if (!isStationary())
  {
    auto& freq = getCollection()->frequencySet(nRoot_);
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

  updateParameters();
}


void SubstitutionProcessCollectionMember::setRootFrequencies(size_t numFreq)
{
  auto& freq = getCollection()->frequencySet(numFreq);
  if (modelToNodes_.size() > 0)
  {
    auto& modi = getCollection()->model(modelToNodes_.begin()->first);

    if (freq.getAlphabet()->getAlphabetType() != modi.getAlphabet()->getAlphabetType())
      throw Exception("SubstitutionProcessCollectionMember::setRootFrequencies. A Frequencies Set cannot be added to a Model Set if it does not have the same alphabet as the models.");
    if (freq.getFrequencies().size() != modi.getNumberOfStates())
      throw Exception("SubstitutionProcessCollectionMember::setRootFrequencies. A Frequencies Set cannot be added to a Model Set if it does not have the same number of states as the models.");
  }

  nRoot_ = numFreq;

  updateParameters();
}

bool SubstitutionProcessCollectionMember::checkOrphanNodes(bool throwEx) const
{
  if (!getParametrizablePhyloTree())
  {
    if (throwEx)
      throw Exception("SubstitutionProcessCollectionMember::checkOrphanNodes(). No Tree");
    
    return true;
  }
  
  vector<unsigned int> ids = getParametrizablePhyloTree()->getAllNodesIndexes();
  unsigned int rootId = getParametrizablePhyloTree()->getNodeIndex(getParametrizablePhyloTree()->getRoot());
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
  if (!getParametrizablePhyloTree())
  {
    if (throwEx)
      throw Exception("SubstitutionProcessCollectionMember::checkUnknownNodes(). No Tree");
    
    return true;
  }
  
  vector<unsigned int> ids = getParametrizablePhyloTree()->getAllNodesIndexes();

  unsigned int rootId = getParametrizablePhyloTree()->getNodeIndex(getParametrizablePhyloTree()->getRoot());

  for (auto& it : modelToNodes_)
  {
    for (auto id : it.second)
    {
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
  for (auto& it : modelToNodes_)
  {
    if (dynamic_pointer_cast<const MixedTransitionModelInterface>(collection().getModel(it.first)))
    {
      return true;
    }
  }
  return false;
}

bool SubstitutionProcessCollectionMember::matchParametersValues(const ParameterList& parameters)
{
  return collection().matchParametersValues(parameters);
}

/**
 * Inheriting from SubstitutionProcess
 */

double SubstitutionProcessCollectionMember::getProbabilityForModel(size_t classIndex) const
{
  if (classIndex >= rateDistribution().getNumberOfCategories())
    throw IndexOutOfBoundsException("SubstitutionProcessCollectionMember::getProbabilityForModel.", classIndex, 0, getRateDistribution()->getNumberOfCategories());
  return rateDistribution().getProbability(classIndex);
}

Vdouble SubstitutionProcessCollectionMember::getClassProbabilities() const
{
  Vdouble vProb;

  for (size_t i = 0; i < rateDistribution().getNumberOfCategories(); i++)
  {
    vProb.push_back(rateDistribution().getProbability(i));
  }

  return vProb;
}

double SubstitutionProcessCollectionMember::getRateForModel(size_t classIndex) const
{
  if (classIndex >= getRateDistribution()->getNumberOfCategories())
    throw IndexOutOfBoundsException("SubstitutionProcessCollectionMember::getRateForModel.", classIndex, 0, getRateDistribution()->getNumberOfCategories());
  return rateDistribution().getCategory(classIndex);
}

