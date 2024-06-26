// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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

std::shared_ptr<const DiscreteDistributionInterface> SubstitutionProcessCollectionMember::getRateDistribution() const
{
  return collection().getRateDistribution(nDist_);
}

std::shared_ptr<DiscreteDistributionInterface> SubstitutionProcessCollectionMember::getRateDistribution()
{
  return collection().getRateDistribution(nDist_);
}

const DiscreteDistributionInterface& SubstitutionProcessCollectionMember::rateDistribution() const
{
  return collection().rateDistribution(nDist_);
}

DiscreteDistributionInterface& SubstitutionProcessCollectionMember::rateDistribution()
{
  return collection().rateDistribution(nDist_);
}

ParameterList SubstitutionProcessCollectionMember::getRateDistributionParameters(bool independent) const
{
  return collection().getRateDistributionParameters(nDist_, independent);
}

ParameterList SubstitutionProcessCollectionMember::getBranchLengthParameters(bool independent) const
{
  if (nTree_ != 0)
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
  std::map<size_t, std::vector<unsigned int>>::const_iterator it;
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
      throw Exception("SubstitutionProcessCollectionMember::setModelScenario: Model " + model->getName() + " in used in scenario" + TextTools::toString(numPath) + " but is unknown from process " + TextTools::toString(nProc_));
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

void SubstitutionProcessCollectionMember::setTreeNumber(size_t nTree, bool check)
{
  if (!collection().hasTreeNumber(nTree))
    throw BadIntException((int)nTree, "SubstitutionProcessCollectionMember::setTreeNumber(). No associated tree.", getAlphabet());

  deleteParameters_(getBranchLengthParameters(true).getParameterNames());

  nTree_ = nTree;

  addParameters_(getBranchLengthParameters(true));

  if (check)
    isFullySetUp();
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
  collection().matchParametersValues(parameters);
  return AbstractSubstitutionProcess::matchParametersValues(parameters);
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
