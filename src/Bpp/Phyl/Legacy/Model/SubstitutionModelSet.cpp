// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Utils/MapTools.h>

#include "../../Model/MixedTransitionModel.h"
#include "SubstitutionModelSet.h"

using namespace bpp;
using namespace std;

SubstitutionModelSet::SubstitutionModelSet(const SubstitutionModelSet& set) :
  AbstractParameterAliasable(set),
  alphabet_             (set.alphabet_),
  nbStates_             (set.nbStates_),
  modelSet_(set.modelSet_.size()),
  rootFrequencies_(set.stationarity_ ? nullptr : set.rootFrequencies_->clone()),
  nodeToModel_          (set.nodeToModel_),
  modelToNodes_         (set.modelToNodes_),
  modelParameters_      (set.modelParameters_),
  stationarity_         (set.stationarity_)
{
  // Duplicate all model objects:
  for (size_t i = 0; i < set.modelSet_.size(); ++i)
  {
    modelSet_[i].reset(set.modelSet_[i]->clone());
  }
}

SubstitutionModelSet& SubstitutionModelSet::operator=(const SubstitutionModelSet& set)
{
  AbstractParameterAliasable::operator=(set);
  alphabet_            = set.alphabet_;
  nbStates_            = set.nbStates_;
  nodeToModel_         = set.nodeToModel_;
  modelToNodes_        = set.modelToNodes_;
  modelParameters_     = set.modelParameters_;
  stationarity_        = set.stationarity_;
  if (set.stationarity_)
    rootFrequencies_ = nullptr;
  else
    rootFrequencies_.reset(set.rootFrequencies_->clone());

  // Duplicate all model objects:
  modelSet_.resize(set.modelSet_.size());
  for (size_t i = 0; i < set.modelSet_.size(); ++i)
  {
    modelSet_[i].reset(set.modelSet_[i]->clone());
  }
  return *this;
}

void SubstitutionModelSet::clear()
{
  resetParameters_();
  nbStates_ = 0;
  modelSet_.clear();
  rootFrequencies_.reset();
  nodeToModel_.clear();
  modelParameters_.clear();
  stationarity_ = true;
}

void SubstitutionModelSet::setRootFrequencies(shared_ptr<FrequencySetInterface> rootFreqs)
{
  if (rootFreqs)
  {
    stationarity_ = false;
    rootFrequencies_ = rootFreqs;
    addParameters_(rootFrequencies_->getParameters());
  }
}

std::vector<int> SubstitutionModelSet::getNodesWithParameter(const std::string& name) const
{
  if (!(hasParameter(name)))
    throw ParameterNotFoundException("SubstitutionModelSet::getNodesWithParameter.", name);

  vector<string> nalias = getAlias(name);
  size_t p = name.rfind("_");
  vector<int> inode = getNodesWithModel(TextTools::to<size_t>(name.substr(p + 1, string::npos)) - 1);

  for (size_t i = 0; i < nalias.size(); i++)
  {
    p = nalias[i].rfind("_");
    size_t pos = TextTools::to<size_t>(nalias[i].substr(p + 1, string::npos));
    if (pos > 0)
    {
      vector<int> ni = getNodesWithModel(pos - 1);
      inode.insert(inode.end(), ni.begin(), ni.end());
    }
  }

  return inode;
}

void SubstitutionModelSet::addModel(shared_ptr<TransitionModelInterface> model, const std::vector<int>& nodesId) // , const vector<string>& newParams)
{
  if (model->alphabet().getAlphabetType() != alphabet_->getAlphabetType())
    throw Exception("SubstitutionModelSet::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same alphabet.");
  if (modelSet_.size() > 0 && model->getNumberOfStates() != nbStates_)
    throw Exception("SubstitutionModelSet::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same number of states.");
  if (modelSet_.size() == 0)
    nbStates_ = model->getNumberOfStates();
  modelSet_.push_back(model);
  size_t thisModelIndex = modelSet_.size() - 1;

  // Associate this model to specified nodes:
  for (size_t i = 0; i < nodesId.size(); i++)
  {
    nodeToModel_[nodesId[i]] = thisModelIndex;
    modelToNodes_[thisModelIndex].push_back(nodesId[i]);
  }

  // Associate parameters:
  string pname;

  vector<string> nplm = model->getParameters().getParameterNames();

  modelParameters_.push_back(model->getParameters());

  for (size_t i  = 0; i < nplm.size(); i++)
  {
    pname = nplm[i];
    Parameter* p = new Parameter(model->getParameters().parameter(pname)); // We work with namespaces here, so model->parameter(pname) does not work.
    p->setName(pname + "_" + TextTools::toString(thisModelIndex + 1));
    addParameter_(p);
  }
}

void SubstitutionModelSet::resetModelToNodeIds()
{
  // reset nodeToModel_
  nodeToModel_.clear();
  // reset modelToNodes_
  modelToNodes_.clear();
}

void SubstitutionModelSet::setNodeToModel(size_t modelIndex, int nodeId)
{
  if (modelIndex > modelSet_.size() - 1)
    throw Exception("SubstitutionModelSet::setNodesToModel. There is no Substitution Model of index " + TextTools::toString(modelIndex));

  nodeToModel_[nodeId] = modelIndex;
  modelToNodes_[modelIndex].push_back(nodeId);
}

void SubstitutionModelSet::replaceModel(size_t modelIndex, shared_ptr<TransitionModelInterface> model)
{
  modelSet_[modelIndex] = model;

  // Erase all parameter references to this model

  ParameterList pl = getNodeParameters();

  for (size_t i = pl.size(); i > 0; i--)
  {
    string pn = pl[i - 1].getName();

    size_t pu = pn.rfind("_");
    int nm = TextTools::toInt(pn.substr(pu + 1, string::npos));

    if (nm == (int)modelIndex + 1)
    {
      vector<string> alpn = getAlias(pn);
      for (unsigned j = 0; j < alpn.size(); j++)
      {
        try
        {
          unaliasParameters(alpn[j], pn);
        }
        catch (Exception& e)
        {
          continue;
        }
      }
      deleteParameter_(pn);
    }
  }

  // Associate new parameters
  string pname;

  vector<string> nplm = model->getParameters().getParameterNames();

  for (size_t i  = 0; i < nplm.size(); i++)
  {
    pname = nplm[i];
    Parameter* p = new Parameter(model->getParameters().parameter(pname)); // We work with namespaces here, so model->parameter(pname) does not work.
    p->setName(pname + "_" + TextTools::toString(modelIndex + 1));
    addParameter_(p);
  }

  // update modelParameters_

  modelParameters_[modelIndex].reset();
  modelParameters_[modelIndex] = *model->getParameters().clone();
}

void SubstitutionModelSet::listModelNames(std::ostream& out) const
{
  for (size_t i = 0; i < modelSet_.size(); i++)
  {
    out << "Model " << i + 1 << ": " << modelSet_[i]->getName() << "\t attached to nodes ";
    for (size_t j = 0; j < modelToNodes_[i].size(); j++)
    {
      out << modelToNodes_[i][j];
    }
    out << endl;
  }
}

void SubstitutionModelSet::fireParameterChanged(const ParameterList& parameters)
{
  // Update root frequencies:
  updateRootFrequencies();

  // Then we update all models in the set:
  for (size_t i = 0; i < modelParameters_.size(); i++)
  {
    for (size_t np = 0; np < modelParameters_[i].size(); np++)
    {
      modelParameters_[i][np].setValue(getParameterValue(modelParameters_[i][np].getName() + "_" + TextTools::toString(i + 1)));
    }
    modelSet_[i]->matchParametersValues(modelParameters_[i]);
  }
}

bool SubstitutionModelSet::checkOrphanModels(bool throwEx) const
{
  vector<size_t> index = MapTools::getValues(nodeToModel_);
  for (size_t i = 0; i < modelSet_.size(); i++)
  {
    if (!VectorTools::contains(index, i))
    {
      if (throwEx)
        throw Exception("SubstitutionModelSet::checkOrphanModels(). Model '" + TextTools::toString(i + 1) + "' is associated to no node.");
      return false;
    }
  }
  return true;
}

bool SubstitutionModelSet::checkOrphanNodes(const Tree& tree, bool throwEx) const
{
  vector<int> ids = tree.getNodesId();
  int rootId = tree.getRootId();
  for (size_t i = 0; i < ids.size(); i++)
  {
    if (ids[i] != rootId && nodeToModel_.find(ids[i]) == nodeToModel_.end())
    {
      if (throwEx)
        throw Exception("SubstitutionModelSet::checkOrphanNodes(). Node '" + TextTools::toString(ids[i]) + "' in tree has no model associated.");
      return false;
    }
  }
  return true;
}

bool SubstitutionModelSet::checkUnknownNodes(const Tree& tree, bool throwEx) const
{
  vector<int> ids = tree.getNodesId();
  int id;
  int rootId = tree.getRootId();
  for (size_t i = 0; i < modelToNodes_.size(); i++)
  {
    for (size_t j = 0; j < modelToNodes_[i].size(); j++)
    {
      id = modelToNodes_[i][j];
      if (id == rootId || !VectorTools::contains(ids, id))
      {
        if (throwEx)
          throw Exception("SubstitutionModelSet::checkUnknownNodes(). Node '" + TextTools::toString(id) + "' is not found in tree or is the root node.");
        return false;
      }
    }
  }
  return true;
}

bool SubstitutionModelSet::hasMixedTransitionModel() const
{
  for (size_t i = 0; i < getNumberOfModels(); ++i)
  {
    if ((dynamic_cast<const MixedTransitionModelInterface*>(getModel(i).get()) != nullptr) && (modelToNodes_[i].size() > 1))
      return true;
  }
  return false;
}
