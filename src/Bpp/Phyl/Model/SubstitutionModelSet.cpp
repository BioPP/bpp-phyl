//
// File: SubstitutionModelSet.cpp
// Created by: Bastien Boussau
//             Julien Dutheil
// Created on: Tue Aug 21 2007
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

#include "SubstitutionModelSet.h"
#include "MixedSubstitutionModel.h"

#include <Bpp/Utils/MapTools.h>

using namespace bpp;
using namespace std;

SubstitutionModelSet::SubstitutionModelSet(const SubstitutionModelSet& set) :
  AbstractParameterAliasable(set),
  alphabet_             (set.alphabet_),
  nbStates_             (set.nbStates_),
  modelSet_(set.modelSet_.size()),
  rootFrequencies_(set.stationarity_ ? 0 : dynamic_cast<FrequenciesSet*>(set.rootFrequencies_->clone())),
  nodeToModel_          (set.nodeToModel_),
  modelToNodes_         (set.modelToNodes_),
  paramToModels_        (set.paramToModels_),
  paramNamesCount_      (set.paramNamesCount_),
  modelParameterNames_  (set.modelParameterNames_),
  modelParameters_      (set.modelParameters_),
  stationarity_         (set.stationarity_)
{
  // Duplicate all model objects:
  for (size_t i = 0; i < set.modelSet_.size(); i++)
    {
      modelSet_[i] = dynamic_cast<SubstitutionModel*>(set.modelSet_[i]->clone());
    }
}

SubstitutionModelSet& SubstitutionModelSet::operator=(const SubstitutionModelSet& set)
{
  AbstractParameterAliasable::operator=(set);
  alphabet_            = set.alphabet_;
  nbStates_            = set.nbStates_;
  nodeToModel_         = set.nodeToModel_;
  modelToNodes_        = set.modelToNodes_;
  paramToModels_       = set.paramToModels_;
  paramNamesCount_     = set.paramNamesCount_;
  modelParameterNames_ = set.modelParameterNames_;
  modelParameters_     = set.modelParameters_;
  stationarity_        = set.stationarity_;
  if (set.stationarity_)
    rootFrequencies_.reset(0);
  else
    rootFrequencies_.reset(dynamic_cast<FrequenciesSet*>(set.rootFrequencies_->clone()));

  // Duplicate all model objects:
  modelSet_.resize(set.modelSet_.size());
  for (size_t i = 0; i < set.modelSet_.size(); i++)
    {
      modelSet_[i] = dynamic_cast<SubstitutionModel*>(set.modelSet_[i]->clone());
    }
  return *this;
}

void SubstitutionModelSet::clear()
{
  resetParameters_();
  nbStates_=0;
  for (size_t i = 0; i < modelSet_.size(); i++)
    {
      delete modelSet_[i];
    }
  modelSet_.clear();
  rootFrequencies_.reset();
  nodeToModel_.clear();
  paramToModels_.clear();
  paramNamesCount_.clear();
  modelParameterNames_.clear();
  modelParameters_.clear();
  stationarity_=true;

}

void SubstitutionModelSet::setRootFrequencies(FrequenciesSet* rootFreqs)
{
  if (rootFreqs){
    stationarity_=false;
    rootFrequencies_.reset(rootFreqs);
    addParameters_(rootFrequencies_->getParameters());
  }
}

std::vector<int> SubstitutionModelSet::getNodesWithParameter(const std::string& name) const
  throw (ParameterNotFoundException)
{
  vector<int> ids;
  size_t offset = stationarity_ ? 0 : rootFrequencies_->getNumberOfParameters();
  for (size_t i = 0; i < paramToModels_.size(); i++)
    {
      if (getParameter_(offset + i).getName() == name)
        {
          for (size_t j = 0; j < paramToModels_[i].size(); j++)
            {
              vector<int> v = modelToNodes_[paramToModels_[i][j]];
              VectorTools::append(ids, v);
            }
          return ids;
        }
    }
  throw ParameterNotFoundException("SubstitutionModelSet::getNodesWithParameter.", name);
}

vector<size_t> SubstitutionModelSet::getModelsWithParameter(const std::string& name) const
  throw (ParameterNotFoundException)
{
  vector<size_t> indices;
  size_t offset = stationarity_ ? 0 : rootFrequencies_->getNumberOfParameters();
  for (size_t i = 0; i < paramToModels_.size(); i++)
    {
      if (getParameter_(offset + i).getName() == name)
        return paramToModels_[i];
    }
  throw ParameterNotFoundException("SubstitutionModelSet::getModelsWithParameter.", name);
}

void SubstitutionModelSet::addModel(SubstitutionModel* model, const std::vector<int>& nodesId, const vector<string>& newParams) throw (Exception)
{
  if (model->getAlphabet()->getAlphabetType() != alphabet_->getAlphabetType())
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
  modelParameters_.push_back(ParameterList());
  for (size_t i  = 0; i < newParams.size(); i++)
    {
      pname = newParams[i];
      modelParameterNames_.push_back(pname);
      vector<size_t> modelsIndex(1, thisModelIndex);
      paramToModels_.push_back(modelsIndex);
      Parameter* p = new Parameter(model->getParameters().getParameter(pname)); // We work with namespaces here, so model->getParameter(pname) does not work.
      modelParameters_[thisModelIndex].addParameter(p->clone());
      p->setName(pname + "_" + TextTools::toString(++paramNamesCount_[pname])); // Change name to unique name in case of shared parameters.
      addParameter_(p);
    }
}

void SubstitutionModelSet::setModel(SubstitutionModel* model, size_t modelIndex) throw (Exception, IndexOutOfBoundsException)
{
  if (model->getAlphabet()->getAlphabetType() != alphabet_->getAlphabetType())
    throw Exception("SubstitutionModelSet::setModel. A Substitution Model cannot be added to a Model Set if it does not have the same alphabet");
  if (modelIndex >= modelSet_.size())
    throw IndexOutOfBoundsException("SubstitutionModelSet::setModel.", modelIndex, 0, modelSet_.size());
  delete modelSet_[modelIndex];
  modelSet_[modelIndex] = model;
}

void SubstitutionModelSet::removeModel(size_t modelIndex) throw (Exception)
{
  modelSet_.erase(modelSet_.begin() + modelIndex);
  // Erase all parameter references to this model and translate other indices...
  for (size_t i = 0; i < paramToModels_.size(); i++)
    {
      for (size_t j = paramToModels_[i].size(); j > 0; j--)
        {
          if (paramToModels_[i][j - 1] == modelIndex)
            {
              paramToModels_[i].erase(paramToModels_[i].begin() + j - 1);
            }
          else if (paramToModels_[i][j - 1] > modelIndex)
            paramToModels_[i][j - 1]--;  // Correct indice due to removal!
        }
    }
  checkOrphanParameters(true);
}

ParameterList SubstitutionModelSet::getModelParameters(size_t modelIndex) const
{
  ParameterList pl;
  size_t offset = stationarity_ ? 0 : rootFrequencies_->getNumberOfParameters(); // Root frequencies are the first parameters! We should ignore them here.
  for (size_t i = 0; i < modelParameterNames_.size(); i++)
    {
      // Check associations:
      const vector<size_t>* modelIndexes = &paramToModels_[i];
      for (size_t j = 0; j < modelIndexes->size(); j++)
        {
          if ((*modelIndexes)[j] == modelIndex)
            {
              pl.addParameter(getParameter_(offset + i));
              break;
            }
        }
    }
  return pl;
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

void SubstitutionModelSet::setParameterToModel(size_t parameterIndex, size_t modelIndex) throw (IndexOutOfBoundsException)
{
  size_t offset = stationarity_ ? 0 : rootFrequencies_->getNumberOfParameters();
  if (parameterIndex < offset)
    throw IndexOutOfBoundsException("SubstitutionModelSet::setParameterToModel. Can't assign a root frequency parameter to a branch model.", parameterIndex, offset - 1, paramToModels_.size() + offset - 1);
  if (parameterIndex >= paramToModels_.size() + offset)
    throw IndexOutOfBoundsException("SubstitutionModelSet::setParameterToModel.", parameterIndex, offset - 1, paramToModels_.size() + offset - 1);
  if (modelIndex >= modelSet_.size())
    throw IndexOutOfBoundsException("SubstitutionModelSet::setParameterToModel.", modelIndex, 0, modelSet_.size() - 1);
  if (VectorTools::contains(paramToModels_[parameterIndex - offset], modelIndex))
    throw Exception("SubstitutionModelSet::setParameterToModel: parameter " + getParameter_(parameterIndex - offset).getName() + " already set to model " + TextTools::toString(modelIndex) + ".");
  paramToModels_[parameterIndex - offset].push_back(modelIndex);
  modelParameters_[modelIndex].addParameter(
                                            modelSet_[modelIndex]->getParameters().getParameter(modelParameterNames_[parameterIndex - offset]));
  //Update value of modified model:
  modelSet_[modelIndex]->matchParametersValues(getParameters().subList(parameterIndex - offset));
}

void SubstitutionModelSet::unsetParameterToModel(size_t parameterIndex, size_t modelIndex) throw (IndexOutOfBoundsException, Exception)
{
  size_t offset = stationarity_ ? 0 : rootFrequencies_->getNumberOfParameters();
  if (parameterIndex < offset)
    throw IndexOutOfBoundsException("SubstitutionModelSet::unsetParameterToModel. Can't unset a root frequency parameter.", parameterIndex, offset - 1, paramToModels_.size() + offset - 1);
  if (parameterIndex >= paramToModels_.size() + offset)
    throw IndexOutOfBoundsException("SubstitutionModelSet::unsetParameterToModel.", parameterIndex, offset - 1, paramToModels_.size() + offset - 1);
  if (modelIndex >= modelSet_.size())
    throw IndexOutOfBoundsException("SubstitutionModelSet::setParameterToModel.", modelIndex, 0, modelSet_.size() - 1);
  if (!VectorTools::contains(paramToModels_[parameterIndex - offset], modelIndex))
    throw Exception("SubstitutionModelSet::unsetParameterToModel: parameter " + getParameter_(parameterIndex).getName() + " is not currently set to model " + TextTools::toString(modelIndex) + ".");
  remove(paramToModels_[parameterIndex - offset].begin(), paramToModels_[parameterIndex - offset].end(), modelIndex);
  modelParameters_[modelIndex].deleteParameter(modelParameterNames_[parameterIndex - offset]);
  checkOrphanModels(true);
  checkOrphanParameters(true);
}

void SubstitutionModelSet::addParameter(const Parameter& parameter, const vector<int>& nodesId) throw (Exception)
{
  modelParameterNames_.push_back(parameter.getName());
  Parameter* p = parameter.clone();
  p->setName(p->getName() + "_" + TextTools::toString(++paramNamesCount_[p->getName()]));
  addParameter_(p);
  // Build model indexes:
  vector<size_t> modelIndexes(nodesId.size());
  for (size_t i = 0; i < nodesId.size(); ++i)
  {
    if (nodeToModel_.find(nodesId[i]) == nodeToModel_.end())
      throw Exception("SubstitutionModelSet::addParameter. This node has no associated model: " + TextTools::toString(nodesId[i]));
    size_t pos = nodeToModel_[nodesId[i]];
    modelParameters_[pos].addParameter(parameter);
    modelIndexes[i] = pos;
  }
  paramToModels_.push_back(modelIndexes);
  //Update model values:
  fireParameterChanged(getParameters().subList(p->getName()));
}

void SubstitutionModelSet::addParameters(const ParameterList& parameters, const vector<int>& nodesId) throw (Exception)
{
  for (size_t i = 0; i < parameters.size(); i++)
    {
      modelParameterNames_.push_back(parameters[i].getName());
    }
  ParameterList pl(parameters);
  for (size_t i = 0; i < pl.size(); i++)
    {
      pl[i].setName(pl[i].getName() + "_" + TextTools::toString(++paramNamesCount_[pl[i].getName()]));
    }
  addParameters_(pl);
  // Build model indexes:
  vector<size_t> modelIndexes(nodesId.size());
  map<size_t, size_t> counts; //Check is a model is affected to several nodes.
  for (size_t i = 0; i < nodesId.size(); i++)
    {
      if (nodeToModel_.find(nodesId[i]) == nodeToModel_.end())
        throw Exception("SubstitutionModelSet::addParameters. This node has no associated model: " + TextTools::toString(nodesId[i]));
      size_t pos = nodeToModel_[nodesId[i]];
      size_t count = counts[pos]++;
      if (count == 0)
        modelParameters_[pos].addParameters(parameters);
      modelIndexes[i] = pos;
    }
  for (size_t i = 0; i < pl.size(); i++)
    {
      paramToModels_.push_back(modelIndexes);
    }
  //Update model values:
  fireParameterChanged(pl);
}

void SubstitutionModelSet::removeParameter(const string& name) throw (ParameterNotFoundException)
{
  size_t offset = stationarity_ ? 0 : rootFrequencies_->getNumberOfParameters();
  for (size_t i = 0; i < modelParameterNames_.size(); i++)
    {
      if (getParameter_(offset + i).getName() == name)
        {
          vector<int> nodesId = getNodesWithParameter(name);
          for (size_t j = 0; j < nodesId.size(); j++)
            {
              size_t pos = nodeToModel_[nodesId[j]];
              string tmp = modelParameterNames_[i];
              if (modelParameters_[pos].hasParameter(tmp))
                modelParameters_[pos].deleteParameter(tmp);
            }
          paramToModels_.erase(paramToModels_.begin() + i);
          modelParameterNames_.erase(modelParameterNames_.begin() + i);
          deleteParameter_(offset + i);
          return;
        }
    }
  throw ParameterNotFoundException("SubstitutionModelSet::removeParameter.", name);
}

void SubstitutionModelSet::fireParameterChanged(const ParameterList& parameters)
{
  // For now, we actualize all parameters... we'll optimize later!
  // Update root frequencies:
  updateRootFrequencies();

  // First we actualize the modelParameters_ array:
  size_t offset = stationarity_ ? 0 : rootFrequencies_->getNumberOfParameters(); // Root frequencies are the first parameters! We should ignore them here.
  for (size_t i = 0; i < modelParameterNames_.size(); i++)
    {
      // Check associations:
      vector<size_t>* modelIndexes = &paramToModels_[i];
      for (size_t j = 0; j < modelIndexes->size(); j++)
        {
          if (!modelParameters_[(*modelIndexes)[j]].hasParameter(modelParameterNames_[i]))
            {
              cerr << "DEBUG: Error, no parameter with name " << modelParameterNames_[i] << " for model " << (*modelIndexes)[j] << endl;
            }
          if (offset + i > getNumberOfParameters()) cerr << "DEBUG: Error, missing parameter " << (offset + i) << "/" << getNumberOfParameters() << endl;
          modelParameters_[(*modelIndexes)[j]].getParameter(modelParameterNames_[i]).setValue(getParameter_(offset + i).getValue());
        }
    }

  // Then we update all models in the set:
  for (size_t i = 0; i < modelParameters_.size(); i++)
    {
      modelSet_[i]->matchParametersValues(modelParameters_[i]);
    }
}

bool SubstitutionModelSet::checkOrphanModels(bool throwEx) const
  throw (Exception)
{
  vector<size_t> index = MapTools::getValues(nodeToModel_);
  for (size_t i = 0; i < modelSet_.size(); i++)
    {
      if (!VectorTools::contains(index, i))
        {
          if (throwEx) throw Exception("SubstitutionModelSet::checkOrphanModels(). Model '" + TextTools::toString(i + 1) + "' is associated to no node.");
          return false;
        }
    }
  return true;
}

bool SubstitutionModelSet::checkOrphanParameters(bool throwEx) const
  throw (Exception)
{
  for (size_t i = 0; i < paramToModels_.size(); i++)
    {
      if (paramToModels_[i].size() == 0)
        {
          if (throwEx) throw Exception("SubstitutionModelSet::checkOrphanParameters(). Parameter '" + getParameter_(i).getName() + "' is associated to no model.");
          return false;
        }
    }
  return true;
}

bool SubstitutionModelSet::checkOrphanNodes(const Tree& tree, bool throwEx) const
  throw (Exception)
{
  vector<int> ids = tree.getNodesId();
  int rootId = tree.getRootId();
  for (size_t i = 0; i < ids.size(); i++)
    {
      if (ids[i] != rootId && nodeToModel_.find(ids[i]) == nodeToModel_.end())
        {
          if (throwEx) throw Exception("SubstitutionModelSet::checkOrphanNodes(). Node '" + TextTools::toString(ids[i]) + "' in tree has no model associated.");
          return false;
        }
    }
  return true;
}

bool SubstitutionModelSet::checkUnknownNodes(const Tree& tree, bool throwEx) const
  throw (Exception)
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
              if (throwEx) throw Exception("SubstitutionModelSet::checkUnknownNodes(). Node '" + TextTools::toString(id) + "' is not found in tree or is the root node.");
              return false;
            }
        }
    }
  return true;
}

bool SubstitutionModelSet::hasMixedSubstitutionModel() const
{
  for (size_t i = 0; i < getNumberOfModels(); i++)
    {
      if (dynamic_cast<const MixedSubstitutionModel*>(getModel(i)) != NULL)
        return true;
    }
  return false;
}

