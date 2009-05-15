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

// From Utils:
#include <Utils/MapTools.h>

using namespace bpp;

/*
void SubstitutionModelSet::makeRandomNodeToModelLinks()
{
  int rand;
  for (int i=0;i<_numberOfNodes; i++) {
    rand = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(_numberOfModels);
    setModelToNode(rand, i);
  }
}*/

SubstitutionModelSet::SubstitutionModelSet(const SubstitutionModelSet & set):
  AbstractParametrizable(set),
  _alphabet             (set._alphabet),
  _nodeToModel          (set._nodeToModel),
  _modelToNodes         (set._modelToNodes),
  _paramToModels        (set._paramToModels),
  _paramNamesCount      (set._paramNamesCount),
  _modelParameterNames  (set._modelParameterNames),
  _modelParameters      (set._modelParameters)
{
  //Duplicate all model objects:
  _modelSet.resize(set._modelSet.size());
  for(unsigned int i = 0; i < set._modelSet.size(); i++)
  {
    _modelSet[i] = dynamic_cast<SubstitutionModel *>(set._modelSet[i]->clone());
  }
  _rootFrequencies = dynamic_cast<FrequenciesSet *>(set._rootFrequencies->clone());
}

SubstitutionModelSet & SubstitutionModelSet::operator=(const SubstitutionModelSet & set)
{
  AbstractParametrizable::operator=(set);
  _alphabet            = set._alphabet;
  _nodeToModel         = set._nodeToModel;
  _modelToNodes        = set._modelToNodes;
  _paramToModels       = set._paramToModels;
  _paramNamesCount     = set._paramNamesCount;
  _modelParameterNames = set._modelParameterNames;
  _modelParameters     = set._modelParameters;
  _rootFrequencies     = dynamic_cast<FrequenciesSet *>(set._rootFrequencies->clone());
  
  //Duplicate all model objects:
  _modelSet.resize(set._modelSet.size());
  for(unsigned int i = 0; i < set._modelSet.size(); i++)
  {
    _modelSet[i] = dynamic_cast<SubstitutionModel *>(set._modelSet[i]->clone());
  }
  return *this;
}

vector<int> SubstitutionModelSet::getNodesWithParameter(const string & name) const
  throw (ParameterNotFoundException)
{
  vector<int> ids;
  unsigned int offset = _rootFrequencies->getNumberOfParameters();
  for(unsigned int i = 0; i < _paramToModels.size(); i++)
  {
    if(getParameter_(offset + i).getName() == name)
    { 
      for(unsigned int j = 0; j < _paramToModels[i].size(); j++)
      {
        vector<int> v = _modelToNodes[_paramToModels[i][j]];
        VectorTools::append(ids, v);
      }
      return ids;
    }
  }
  throw ParameterNotFoundException("SubstitutionModelSet::getNodesWithParameter.", name);
}

vector<unsigned int> SubstitutionModelSet::getModelsWithParameter(const string & name) const
  throw (ParameterNotFoundException)
{
  vector<unsigned int> indices;
  unsigned int offset = _rootFrequencies->getNumberOfParameters();
  for(unsigned int i = 0; i < _paramToModels.size(); i++)
  {
    if(getParameter_(offset + i).getName() == name)
      return _paramToModels[i];
  }
  throw ParameterNotFoundException("SubstitutionModelSet::getModelsWithParameter.", name);
}

void SubstitutionModelSet::addModel(SubstitutionModel *model, const vector<int> & nodesId, const vector<string> & newParams) throw (Exception)
{
  if(model->getAlphabet()->getAlphabetType() != _alphabet->getAlphabetType())
    throw Exception("SubstitutionModelSet::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same alphabet.");
  if(_modelSet.size() > 0 && model->getNumberOfStates() != _modelSet[0]->getNumberOfStates())
    throw Exception("SubstitutionModelSet::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same number of states.");
  _modelSet.push_back(model);
  unsigned int thisModelIndex = _modelSet.size() - 1;

  //Associate this model to specified nodes:
  for(unsigned int i = 0; i < nodesId.size(); i++)
  {
    _nodeToModel[nodesId[i]] = thisModelIndex;
    _modelToNodes[thisModelIndex].push_back(nodesId[i]);
  }

  //Associate parameters:
  string pname;
  _modelParameters.push_back(ParameterList());
  for(unsigned int i  = 0; i < newParams.size(); i++)
  {
    pname = newParams[i];
    _modelParameterNames.push_back(pname);
    vector<unsigned int> modelsIndex(1, thisModelIndex);
    _paramToModels.push_back(modelsIndex);
    Parameter p(model->getParameters().getParameter(pname)); //We work with namespaces here, so model->getParameter(pname) does not work.
    _modelParameters[thisModelIndex].addParameter(p);
    p.setName(pname + "_" + TextTools::toString(++_paramNamesCount[pname])); //Change name to unique name in case of shared parameters.
    addParameter_(p);
  }
}

void SubstitutionModelSet::setModel(SubstitutionModel * model, unsigned int modelIndex) throw (Exception, IndexOutOfBoundsException)
{
  if(model->getAlphabet()->getAlphabetType() != _alphabet->getAlphabetType())
    throw Exception("SubstitutionModelSet::setModel. A Substitution Model cannot be added to a Model Set if it does not have the same alphabet");
  if(modelIndex >= _modelSet.size())
    throw IndexOutOfBoundsException("SubstitutionModelSet::setModel.", modelIndex, 0, _modelSet.size());
  delete _modelSet[modelIndex];
  _modelSet[modelIndex] = model;
}
   
void SubstitutionModelSet::removeModel(unsigned int modelIndex) throw (Exception)
{
  _modelSet.erase(_modelSet.begin() + modelIndex);
  //Erase all parameter references to this model and translate other indices...
  for(unsigned int i = 0; i < _paramToModels.size(); i++)
  {
    for(unsigned int j = _paramToModels[i].size(); j > 0; j--)
    {
      if(_paramToModels[i][j-1] == modelIndex)
        _paramToModels[i].erase(_paramToModels[i].begin() + j - 1);
      else if(_paramToModels[i][j-1] > modelIndex)
        _paramToModels[i][j-1]--; //Correct indice due to removal!
    }
  }
  checkOrphanParameters(true);
}

void SubstitutionModelSet::listModelNames(ostream & out) const
{
  for(unsigned int i = 0; i < _modelSet.size(); i++)
  {
    out << "Model "<< i + 1 << ": " << _modelSet[i]->getName() << "\t attached to nodes ";
    for(unsigned int j = 0; j < _modelToNodes[i].size(); j++)
      out << _modelToNodes[i][j];
    out << endl;
  }
}

void SubstitutionModelSet::setParameterToModel(unsigned int parameterIndex, unsigned int modelIndex) throw (IndexOutOfBoundsException)
{
  unsigned int offset = _rootFrequencies->getNumberOfParameters();
  if(parameterIndex < offset) throw IndexOutOfBoundsException("SubstitutionModelSet::setParameterToModel. Can't assign a root frequency parameter to a branch model.", parameterIndex, offset - 1, _paramToModels.size() + offset - 1);
  if(parameterIndex >= _paramToModels.size() + offset) throw IndexOutOfBoundsException("SubstitutionModelSet::setParameterToModel.", parameterIndex, offset - 1, _paramToModels.size() + offset - 1);
  if(modelIndex >= _modelSet.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::setParameterToModel.", modelIndex, 0, _modelSet.size() - 1);
  if(VectorTools::contains(_paramToModels[parameterIndex - offset], modelIndex)) throw Exception("SubstitutionModelSet::setParameterToModel: parameter " + getParameter_(parameterIndex).getName() + " already set to model " + TextTools::toString(modelIndex) + ".");
  _paramToModels[parameterIndex - offset].push_back(modelIndex);
  _modelParameters[modelIndex].addParameter(
      _modelSet[modelIndex]->getParameters().getParameter(_modelParameterNames[parameterIndex - offset]));
}

void SubstitutionModelSet::unsetParameterToModel(unsigned int parameterIndex, unsigned int modelIndex) throw (IndexOutOfBoundsException, Exception)
{
  unsigned int offset = _rootFrequencies->getNumberOfParameters();
  if(parameterIndex < offset) throw IndexOutOfBoundsException("SubstitutionModelSet::unsetParameterToModel. Can't unset a root frequency parameter.", parameterIndex, offset - 1, _paramToModels.size() + offset - 1);
  if(parameterIndex >= _paramToModels.size() + offset) throw IndexOutOfBoundsException("SubstitutionModelSet::unsetParameterToModel.", parameterIndex, offset - 1, _paramToModels.size() + offset - 1);
  if(modelIndex >= _modelSet.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::setParameterToModel.", modelIndex, 0, _modelSet.size() - 1);
  if(!VectorTools::contains(_paramToModels[parameterIndex - offset], modelIndex)) throw Exception("SubstitutionModelSet::unsetParameterToModel: parameter " + getParameter_(parameterIndex).getName() + " is not currently set to model " + TextTools::toString(modelIndex) + ".");
  remove(_paramToModels[parameterIndex - offset].begin(), _paramToModels[parameterIndex - offset].end(), modelIndex);
  _modelParameters[modelIndex].deleteParameter(_modelParameterNames[parameterIndex - offset]);
  checkOrphanModels(true);
  checkOrphanParameters(true);
}

void SubstitutionModelSet::addParameter(const Parameter & parameter, const vector<int> & nodesId) throw (Exception)
{
  _modelParameterNames.push_back(parameter.getName());
  Parameter p(parameter);
  p.setName(p.getName() + "_" + TextTools::toString(++_paramNamesCount[p.getName()]));
  addParameter_(p);
  //Build model indexes:
  vector<unsigned int> modelIndexes(nodesId.size());
  for(unsigned int i = 0; i < nodesId.size(); i++)
  {
    unsigned int pos = _nodeToModel[nodesId[i]];
    _modelParameters[pos].addParameter(parameter);
    modelIndexes[i] = pos;
  }
  _paramToModels.push_back(modelIndexes);
}

void SubstitutionModelSet::addParameters(const ParameterList & parameters, const vector<int> & nodesId) throw (Exception)
{
  for(unsigned int i = 0; i < parameters.size(); i++)
  {
    _modelParameterNames.push_back(parameters[i]->getName());
  }
  ParameterList pl(parameters);
  for(unsigned int i = 0; i < pl.size(); i++)
  {
    pl[i]->setName(pl[i]->getName() + "_" + TextTools::toString(++_paramNamesCount[pl[i]->getName()]));
  }
  addParameters_(pl);
  //Build model indexes:
  vector<unsigned int> modelIndexes(nodesId.size());
  for(unsigned int i = 0; i < nodesId.size(); i++)
  {
    unsigned int pos = _nodeToModel[nodesId[i]];
    _modelParameters[pos].addParameters(parameters);
    modelIndexes[i] = pos;
  }
  for(unsigned int i = 0; i < pl.size(); i++)
  {
    _paramToModels.push_back(modelIndexes);
  }
}

void SubstitutionModelSet::removeParameter(const string & name) throw (ParameterNotFoundException)
{
  unsigned int offset = _rootFrequencies->getNumberOfParameters();
  for (unsigned int i = 0; i < _modelParameterNames.size(); i++)
  {
    if(getParameter_(offset + i).getName() == name)
    {
      vector<int> nodesId = getNodesWithParameter(name);
      for (unsigned int j = 0; j < nodesId.size(); j++)
      {
        unsigned int pos = _nodeToModel[nodesId[j]];
        string tmp = _modelParameterNames[i];
        if (_modelParameters[pos].hasParameter(tmp))
          _modelParameters[pos].deleteParameter(tmp);
      }
      _paramToModels.erase(_paramToModels.begin() + i);
      _modelParameterNames.erase(_modelParameterNames.begin() + i);
      deleteParameter_(offset + i);
      return;
    }
  }
  throw ParameterNotFoundException("SubstitutionModelSet::removeParameter.", name);
}

void SubstitutionModelSet::fireParameterChanged(const ParameterList & parameters)
{
  //For now, we actualize all parameters... we'll optimize later!
  //Update root frequencies:
  updateRootFrequencies();

  //First we actualize the _modelParameters array:
  unsigned int offset = _rootFrequencies->getNumberOfParameters(); //Root frequencies are the first parameters! We should ignore them here.
  for(unsigned int i = 0; i < _modelParameterNames.size(); i++)
  {
    //Check associations:
    vector<unsigned int> * modelIndexes = & _paramToModels[i];
    for(unsigned int j = 0; j < modelIndexes->size(); j++)
    {
      if (!_modelParameters[(*modelIndexes)[j]].hasParameter(_modelParameterNames[i]))
      {
        cerr << "DEBUG: Error, no parameter with name " << _modelParameterNames[i] << " for model " << (*modelIndexes)[j] << endl;
      }
      if(offset + i > getNumberOfParameters()) cerr << "DEBUG: Error, missing parameter " << (offset + i) << "/" << getNumberOfParameters() << endl;
      _modelParameters[(*modelIndexes)[j]].getParameter(_modelParameterNames[i]).setValue(getParameter_(offset + i).getValue());
    }
  }

  //Then we update all models in the set:
  for(unsigned int i = 0; i < _modelParameters.size(); i++)
  {
    _modelSet[i]->matchParametersValues(_modelParameters[i]);
  }
}

bool SubstitutionModelSet::checkOrphanModels(bool throwEx) const
throw (Exception)
{
  vector<unsigned int> index = MapTools::getValues(_nodeToModel);
  for(unsigned int i = 0; i < _modelSet.size(); i++)
  {
    if(!VectorTools::contains(index, i))
    {
      if(throwEx) throw Exception("SubstitutionModelSet::checkOrphanModels(). Model '" + TextTools::toString(i+1) + "' is associated to no node.");
      return false;
    }
  }
  return true;
}

bool SubstitutionModelSet::checkOrphanParameters(bool throwEx) const
throw (Exception)
{
  for(unsigned int i = 0; i < _paramToModels.size(); i++)
  {
    if(_paramToModels[i].size() == 0)
    {
      if(throwEx) throw Exception("SubstitutionModelSet::checkOrphanParameters(). Parameter '" + getParameter_(i).getName() + "' is associated to no model.");
      return false;
    }
  }
  return true;
}

bool SubstitutionModelSet::checkOrphanNodes(const Tree & tree, bool throwEx) const
throw (Exception)
{
  vector<int> ids = tree.getNodesId();
  int rootId = tree.getRootId();
  for(unsigned int i = 0; i < ids.size(); i++)
  {
    if(ids[i] != rootId && _nodeToModel.find(ids[i]) == _nodeToModel.end())
    {
      if(throwEx) throw Exception("SubstitutionModelSet::checkOrphanNodes(). Node '" + TextTools::toString(ids[i]) + "' in tree has no model associated.");
      return false;
    }
  }
  return true;
}

bool SubstitutionModelSet::checkUnknownNodes(const Tree & tree, bool throwEx) const
throw (Exception)
{
  vector<int> ids = tree.getNodesId();
  int id;
  int rootId = tree.getRootId();
  for(unsigned int i = 0; i < _modelToNodes.size(); i++)
  {
    for(unsigned int j = 0; j < _modelToNodes[i].size(); j++)
    {
      id = _modelToNodes[i][j];
      if(id == rootId || !VectorTools::contains(ids, id))
      {
        if(throwEx) throw Exception("SubstitutionModelSet::checkUnknownNodes(). Node '" + TextTools::toString(id) + "' is not found in tree or is the root node.");
        return false;
      }
    }
  }
  return true;
}

