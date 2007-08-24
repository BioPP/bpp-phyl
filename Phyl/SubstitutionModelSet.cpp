//
// File: SubstitutionModelSet.cpp
// Created by: Bastien Bousseau
//             Julien Duteil
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
  _nodeToModel          (set._nodeToModel),
  _paramToModels        (set._paramToModels),
  _modelParameterNames  (set._modelParameterNames),
  _modelParameters      (set._modelParameters),
  _alphabet             (set._alphabet),
  _freqConstraint       (set._freqConstraint),
  _rootFrequencies      (set._rootFrequencies)
{
  //Duplicate all model objects:
  _modelSet.resize(set._modelSet.size());
  for(unsigned int i = 0; i < set._modelSet.size(); i++)
  {
    _modelSet[i] = set._modelSet[i]->clone();
  }
}

SubstitutionModelSet & SubstitutionModelSet::operator=(const SubstitutionModelSet & set)
{
  AbstractParametrizable::operator=(set);
  _nodeToModel         = set._nodeToModel;
  _paramToModels       = set._paramToModels;
  _modelParameterNames = set._modelParameterNames;
  _modelParameters     = set._modelParameters;
  _alphabet            = set._alphabet;
  _freqConstraint      = set._freqConstraint;
  _rootFrequencies     = set._rootFrequencies;
  
  //Duplicate all model objects:
  _modelSet.resize(set._modelSet.size());
  for(unsigned int i = 0; i < set._modelSet.size(); i++)
  {
    _modelSet[i] = set._modelSet[i]->clone();
  }
  return *this;
}

void SubstitutionModelSet::addModel(SubstitutionModel * model, const vector<int> & nodesId, const vector<string> & newParams) throw (Exception)
{
  if(model->getAlphabet()->getAlphabetType() != _alphabet->getAlphabetType())
    throw Exception("SubstitutionModelSet::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same alphabet");
  _modelSet.push_back(model);
  unsigned int thisModelIndex = _modelSet.size() - 1;
  //Associate this model to specified nodes:
  for(unsigned int i = 0; i < nodesId.size(); i++)
  {
    _nodeToModel[nodesId[i]] = thisModelIndex;
  }
  //Associate parameters:
  string pname;
  _modelParameters.push_back(ParameterList());
  for(unsigned int i  = 0; i < newParams.size(); i++)
  {
    pname = newParams[i];
    _modelParameterNames.push_back(pname);
    Parameter p(model->getParameter(pname));
    p.setName(pname + TextTools::toString(thisModelIndex)); //Change name to unique name in case of sared parameters.
    _parameters.addParameter(p);
    _modelParameters[thisModelIndex].addParameter(p);
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
   
const SubstitutionModel * SubstitutionModelSet::getModelForNode(int nodeId) const
{
  if(_nodeToModel.find(nodeId) == _nodeToModel.end())
    return NULL;
  else
    return _modelSet[_nodeToModel[nodeId]];
}

SubstitutionModel * SubstitutionModelSet::getModelForNode(int nodeId)
{
  if(_nodeToModel.find(nodeId) == _nodeToModel.end())
    return NULL;
  else
    return _modelSet[_nodeToModel[nodeId]];
}

void SubstitutionModelSet::listModelNames(ostream & out) const
{
  for(unsigned int i = 0; i < _modelSet.size(); i++)
    out << "Model "<< i + 1 << " : " << _modelSet[i]->getName() << endl;
}

void SubstitutionModelSet::addParameter(const Parameter & parameter, const vector<int> & nodesId) throw (Exception)
{
  _modelParameterNames.push_back(parameter.getName());
  Parameter p(parameter);
  p.setName(p.getName() + TextTools::toString(_modelSet.size()));
  _parameters.addParameter(p);
  //Build model indexes:
  vector<unsigned int> modelIndexes(nodesId.size());
  for(unsigned int i = 0; i < nodesId.size(); i++)
  {
    modelIndexes[i] = _nodeToModel[nodesId[i]];
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
  string suffix = TextTools::toString(_modelSet.size());
  for(unsigned int i = 0; i < pl.size(); i++)
  {
    pl[i]->setName(pl[i]->getName() + suffix);
  }
  _parameters.addParameters(pl);
  //Build model indexes:
  vector<unsigned int> modelIndexes(nodesId.size());
  for(unsigned int i = 0; i < nodesId.size(); i++)
  {
    modelIndexes[i] = _nodeToModel[nodesId[i]];
  }
  for(unsigned int i = 0; i < pl.size(); i++)
  {
    _paramToModels.push_back(modelIndexes);
  }
}

void SubstitutionModelSet::updateRootFrequencies()
{
  for(unsigned int i = 0; i < _alphabet->getSize(); i++)
  {
    _rootFrequencies[i] = _parameters.getParameter("RootFreq" + TextTools::toString(i))->getValue();
  }
}

void SubstitutionModelSet::fireParameterChanged(const ParameterList & parameters)
{
  //For now, we actualize all parameters... we'll optimize later!
  
  //Update root frequencies:
  updateRootFrequencies();

  //First we actualize the _modelParameters array:
  for(unsigned int i = _alphabet->getSize(); i < _parameters.size(); i++)
  {
    //Check associations:
    vector<unsigned int> * modelIndexes = & _paramToModels[i];
    for(unsigned int j = 0; j < modelIndexes->size(); j++)
    {
      _modelParameters[(*modelIndexes)[j]].getParameter(_modelParameterNames[i])->setValue(_parameters[i]->getValue());
    }
  }

  //Then we update all models in the set:
  for(unsigned int i = 0; i < _modelParameters.size(); i++)
  {
    _modelSet[i]->matchParametersValues(_modelParameters[i]);
  }
}

