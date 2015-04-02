//
// File: NonHomogeneousSubstitutionProcess.cpp
// Created by: Bastien Boussau
//             Julien Dutheil
//             Laurent Guéguen
// Created on: vendredi 21 juin 2013, à 12h 16
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

#include "NonHomogeneousSubstitutionProcess.h"
#include "../Model/MixedSubstitutionModel.h"

#include <Bpp/Utils/MapTools.h>

using namespace bpp;
using namespace std;

NonHomogeneousSubstitutionProcess::NonHomogeneousSubstitutionProcess(const NonHomogeneousSubstitutionProcess& set) :
  AbstractParameterAliasable(set),
  AbstractSubstitutionProcess(set),
  modelSet_(set.modelSet_.size()),
  rootFrequencies_(set.stationarity_ ? 0 : dynamic_cast<FrequenciesSet*>(set.rootFrequencies_->clone())),
  rDist_                (dynamic_cast<DiscreteDistribution*>(set.rDist_->clone())),
  nodeToModel_          (set.nodeToModel_),
  modelToNodes_         (set.modelToNodes_),
  modelParameters_      (set.modelParameters_),
  stationarity_         (set.stationarity_),
  computingTree_()
{
  computingTree_.reset(new ComputingTree(*pTree_.get(), *rDist_.get()));

  // Duplicate all model objects:
  for (size_t i = 0; i < set.modelSet_.size(); i++)
  {
    modelSet_[i]=dynamic_cast<SubstitutionModel*>(set.modelSet_[i]->clone());
    computingTree_->addModel(modelSet_[i],set.modelToNodes_[i]);
  }

  computingTree_->checkModelOnEachNode();
}

NonHomogeneousSubstitutionProcess& NonHomogeneousSubstitutionProcess::operator=(const NonHomogeneousSubstitutionProcess& set)
{
  clear();
  
  AbstractParameterAliasable::operator=(set);
  AbstractSubstitutionProcess::operator=(set);
  nodeToModel_         = set.nodeToModel_;
  modelToNodes_        = set.modelToNodes_;
  modelParameters_     = set.modelParameters_;
  stationarity_        = set.stationarity_;

  if (set.stationarity_)
    rootFrequencies_.reset(0);
  else
    rootFrequencies_.reset(dynamic_cast<FrequenciesSet*>(set.rootFrequencies_->clone()));

  rDist_.reset(dynamic_cast<DiscreteDistribution*>(set.rDist_->clone()));
  
  computingTree_.reset(new ComputingTree(*pTree_.get(), *rDist_.get()));

  // Duplicate all model objects:

  modelSet_.resize(set.modelSet_.size());

  for (size_t i = 0; i < set.modelSet_.size(); i++)
  {
    modelSet_[i]=dynamic_cast<SubstitutionModel*>(set.modelSet_[i]->clone());
    computingTree_->addModel(modelSet_[i],set.modelToNodes_[i]);
  }

  computingTree_->checkModelOnEachNode();
  
  return *this;
}

void NonHomogeneousSubstitutionProcess::clear()
{
  resetParameters_();

  for (size_t i = 0; i < modelSet_.size(); i++)
    if (modelSet_[i])
      delete modelSet_[i];

  modelSet_.clear();
  rootFrequencies_.reset();
  rDist_.reset();
  nodeToModel_.clear();
  modelParameters_.clear();
  computingTree_.release();
  
  stationarity_=true;
}

void NonHomogeneousSubstitutionProcess::setRootFrequencies(FrequenciesSet* rootFreqs)
{
  if (rootFreqs){
    stationarity_=false;
    rootFrequencies_.reset(rootFreqs);
    addParameters_(rootFrequencies_->getIndependentParameters());
  }
}


void NonHomogeneousSubstitutionProcess::setModelToNode(size_t modelIndex, int nodeNumber)
{
  if (modelIndex >= nodeToModel_.size()) throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::setModelToNode.", modelIndex, 0, nodeToModel_.size() - 1);
  nodeToModel_[nodeNumber] = modelIndex;
  
  vector<int> vNod;
  vNod.push_back(nodeNumber);

  computingTree_->addModel(modelSet_[modelIndex], vNod);
}

 
void NonHomogeneousSubstitutionProcess::addModel(SubstitutionModel* model, const std::vector<int>& nodesId)
{
  if (modelSet_.size() > 0 && model->getAlphabet()->getAlphabetType() != modelSet_[0]->getAlphabet()->getAlphabetType())
    throw Exception("NonHomogeneousSubstitutionProcess::addModel. A Substitution Model cannot be added to a Substituion Process if it does not have the same alphabet.");
  if (modelSet_.size() > 0 && model->getNumberOfStates() != modelSet_[0]->getNumberOfStates())
    throw Exception("NonHomogeneousSubstitutionProcess::addModel. A Substitution Model cannot be added to a Substitution Process if it does not have the same number of states.");

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
  ParameterList pl=model->getIndependentParameters();
  modelParameters_.push_back(pl);
   
  for (size_t i  = 0; i < pl.size(); i++)
    {
      Parameter* p = pl[i].clone();
      p->setName(p->getName() + "_" + TextTools::toString(modelParameters_.size()));
      addParameter_(p);
    }

  computingTree_->addModel(model, nodesId);
  computingTree_->checkModelOnEachNode();
}


void NonHomogeneousSubstitutionProcess::setModel(SubstitutionModel* model, size_t modelIndex)
{
  if (modelSet_.size() > 0 && model->getAlphabet()->getAlphabetType() != modelSet_[0]->getAlphabet()->getAlphabetType())
    throw Exception("NonHomogeneousSubstitutionProcess::addModel. A Substitution Model cannot be added to a Substituion Process if it does not have the same alphabet.");
  if (modelSet_.size() > 0 && model->getNumberOfStates() != modelSet_[0]->getNumberOfStates())
    throw Exception("NonHomogeneousSubstitutionProcess::addModel. A Substitution Model cannot be added to a Substitution Process if it does not have the same number of states.");

  if (modelIndex >= modelSet_.size())
    throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::setModel.", modelIndex, 0, modelSet_.size());
  
  if (modelSet_[modelIndex])
    delete modelSet_[modelIndex];
  modelSet_[modelIndex]=model;

  // Change associate parameters
  ParameterList& pl1=modelParameters_[modelIndex];
  for (size_t i=0; i<pl1.size(); i++){
    string pn=pl1[i].getName()+ "_" + TextTools::toString(modelIndex+1);
    deleteParameter_(pn);
  }
  string pname;
  ParameterList pl=model->getIndependentParameters();
  modelParameters_[modelIndex]=pl;
  
  for (size_t i  = 0; i < pl.size(); i++)
    {
      Parameter* p = pl[i].clone();
      p->setName(p->getName() + "_" + TextTools::toString(modelIndex+1));
      addParameter_(p);
    }

  computingTree_->addModel(model, modelToNodes_[modelIndex]);
  computingTree_->checkModelOnEachNode();
}

void NonHomogeneousSubstitutionProcess::listModelNames(std::ostream& out) const
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

void NonHomogeneousSubstitutionProcess::fireParameterChanged(const ParameterList& parameters)
{
  // Update root frequencies:
  updateRootFrequencies();

  //Update rate distribution:
  rDist_->matchParametersValues(parameters);


  // Then we update all models in the set:
  for (size_t i = 0; i < modelParameters_.size(); i++)
    {
      for (size_t np = 0 ; np< modelParameters_[i].size() ; np++)
        {
          modelParameters_[i][np].setValue(getParameterValue(modelParameters_[i][np].getName()+"_"+TextTools::toString(i+1)));
        }
      if (modelSet_[i]->matchParametersValues(modelParameters_[i]))
        computingTree_->update(modelToNodes_[i]);
    }

  AbstractSubstitutionProcess::fireParameterChanged(parameters);
}


ParameterList NonHomogeneousSubstitutionProcess::getSubstitutionModelParameters(bool independent) const
{
  ParameterList pl;
  
  // Then we update all models in the set:
  for (size_t i = 0; i < modelParameters_.size(); i++)
    {
      for (size_t np = 0 ; np< modelParameters_[i].size() ; np++)
        {
          Parameter p(modelParameters_[i][np]);
          p.setName(p.getName()+"_"+TextTools::toString(i+1));
          pl.addParameter(p);
        }
    }

  return pl;
}

bool NonHomogeneousSubstitutionProcess::checkOrphanNodes(bool throwEx) const
{
  vector<int> ids = getTree().getNodesId();
  int rootId = getTree().getRootId();
  for (size_t i = 0; i < ids.size(); i++)
    {
      if (ids[i] != rootId && nodeToModel_.find(ids[i]) == nodeToModel_.end())
        {
          if (throwEx) throw Exception("NonHomogeneousSubstitutionProcess::checkOrphanNodes(). Node '" + TextTools::toString(ids[i]) + "' in tree has no model associated.");
          return false;
        }
    }
  return true;
}

bool NonHomogeneousSubstitutionProcess::checkUnknownNodes(bool throwEx) const
{
  vector<int> ids = getTree().getNodesId();
  int id;
  int rootId = getTree().getRootId();
  std::map<size_t, std::vector<int> >::const_iterator it;
  
  for (it=modelToNodes_.begin(); it!=modelToNodes_.end(); it++)
    {
      for (size_t j = 0; j < it->second.size(); j++)
        {
          id = it->second[j];
          if (id == rootId || !VectorTools::contains(ids, id))
            {
              if (throwEx)
                throw Exception("NonHomogeneousSubstitutionProcess::checkUnknownNodes(). Node '" + TextTools::toString(id) + "' is not found in tree or is the root node.");
              return false;
            }
        }
    }
  return true;
}

bool NonHomogeneousSubstitutionProcess::hasMixedSubstitutionModel() const
{
  for (size_t i = 0; i < getNumberOfModels(); i++)
    {
      if (dynamic_cast<const MixedSubstitutionModel*>(getModel(i)) != NULL)
        return true;
    }
  return false;
}

/*
 * Inheriting from SubstitutionProcess
 */
  
bool NonHomogeneousSubstitutionProcess::isCompatibleWith(const SiteContainer& data) const
{
  if (modelSet_.size() > 0) 
    return data.getAlphabet()->getAlphabetType() == modelSet_[0]->getAlphabet()->getAlphabetType();
  else
    return true;
}


bool NonHomogeneousSubstitutionProcess::hasDerivableParameter(const std::string& name) const
{
  // Up to now (!), only branch length parameters are derivable
  return (name.substr(0,5)=="BrLen");
}

NonHomogeneousSubstitutionProcess* NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(
                                                                                                                SubstitutionModel* model,
                                                                                                                DiscreteDistribution* rdist,
                                                                                                                FrequenciesSet* rootFreqs,
                                                                                                                ParametrizableTree* tree
                                                                                                                )
{
  // Check alphabet:
  if  (rootFreqs && model->getAlphabet()->getAlphabetType() != rootFreqs->getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("NonHomogeneousSubstitutionProcess::createHomogeneousModelSet()", model->getAlphabet(), rootFreqs->getAlphabet());
  if (dynamic_cast<MixedSubstitutionModel*>(model) != NULL)
    throw Exception("createHomogeneousSubstitutionProcess not yet programmed for mixture models.");

  NonHomogeneousSubstitutionProcess*  modelSet = new NonHomogeneousSubstitutionProcess(rdist, tree, rootFreqs);

  // We assign this model to all nodes in the tree (excepted root node), and link all parameters with it.
  vector<int> ids = tree->getTree().getNodesId();
  int rootId = tree->getTree().getRootId();
  unsigned int pos = 0;
  for (unsigned int i = 0; i < ids.size(); i++)
    {
      if (ids[i] == rootId)
        {
          pos = i;
          break;
        }
    }
  ids.erase(ids.begin() + pos);
  modelSet->addModel(model, ids);
  return modelSet;
}

NonHomogeneousSubstitutionProcess* NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(
                                                                                                                   SubstitutionModel* model,
                                                                                                                   DiscreteDistribution* rdist,
                                                                                                                   FrequenciesSet* rootFreqs,
                                                                                                                   ParametrizableTree* tree,
                                                                                                                   const vector<string>& globalParameterNames
                                                                                                                   )
{
  // Check alphabet:
  if (rootFreqs && model->getAlphabet()->getAlphabetType() != rootFreqs->getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("NonHomogeneousSubstitutionProcess::createNonHomogeneousModelSet()", model->getAlphabet(), rootFreqs->getAlphabet());
  if (dynamic_cast<MixedSubstitutionModel*>(model) != NULL)
    throw Exception("createNonHomogeneousSubstitutionProcess not yet programmed for mixture models.");

  ParameterList globalParameters, branchParameters;
  globalParameters = model->getIndependentParameters();

  vector<string> globalParameterNames2; // vector of the complete names of global parameters

  // First get correct parameter names
  size_t i, j;

  for (i = 0; i < globalParameterNames.size(); i++)
    {
      if (globalParameterNames[i].find("*") != string::npos)
        {
          for (j = 0; j < globalParameters.size(); j++)
            {
              StringTokenizer stj(globalParameterNames[i], "*", true, false);
              size_t pos1, pos2;
              string parn = globalParameters[j].getName();
              bool flag(true);
              string g = stj.nextToken();
              pos1 = parn.find(g);
              if (pos1 != 0)
                flag = false;
              pos1 += g.length();
              while (flag && stj.hasMoreToken())
                {
                  g = stj.nextToken();
                  pos2 = parn.find(g, pos1);
                  if (pos2 == string::npos)
                    {
                      flag = false;
                      break;
                    }
                  pos1 = pos2 + g.length();
                }
              if (flag &&
                  ((g.length() == 0) || (pos1 == parn.length()) || (parn.rfind(g) == parn.length() - g.length())))
                globalParameterNames2.push_back(parn);
            }
        }
      else if (!globalParameters.hasParameter(globalParameterNames[i]))
        throw Exception("NonHomogeneousSubstitutionProcess::createNonHomogeneousModelSet. Parameter '" + globalParameterNames[i] + "' is not valid.");
      else
        globalParameterNames2.push_back(globalParameterNames[i]);
    }

  // remove non global parameters
  for (i = globalParameters.size(); i > 0; i--)
    {
      if (find(globalParameterNames2.begin(), globalParameterNames2.end(), globalParameters[i - 1].getName()) == globalParameterNames2.end())
        {
          // not a global parameter:
          branchParameters.addParameter(globalParameters[i - 1]);
          globalParameters.deleteParameter(i - 1);
        }
    }

  // bool mixed = (dynamic_cast<MixedSubstitutionModel*>(model) != NULL);
  NonHomogeneousSubstitutionProcess*  modelSet;
  // if (mixed)
  // {
  //   modelSet = new MixedNonHomogeneousSubstitutionProcess(model->getAlphabet());
  //   // Remove the "relproba" parameters from the branch parameters and put them in the global parameters, for the hypernodes
  //   for (i = branchParameters.size(); i > 0; i--)
  //   {
  //     if (branchParameters[i - 1].getName().find("relproba") != string::npos)
  //     {
  //       globalParameters.addParameter(branchParameters[i - 1]);
  //       branchParameters.deleteParameter(i - 1);
  //     }
  //   }
  // }
  // else
  modelSet = new NonHomogeneousSubstitutionProcess(rdist, tree, rootFreqs);

  // We assign a copy of this model to all nodes in the tree (excepted root node), and link all parameters with it.
  vector<int> ids = tree->getTree().getNodesId();
  int rootId = tree->getTree().getRootId();
  size_t pos = 0;
  for (i = 0; i < ids.size(); i++)
    {
      if (ids[i] == rootId)
        {
          pos = i;
          break;
        }
    }

  ids.erase(ids.begin() + pos);
  for (i = 0; i < ids.size(); i++)
    {
      modelSet->addModel(dynamic_cast<SubstitutionModel*>(model->clone()), vector<int>(1, ids[i]));
    }

  // Now alias all global parameters on all nodes:
  for (i=0; i < globalParameters.size(); i++)
    {
      string pname=globalParameters[i].getName();

      for (size_t nn = 1; nn < ids.size(); nn++)
        modelSet->aliasParameters(pname+"_1",pname+"_"+TextTools::toString(nn+1));
    }
  
  // Defines the hypernodes if mixed
  // if (mixed)
  // {
  //   MixedNonHomogeneousSubstitutionProcess* pMSMS = dynamic_cast<MixedNonHomogeneousSubstitutionProcess*>(modelSet);
  //   MixedSubstitutionModel* pMSM = dynamic_cast<MixedSubstitutionModel*>(model);

  //   size_t nbm = pMSM->getNumberOfModels();
  //   for (i = 0; i < nbm; i++)
  //   {
  //     pMSMS->addEmptyHyperNode();
  //     for (j = 0; j < ids.size(); j++)
  //     {
  //       pMSMS->addToHyperNode(j, vector<int>(1, static_cast<int>(i)));
  //     }
  //   }
  //   pMSMS->computeHyperNodesProbabilities();
  // }

  // delete model; // delete template model.
  return modelSet;
}

