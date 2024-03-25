// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Utils/MapTools.h>
#include <algorithm>

#include "../Model/MixedTransitionModel.h"
#include "NonHomogeneousSubstitutionProcess.h"
#include "RateAcrossSitesSubstitutionProcess.h"
#include "SimpleSubstitutionProcess.h"

using namespace bpp;
using namespace std;

NonHomogeneousSubstitutionProcess::NonHomogeneousSubstitutionProcess(const NonHomogeneousSubstitutionProcess& set) :
  AbstractParameterAliasable(set),
  AbstractAutonomousSubstitutionProcess(set),
  modelSet_(set.modelSet_.size()),
  rDist_                (set.rDist_ ? dynamic_cast<DiscreteDistributionInterface*>(set.rDist_->clone()) : 0),
  nodeToModel_          (set.nodeToModel_),
  modelToNodes_         (set.modelToNodes_),
  modelParameters_      (set.modelParameters_)
{
  // Duplicate all model objects:
  for (size_t i = 0; i < set.modelSet_.size(); i++)
  {
    modelSet_[i] = shared_ptr<BranchModelInterface>(set.modelSet_[i]->clone());
  }

  if (modelScenario_)
    for (size_t i = 0; i < modelSet_.size(); i++)
    {
      modelScenario_->changeModel(dynamic_pointer_cast<MixedTransitionModelInterface>(set.modelSet_[i]), dynamic_pointer_cast<MixedTransitionModelInterface>(modelSet_[i]));
    }
}

NonHomogeneousSubstitutionProcess& NonHomogeneousSubstitutionProcess::operator=(const NonHomogeneousSubstitutionProcess& set)
{
  clear();

  AbstractParameterAliasable::operator=(set);
  AbstractAutonomousSubstitutionProcess::operator=(set);
  nodeToModel_         = set.nodeToModel_;
  modelToNodes_        = set.modelToNodes_;
  modelParameters_     = set.modelParameters_;

  rDist_.reset(rDist_ ? dynamic_cast<DiscreteDistributionInterface*>(set.rDist_->clone()) : 0);

  // Duplicate all model objects:

  modelSet_.resize(set.modelSet_.size());

  for (size_t i = 0; i < set.modelSet_.size(); i++)
  {
    modelSet_[i] = shared_ptr<BranchModelInterface>(set.modelSet_[i]->clone());
  }

  if (modelScenario_)
    for (size_t i = 0; i < modelSet_.size(); i++)
    {
      modelScenario_->changeModel(dynamic_pointer_cast<MixedTransitionModelInterface>(set.modelSet_[i]), dynamic_pointer_cast<MixedTransitionModelInterface>(modelSet_[i]));
    }

  return *this;
}

void NonHomogeneousSubstitutionProcess::clear()
{
  resetParameters_();

  modelSet_.clear();
  rDist_.reset();
  nodeToModel_.clear();
  modelParameters_.clear();
}

void NonHomogeneousSubstitutionProcess::setModelToNode(size_t modelIndex, unsigned int nodeNumber)
{
  if (modelIndex >= nodeToModel_.size())
    throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::setModelToNode.", modelIndex, 0, nodeToModel_.size() - 1);
  nodeToModel_[nodeNumber] = modelIndex;

  vector<unsigned int> vNod;
  vNod.push_back(nodeNumber);
}


void NonHomogeneousSubstitutionProcess::addModel(shared_ptr<BranchModelInterface> model, const vector<unsigned int>& nodesId)
{
  if (modelSet_.size() > 0 && model->getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw Exception("NonHomogeneousSubstitutionProcess::addModel. A Substitution Model cannot be added to a Substitution Process if it does not have the same alphabet.");
  if (modelSet_.size() > 0 && model->getNumberOfStates() != getNumberOfStates())
    throw Exception("NonHomogeneousSubstitutionProcess::addModel. A Substitution Model cannot be added to a Substitution Process if it does not have the same number of states.");

  modelSet_.push_back(model);
  size_t thisModelIndex = modelSet_.size() - 1;


  // Associate this model to specified nodes:
  modelToNodes_[thisModelIndex] = nodesId;
  for (size_t i = 0; i < nodesId.size(); i++)
  {
    nodeToModel_[nodesId[i]] = thisModelIndex;
  }

  // Associate parameters:
  string pname;
  ParameterList pl = model->getIndependentParameters();
  modelParameters_.push_back(pl);

  for (size_t i  = 0; i < pl.size(); i++)
  {
    Parameter* p = pl[i].clone();
    p->setName(p->getName() + "_" + TextTools::toString(modelSet_.size()));
    addParameter_(p);
  }
}

void NonHomogeneousSubstitutionProcess::setModel(shared_ptr<BranchModelInterface> model, size_t modelIndex)
{
  if (modelSet_.size() > 0 && model->getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw Exception("NonHomogeneousSubstitutionProcess::setModel. A Substitution Model cannot be added to a Substitution Process if it does not have the same alphabet.");
  if (modelSet_.size() > 0 && model->getNumberOfStates() != getNumberOfStates())
    throw Exception("NonHomogeneousSubstitutionProcess::setModel. A Substitution Model cannot be added to a Substitution Process if it does not have the same number of states.");

  if (modelIndex >= modelSet_.size())
    throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::setModel.", modelIndex, 0, modelSet_.size());

  modelSet_[modelIndex] = model;

  // Change associate parameters
  ParameterList& pl1 = modelParameters_[modelIndex];
  for (size_t i = 0; i < pl1.size(); i++)
  {
    string pn = pl1[i].getName() + "_" + TextTools::toString(modelIndex + 1);
    deleteParameter_(pn);
  }
  string pname;
  ParameterList pl = model->getIndependentParameters();
  modelParameters_[modelIndex] = pl;

  for (size_t i  = 0; i < pl.size(); i++)
  {
    Parameter* p = pl[i].clone();
    p->setName(p->getName() + "_" + TextTools::toString(modelIndex + 1));
    addParameter_(p);
  }
}

void NonHomogeneousSubstitutionProcess::listModelNames(ostream& out) const
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
  // Update rate distribution:
  if (rDist_)
    rDist_->matchParametersValues(parameters);


  // Then we update all models in the set:
  for (size_t i = 0; i < modelParameters_.size(); i++)
  {
    for (size_t np = 0; np < modelParameters_[i].size(); np++)
    {
      modelParameters_[i][np].setValue(getParameterValue(modelParameters_[i][np].getName() + "_" + TextTools::toString(i + 1)));
    }
    modelSet_[i]->matchParametersValues(modelParameters_[i]);
  }

  AbstractAutonomousSubstitutionProcess::fireParameterChanged(parameters);
}


ParameterList NonHomogeneousSubstitutionProcess::getSubstitutionModelParameters(bool independent) const
{
  ParameterList pl;

  // Then we update all models in the set:
  for (size_t i = 0; i < modelParameters_.size(); i++)
  {
    for (size_t np = 0; np < modelParameters_[i].size(); np++)
    {
      if (!independent || hasIndependentParameter(modelParameters_[i][np].getName() + "_" + TextTools::toString(i + 1)))
      {
        Parameter p(modelParameters_[i][np]);
        p.setName(p.getName() + "_" + TextTools::toString(i + 1));
        pl.addParameter(p);
      }
    }
  }

  return pl;
}

bool NonHomogeneousSubstitutionProcess::checkOrphanNodes(bool throwEx) const
{
  if (!getParametrizablePhyloTree())
  {
    if (throwEx)
      throw Exception("NonHomogeneousSubstitutionProcess::checkOrphanNodes(). No tree provided.");
    return false;
  }

  vector<unsigned int> ids = getParametrizablePhyloTree()->getAllNodesIndexes();
  unsigned int rootId = getParametrizablePhyloTree()->getNodeIndex(getParametrizablePhyloTree()->getRoot());
  for (size_t i = 0; i < ids.size(); i++)
  {
    if (ids[i] != rootId && nodeToModel_.find(ids[i]) == nodeToModel_.end())
    {
      if (throwEx)
        throw Exception("NonHomogeneousSubstitutionProcess::checkOrphanNodes(). Node '" + TextTools::toString(ids[i]) + "' in tree has no model associated.");
      return false;
    }
  }
  return true;
}

bool NonHomogeneousSubstitutionProcess::checkUnknownNodes(bool throwEx) const
{
  if (!getParametrizablePhyloTree())
  {
    if (throwEx)
      throw Exception("NonHomogeneousSubstitutionProcess::checkUnknownNodes(). No tree provided.");
    return false;
  }
  vector<unsigned int> ids = getParametrizablePhyloTree()->getAllNodesIndexes();
  unsigned int id;
  unsigned int rootId = getParametrizablePhyloTree()->getNodeIndex(getParametrizablePhyloTree()->getRoot());

  map<size_t, vector<unsigned int>>::const_iterator it;

  for (it = modelToNodes_.begin(); it != modelToNodes_.end(); it++)
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

bool NonHomogeneousSubstitutionProcess::hasMixedTransitionModel() const
{
  for (size_t i = 1; i <= getNumberOfModels(); i++)
  {
    if (dynamic_pointer_cast<const MixedTransitionModelInterface>(getModel(i)) != nullptr)
      return true;
  }
  return false;
}


void NonHomogeneousSubstitutionProcess::setModelScenario(shared_ptr<ModelScenario> modelscenario)
{
  auto vmod = modelscenario->getModels();

  for (auto& mod:vmod)
  {
    if (find(modelSet_.begin(), modelSet_.end(), mod) == modelSet_.end())
      throw Exception("NonHomogeneousSubstitutionProcess::setModelScenario: unknown model " + mod->getName());
  }

  modelScenario_ = modelscenario;
}


unique_ptr<AutonomousSubstitutionProcessInterface> NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(
    shared_ptr<BranchModelInterface> model,
    shared_ptr<DiscreteDistributionInterface> rdist,
    shared_ptr<PhyloTree> tree,
    shared_ptr<FrequencySetInterface> rootFreqs,
    shared_ptr<ModelScenario> scenario)
{
  if (!tree)
    throw Exception("NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess: missing tree.");

  // Check alphabet:
  if  (rootFreqs && model->alphabet().getAlphabetType() != rootFreqs->alphabet().getAlphabetType())
    throw AlphabetMismatchException("NonHomogeneousSubstitutionProcess::createHomogeneousModelSet()", model->getAlphabet().get(), rootFreqs->getAlphabet().get());

  unique_ptr<AutonomousSubstitutionProcessInterface> modelSet;

  if (!rdist)
    modelSet.reset(new SimpleSubstitutionProcess(model, tree));
  else
    modelSet.reset(new RateAcrossSitesSubstitutionProcess(model, rdist, tree));

  if (rootFreqs)
    modelSet->setRootFrequencySet(rootFreqs);

  if (scenario)
    modelSet->setModelScenario(scenario);

  return modelSet;
}

unique_ptr<NonHomogeneousSubstitutionProcess> NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(
    shared_ptr<BranchModelInterface> model,
    shared_ptr<DiscreteDistributionInterface> rdist,
    shared_ptr<PhyloTree> tree,
    shared_ptr<FrequencySetInterface> rootFreqs,
    const vector<string>& globalParameterNames,
    shared_ptr<ModelScenario> scenario)
{
  if (!tree)
    throw Exception("NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess: missing tree.");

  // Check alphabet:
  if (rootFreqs && model->alphabet().getAlphabetType() != rootFreqs->alphabet().getAlphabetType())
    throw AlphabetMismatchException("NonHomogeneousSubstitutionProcess::createNonHomogeneousModelSet()", model->getAlphabet().get(), rootFreqs->getAlphabet().get());
  if (dynamic_pointer_cast<MixedTransitionModelInterface>(model) != nullptr)
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

  auto modelSet = make_unique<NonHomogeneousSubstitutionProcess>(rdist, tree, rootFreqs);

  // We assign a copy of this model to all nodes in the tree (excepted root node), and link all parameters with it.
  vector<unsigned int> ids = tree->getAllNodesIndexes();
  unsigned int rootId = tree->getRootIndex();
  size_t pos = 0;
  for (i = 0; i < ids.size(); ++i)
  {
    if (ids[i] == rootId)
    {
      pos = i;
      break;
    }
  }

  ids.erase(ids.begin() + (long)pos);
  sort(ids.begin(), ids.end());

  for (i = 0; i < ids.size(); i++)
  {
    modelSet->addModel(shared_ptr<BranchModelInterface>(model->clone()), vector<unsigned int>(1, ids[i]));
  }

  // Now alias all global parameters on all nodes:
  for (i = 0; i < globalParameters.size(); i++)
  {
    string pname = globalParameters[i].getName();

    for (size_t nn = 1; nn < ids.size(); nn++)
    {
      modelSet->aliasParameters(pname + "_1", pname + "_" + TextTools::toString(nn + 1));
    }
  }

  if (scenario)
    throw Exception("NonHomogeneousSubstitutionProcess::createNonHomogeneousModelSet : setModelScenario(scenario) to be implemented.");

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
