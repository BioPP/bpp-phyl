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
  if (!(hasParameter(name)))
    throw ParameterNotFoundException("SubstitutionModelSet::getNodesWithParameter.", name);
    
  vector<string> nalias=getAlias(name);
  size_t p=name.rfind("_");
  vector<int> inode = getNodesWithModel(TextTools::to<size_t>(name.substr(p+1,string::npos)) - 1);
  
  for (size_t i = 0; i < nalias.size(); i++)
    {
      p=nalias[i].rfind("_");
      size_t pos=TextTools::to<size_t>(nalias[i].substr(p+1,string::npos));
      if (pos>0){
        vector<int> ni = getNodesWithModel(pos-1);
        inode.insert(inode.end(),ni.begin(),ni.end());
      }
    }
  
  return inode;
}

void SubstitutionModelSet::addModel(SubstitutionModel* model, const std::vector<int>& nodesId)//, const vector<string>& newParams) throw (Exception)
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

  vector<string> nplm=model->getParameters().getParameterNames();
  
  modelParameters_.push_back(*model->getParameters().clone());
  
  for (size_t i  = 0; i < nplm.size(); i++)
    {
      pname = nplm[i];
      Parameter* p = new Parameter(model->getParameters().getParameter(pname)); // We work with namespaces here, so model->getParameter(pname) does not work.
      p->setName(pname + "_" + TextTools::toString(thisModelIndex+1));
      addParameter_(p);
    }
}

void SubstitutionModelSet::replaceModel(size_t modelIndex, SubstitutionModel* model) throw (Exception)
{
  delete modelSet_[modelIndex];
  modelSet_[modelIndex]=model;
  
  // Erase all parameter references to this model

  ParameterList pl=getNodeParameters();

  for (size_t i = pl.size(); i>0; i--)
    {
      string pn=pl[i-1].getName();
      
      size_t pu=pn.rfind("_");
      int nm=TextTools::toInt(pn.substr(pu+1,string::npos));

      if (nm==(int)modelIndex+1){
        vector<string> alpn=getAlias(pn);
        for (unsigned j=0; j<alpn.size(); j++)
          try {
            unaliasParameters(alpn[j],pn);
          }
          catch (Exception& e)
            {
              continue;
            }
        deleteParameter_(pn);
      }
    }

  // Associate new parameters
  string pname;

  vector<string> nplm=model->getParameters().getParameterNames();
  
  for (size_t i  = 0; i < nplm.size(); i++)
  {
    pname = nplm[i];
    Parameter* p = new Parameter(model->getParameters().getParameter(pname)); // We work with namespaces here, so model->getParameter(pname) does not work.
    p->setName(pname + "_" + TextTools::toString(modelIndex+1));
    addParameter_(p);
  }

  // update modelParameters_

  modelParameters_[modelIndex].reset();
  modelParameters_[modelIndex]=*model->getParameters().clone();
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
  AbstractParameterAliasable::fireParameterChanged(parameters);
  
  // Update root frequencies:
  updateRootFrequencies();

  // Then we update all models in the set:
  for (size_t i = 0; i < modelParameters_.size(); i++)
    {
      for (size_t np = 0 ; np< modelParameters_[i].size() ; np++)
        {
          modelParameters_[i][np].setValue(getParameterValue(modelParameters_[i][np].getName()+"_"+TextTools::toString(i+1)));
        }
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

