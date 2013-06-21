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
  pSubProColl_(set.pSubProColl),
  nodeToModel_(set.nodeToModel_),
  modelToNodes_(set.modelToNodes_),
  nTree_(set.nTree_),
  nRate_(set.nRate_),
  stationarity_(set.stationarity_),
  nRoot_(set.nRoot_)
{
}

SubstitutionModelSet& SubstitutionModelSet::operator=(const SubstitutionModelSet& set)
{
  pSubProColl_ = set.pSubProColl;
  nodeToModel_ = set.nodeToModel_;
  modelToNodes_ = set.modelToNodes_;
  nTree_ = set.nTree_;
  nRate_ = set.nRate_;
  stationarity_ = set.stationarity_;
  nRoot_ = set.nRoot_;

  return *this;
}

void SubstitutionModelSet::clear()
{
  nodeToModel_.clear();
  modelToNodes_.clear();
  stationarity_=true;
}

void SubstitutionModelSet::addModel(unsigned int numModel, const std::vector<int>& nodesId)
{
  SubstitutionModel* nmod=pSubProColl_->getModel(numModel);

  if (modelToNodes_.size()>0){
    SubstitutionModel* pmodi=pSubProColl_->getModel(modelToNodes_.begin()->first);
    if (nmod->getAlphabet()->getAlphabetType() !=  pmodi->getAlphabet()->getAlphabetType())
      throw Exception("SubstitutionModelSet::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same alphabet.");
    if (nmod->getNumberOfStates() != pmodi->getNumberOfStates())
      throw Exception("SubstitutionModelSet::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same number of states.");
  }
  else if (!stationarity_) {
    FrequenciesSet* pfreq=pSubProColl_->getFrequencies(nRoot_);
    if (pfreq->getAlphabet()->getAlphabetType() != nmod->getAlphabet()->getAlphabetType())
      throw Exception("SubstitutionModelSet::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same alphabet as the root frequencies.");
    if (pfreq->getFrequencies->size() != nmod->getNumberOfStates())
      throw Exception("SubstitutionModelSet::addModel. A Substitution Model cannot be added to a Model Set if it does not have the same number of states as the root frequencies.");
  }
  
  // Associate this model to specified nodes:
  for (size_t i = 0; i < nodesId.size(); i++)
    {
      nodeToModel_[nodesId[i]] = thisModelIndex;
      modelToNodes_[thisModelIndex].push_back(nodesId[i]);
    }
}

void SubstitutionModelSet::addRootFrequencies(unsigned int numFreq)
{
  FrequenciesSet* nfreq=pSubProColl_->getFrequencies(numFreq);
  if (modelToNodes_.size()>0){
    SubstitutionModel* pmodi=pSubProColl_->getModel(modelToNodes_.begin()->first);

    if (pfreq->getAlphabet()->getAlphabetType() != pmodi->getAlphabet()->getAlphabetType())
      throw Exception("SubstitutionModelSet::addRootFrequencies. A Frequencies Set cannot be added to a Model Set if it does not have the same alphabet as the models.");
    if (pfreq->getFrequencies->size() != pmodi->getNumberOfStates())
      throw Exception("SubstitutionModelSet::addRootFrequencies. A Frequencies Set cannot be added to a Model Set if it does not have the same number of states as the models.");
  }
  
  stationarity_=false;
  nRoot_=numFreq;
}

void SubstitutionModelSet::setTree(unsigned int numTree)
{
  nTree_=numTree;
}

void SubstitutionModelSet::setRate(unsigned int numRate)
{
  nRate_=numRate;
}



// void SubstitutionModelSet::removeModel(size_t modelIndex) throw (Exception)
// {
//   modelSet_.erase(modelSet_.begin() + modelIndex);
  
//   // Erase all parameter references to this model and translate other indices...

//   ParameterList pl=getParameters();

//   for (size_t i = pl.size(); i>0; i--)
//     {
//       string pn=pl[i-1].getName();

//       size_t pu=pn.rfind("_");
//       int nm=TextTools::toInt(pn.substr(pu+1,string::npos));
//       if (nm==(int)modelIndex+1){
//         vector<string> alpn=getAlias(pn);
//         for (unsigned j=0; j<alpn.size(); j++)
//           try {
//             unaliasParameters(alpn[j],pn);
//           }
//           catch (Exception& e)
//             {
//               continue;
//             }
//         deleteParameter_(i-1);
//       }
//     }

//   modelParameters_[modelIndex].reset();
// }

// void SubstitutionModelSet::listModelNames(std::ostream& out) const
// {
//   for (size_t i = 0; i < modelSet_.size(); i++)
//     {
//       out << "Model " << i + 1 << ": " << modelSet_[i]->getName() << "\t attached to nodes ";
//       for (size_t j = 0; j < modelToNodes_[i].size(); j++)
//         {
//           out << modelToNodes_[i][j];
//         }
//       out << endl;
//     }
// }


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
  std::map<size_t, std::vector<int> >::const_iterator it;
  for (it= modelToNodes_.begin(); it != modelToNodes_.end() ; it++)
    {
      if (dynamic_cast<const MixedSubstitutionModel*>(pSubProColl_->getModel(it->first)) != NULL)
        return true;
    }
  return false;
}

