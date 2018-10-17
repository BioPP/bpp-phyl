//
// File: SubstitutionModelSetTools.cpp
// Created by: Julien Dutheil
// Created on: tue Sep 17 16:57 2007
//

/*
   Copyright or © or Copr. CNRS, (November 16, 2004)

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

#include "SubstitutionModelSetTools.h"
#include "MixedSubstitutionModelSet.h"
#include "MixedSubstitutionModel.h"

using namespace bpp;

using namespace std;

SubstitutionModelSet* SubstitutionModelSetTools::createHomogeneousModelSet(
  TransitionModel* model,
  FrequenciesSet* rootFreqs,
  const Tree* tree
  )
{
  // Check alphabet:
  if (model->getAlphabet()->getAlphabetType() != rootFreqs->getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("SubstitutionModelSetTools::createHomogeneousModelSet()", model->getAlphabet(), rootFreqs->getAlphabet());
  if (dynamic_cast<MixedSubstitutionModel*>(model) != NULL)
    throw Exception("createHomogeneousModelSet non yet programmed for mixture models.");

  SubstitutionModelSet*  modelSet = new SubstitutionModelSet(model->getAlphabet());

  modelSet->setRootFrequencies(rootFreqs);
  // We assign this model to all nodes in the tree (excepted root node), and link all parameters with it.
  vector<int> ids = tree->getNodesId();
  int rootId = tree->getRootId();
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

SubstitutionModelSet* SubstitutionModelSetTools::createNonHomogeneousModelSet(
  TransitionModel* model,
  FrequenciesSet* rootFreqs,
  const Tree* tree,
  const std::map<std::string, std::string>& aliasFreqNames,
  const std::map<std::string, std::vector<Vint> >& globalParameterNames
  )
{
  // Check alphabet:
  if (rootFreqs && model->getAlphabet()->getAlphabetType() != rootFreqs->getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("SubstitutionModelSetTools::createNonHomogeneousModelSet()", model->getAlphabet(), rootFreqs->getAlphabet());
  ParameterList globalParameters;
  globalParameters = model->getParameters();
  
  vector<string> modelParamNames=globalParameters.getParameterNames();
  
  map<string, vector<Vint> > globalParameterNames2;
  
  // First get correct parameter names

  for (const auto& name: globalParameterNames)
  {
    vector<string> complName=ApplicationTools::matchingParameters(name.first, modelParamNames);

    if (complName.size()==0)
      throw Exception("SubstitutionModelSetTools::createNonHomogeneousModelSet. Parameter '" + name.first + "' is not valid.");
    else
      for (size_t j=0; j<complName.size(); j++)
        globalParameterNames2[complName[j]]=name.second;
  }

  SubstitutionModelSet*  modelSet;
  modelSet = new SubstitutionModelSet(model->getAlphabet());

  if (rootFreqs)
    modelSet->setRootFrequencies(rootFreqs);
  
  // We assign a copy of this model to all nodes in the tree (excepted root node), and link all parameters with it.
  vector<int> ids = tree->getNodesId();
  int rootId = tree->getRootId();
  size_t pos = 0;
  for (size_t i = 0; i < ids.size(); i++)
  {
    if (ids[i] == rootId)
    {
      pos = i;
      break;
    }
  }

  ids.erase(ids.begin() + static_cast<ptrdiff_t>(pos));
  for (size_t i = 0; i < ids.size(); i++)
  {
    modelSet->addModel(dynamic_cast<TransitionModel*>(model->clone()), vector<int>(1, ids[i]));
  }

  // Now alias all global parameters on all nodes:
  for (size_t nn=0;nn<globalParameters.size();nn++)
  {
    const Parameter& param=globalParameters[nn];
    
    string pname=param.getName();

    if (globalParameterNames2.find(pname)!=globalParameterNames2.end())
    {
      const vector<Vint>& vvids(globalParameterNames2[pname]);

      if (vvids.size()==0)
      {
        size_t fmid=modelSet->getModelIndexForNode(ids[0])+1;
        for (size_t i = 1; i < ids.size(); i++)
          modelSet->aliasParameters(pname+"_"+TextTools::toString(fmid),pname+"_"+TextTools::toString(modelSet->getModelIndexForNode(ids[i])+1));
      }
      else
        for (const auto& vids:vvids)
        {
          size_t fmid=modelSet->getModelIndexForNode(vids[0]+1);
          for (size_t i=1;i<vids.size();i++)
            modelSet->aliasParameters(pname+"_"+TextTools::toString(fmid),pname+"_"+TextTools::toString(modelSet->getModelIndexForNode(vids[i])+1));
        }
    }
  }
  
  // and alias on the root
  std::map<std::string, std::string>::const_iterator it;

  for (it = aliasFreqNames.begin(); it != aliasFreqNames.end(); it++)
  {
    if (globalParameterNames2.find(it->second)!=globalParameterNames2.end())
      modelSet->aliasParameters(it->second + "_1",it->first);
  }

  delete model; // delete template model.
  return modelSet;
}

