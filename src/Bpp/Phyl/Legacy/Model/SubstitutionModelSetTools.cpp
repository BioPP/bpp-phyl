// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../../Model/MixedTransitionModel.h"
#include "MixedSubstitutionModelSet.h"
#include "SubstitutionModelSetTools.h"

using namespace bpp;

using namespace std;

unique_ptr<SubstitutionModelSet> SubstitutionModelSetTools::createHomogeneousModelSet(
    shared_ptr<TransitionModelInterface> model,
    shared_ptr<FrequencySetInterface> rootFreqs,
    const Tree& tree
    )
{
  // Check alphabet:
  if (model->getAlphabet()->getAlphabetType() != rootFreqs->getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("SubstitutionModelSetTools::createHomogeneousModelSet()", model->getAlphabet(), rootFreqs->getAlphabet());
  if (dynamic_cast<MixedTransitionModelInterface*>(model.get()) != nullptr)
    throw Exception("createHomogeneousModelSet non yet programmed for mixture models.");

  auto modelSet = make_unique<SubstitutionModelSet>(model->getAlphabet());

  modelSet->setRootFrequencies(rootFreqs);
  // We assign this model to all nodes in the tree (excepted root node), and link all parameters with it.
  vector<int> ids = tree.getNodesId();
  int rootId = tree.getRootId();
  unsigned int pos = 0;
  for (unsigned int i = 0; i < ids.size(); ++i)
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

unique_ptr<SubstitutionModelSet> SubstitutionModelSetTools::createNonHomogeneousModelSet(
    shared_ptr<TransitionModelInterface> model,
    shared_ptr<FrequencySetInterface> rootFreqs,
    const Tree& tree,
    const map<string, string>& aliasFreqNames,
    const map<string, vector<Vint>>& globalParameterNames
    )
{
  // Check alphabet:
  if (rootFreqs && model->alphabet().getAlphabetType() != rootFreqs->alphabet().getAlphabetType())
    throw AlphabetMismatchException("SubstitutionModelSetTools::createNonHomogeneousModelSet()", model->getAlphabet(), rootFreqs->getAlphabet());
  ParameterList globalParameters;
  globalParameters = model->getParameters();

  vector<string> modelParamNames = globalParameters.getParameterNames();

  map<string, vector<Vint>> globalParameterNames2;

  // First get correct parameter names

  for (const auto& name : globalParameterNames)
  {
    vector<string> complName = ApplicationTools::matchingParameters(name.first, modelParamNames);

    if (complName.size() == 0)
      throw Exception("SubstitutionModelSetTools::createNonHomogeneousModelSet. Parameter '" + name.first + "' is not valid.");
    else
      for (auto& cni : complName)
      {
        globalParameterNames2[cni] = name.second;
      }
  }

  auto modelSet = make_unique<SubstitutionModelSet>(model->getAlphabet());

  if (rootFreqs)
    modelSet->setRootFrequencies(rootFreqs);

  // We assign a copy of this model to all nodes in the tree (excepted root node), and link all parameters with it.
  vector<int> ids = tree.getNodesId();
  int rootId = tree.getRootId();
  size_t pos = 0;
  for (size_t i = 0; i < ids.size(); ++i)
  {
    if (ids[i] == rootId)
    {
      pos = i;
      break;
    }
  }

  ids.erase(ids.begin() + static_cast<ptrdiff_t>(pos));
  for (size_t i = 0; i < ids.size(); ++i)
  {
    modelSet->addModel(shared_ptr<TransitionModelInterface>(model->clone()), vector<int>(1, ids[i]));
  }

  // Now alias all global parameters on all nodes:
  for (size_t nn = 0; nn < globalParameters.size(); ++nn)
  {
    const Parameter& param = globalParameters[nn];

    string pname = param.getName();

    if (globalParameterNames2.find(pname) != globalParameterNames2.end())
    {
      const vector<Vint>& vvids(globalParameterNames2[pname]);

      if (vvids.size() == 0)
      {
        size_t fmid = modelSet->getModelIndexForNode(ids[0]) + 1;
        for (size_t i = 1; i < ids.size(); ++i)
        {
          modelSet->aliasParameters(pname + "_" + TextTools::toString(fmid), pname + "_" + TextTools::toString(modelSet->getModelIndexForNode(ids[i]) + 1));
        }
      }
      else
        for (const auto& vids : vvids)
        {
          size_t fmid = modelSet->getModelIndexForNode(vids[0] + 1);
          for (size_t i = 1; i < vids.size(); i++)
          {
            modelSet->aliasParameters(pname + "_" + TextTools::toString(fmid), pname + "_" + TextTools::toString(modelSet->getModelIndexForNode(vids[i]) + 1));
          }
        }
    }
  }

  // and alias on the root
  for (auto& it : aliasFreqNames)
  {
    if (globalParameterNames2.find(it.second) != globalParameterNames2.end())
      modelSet->aliasParameters(it.second + "_1", it.first);
  }

  return modelSet;
}
