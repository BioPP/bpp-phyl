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
  SubstitutionModel* model,
  FrequenciesSet* rootFreqs,
  const Tree* tree
  ) throw (AlphabetException, Exception)
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
  modelSet->addModel(model, ids, model->getParameters().getParameterNames());

  return modelSet;
}

SubstitutionModelSet* SubstitutionModelSetTools::createNonHomogeneousModelSet(
  SubstitutionModel* model,
  FrequenciesSet* rootFreqs,
  const Tree* tree,
  const vector<string>& globalParameterNames
  ) throw (AlphabetException, Exception)
{
  // Check alphabet:
  if (rootFreqs && model->getAlphabet()->getAlphabetType() != rootFreqs->getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("SubstitutionModelSetTools::createNonHomogeneousModelSet()", model->getAlphabet(), rootFreqs->getAlphabet());
  ParameterList globalParameters, branchParameters;
  globalParameters = model->getParameters();

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
      throw Exception("SubstitutionModelSetTools::createNonHomogeneousModelSet. Parameter '" + globalParameterNames[i] + "' is not valid.");
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

  bool mixed = (dynamic_cast<MixedSubstitutionModel*>(model) != NULL);
  SubstitutionModelSet*  modelSet;
  if (mixed)
  {
    modelSet = new MixedSubstitutionModelSet(model->getAlphabet());
    // Remove the "relproba" parameters from the branch parameters and put them in the global parameters, for the hypernodes
    for (i = branchParameters.size(); i > 0; i--)
    {
      if (branchParameters[i - 1].getName().find("relproba") != string::npos)
      {
        globalParameters.addParameter(branchParameters[i - 1]);
        branchParameters.deleteParameter(i - 1);
      }
    }
  }
  else
    modelSet = new SubstitutionModelSet(model->getAlphabet());

  if (rootFreqs)
    modelSet->setRootFrequencies(rootFreqs);

  // We assign a copy of this model to all nodes in the tree (excepted root node), and link all parameters with it.
  vector<int> ids = tree->getNodesId();
  int rootId = tree->getRootId();
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
    modelSet->addModel(dynamic_cast<SubstitutionModel*>(model->clone()), vector<int>(1, ids[i]), branchParameters.getParameterNames());
  }

  // Now add global parameters to all nodes:
  modelSet->addParameters(globalParameters, ids);

  // Defines the hypernodes if mixed
  if (mixed)
  {
    MixedSubstitutionModelSet* pMSMS = dynamic_cast<MixedSubstitutionModelSet*>(modelSet);
    MixedSubstitutionModel* pMSM = dynamic_cast<MixedSubstitutionModel*>(model);

    size_t nbm = pMSM->getNumberOfModels();
    for (i = 0; i < nbm; i++)
    {
      pMSMS->addEmptyHyperNode();
      for (j = 0; j < ids.size(); j++)
      {
        pMSMS->addToHyperNode(j, vector<int>(1, static_cast<int>(i)));
      }
    }
    pMSMS->computeHyperNodesProbabilities();
  }

  delete model; // delete template model.
  return modelSet;
}

