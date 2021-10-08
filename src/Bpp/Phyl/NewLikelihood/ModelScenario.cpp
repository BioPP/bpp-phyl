//
// File: ModelScenario.cpp
// Created by: Laurent Guéguen
// Created on: mercredi 27 novembre 2019, à 09h 09
// From: MixedSubstitutionModelSet
//

/*
   Copyright or <A9> or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "ModelScenario.h"
#include <numeric>

using namespace bpp;
using namespace std;

bool ModelScenario::complete()
{
  ModelPath nhn;
  for (const auto& mp : vModelPaths_)
  {
    nhn += *mp;
  }

  auto rest = std::make_shared<ModelPath> ();
  auto models = nhn.getModels();
  for (const auto& model:models)
  {
    Vuint v((uint)model->getNumberOfModels());
    std::iota(v.begin(), v.end(), 0);
    rest->setModel(model, v);
  }

  (*rest) -= nhn;

  if (rest->size() != 0)
  {
    addModelPath(rest);
    return true;
  }

  return false;
}

void ModelScenario::changeModel(std::shared_ptr<MixedTransitionModel> m1,
                                std::shared_ptr<MixedTransitionModel> m2)
{
  for (auto& mp: vModelPaths_)
  {
    mp->changeModel(m1, m2);
  }
}

bool ModelScenario::hasExclusivePaths() const
{
  ModelPath tthn;

  size_t nhn = getNumberOfModelPaths();
  for (size_t i = 0; i < nhn; i++)
  {
    if (tthn.intersects(*getModelPath(i)))
      return false;
    else
      tthn += (*getModelPath(i));
  }

  return true;
}

// void ModelScenario::fireParameterChanged(const ParameterList& parameters)
// {
//   SubstitutionModelSet::fireParameterChanged(parameters);

//   // should be restricted only when probability related parameters are changed
//   computeModelPathsProbabilities();
// }


void ModelScenario::computeModelPathsProbabilities()
{
  return;

  size_t nbh = getNumberOfModelPaths();

  if (nbh == 0)
    return;

  // Compute the probabilities of the hypernodes from the lead mixed
  // model of the first ModelPath

  std::shared_ptr<MixedTransitionModel> pfSM(vModelPaths_[0]->getLeadModel());

  if (pfSM == 0)
    throw Exception("ModelScenario::computeModelPathsProbabilities: missing lead Model.");

  for (size_t nh = 0; nh < nbh; nh++)
  {
    ModelPath& h = *getModelPath(nh);
    if (h.hasModel(pfSM))
    {
      const ModelPath::PathNode& fnd = h.getPathNode(pfSM);

      double fprob = 0;
      for (const auto& fn:fnd)
      {
        fprob += pfSM->getNProbability(static_cast<size_t>(fn));
      }

      h.setProbability(fprob);
    }
    else
      throw Exception("ModelScenario::computeModelPathsProbabilities : reference model " + pfSM->getName() + " is missing in ModelPath " + TextTools::toString(nh));
  }

  // Sets the new probabilities & rates of the mixmodels

  auto models = getModels();

  for (auto model:getModels())
  {
    if (model != pfSM)
    {
      for (size_t nh = 0; nh < nbh; nh++)
      {
        ModelPath& h = *getModelPath(nh);
        if (!h.hasModel(model))
          throw Exception("ModelScenario::computeModelPathsProbabilities : reference model " + model->getName() + " is missing in ModelPath " + TextTools::toString(nh));

        const ModelPath::PathNode& fnd = h.getPathNode(model);
        double prob = 0;
        for (auto& fn:fnd)
        {
          prob += model->getNProbability(static_cast<size_t>(fn));
        }

        // sets the real probabilities
        for (auto& fn:fnd)
        {
          model->setNProbability(static_cast<size_t>(fn), h.getProbability() * model->getNProbability(static_cast<size_t>(fn)) / prob);
        }
      }

      // normalizes Vrates with the real probabilities

      model->normalizeVRates();
    }

    // sets the conditional probabilities

    for (size_t nh = 0; nh < nbh; nh++)
    {
      ModelPath& h = *getModelPath(nh);
      const ModelPath::PathNode& fnd = h.getPathNode(model);
      for (auto& fn:fnd)
      {
        model->setNProbability(static_cast<size_t>(fn),
                               model->getNProbability(static_cast<size_t>(fn)) / h.getProbability());
      }
    }
  }
}

std::string ModelScenario::to_string() const
{
  string output;

  for (const auto& mp:vModelPaths_)
  {
    output += "<" + mp->to_string() + ">";
  }

  return output;
}

std::vector<std::shared_ptr<MixedTransitionModel> > ModelScenario::getModels() const
{
  std::vector<std::shared_ptr<MixedTransitionModel> > models, models2;

  for (const auto& mp:vModelPaths_)
  {
    auto vmodel = mp->getModels();
    for (auto& model:vmodel)
    {
      if (std::find(models.begin(), models.end(), model) == models.end())
        models.push_back(model);
      else if (std::find(models2.begin(), models2.end(), model) == models2.end())
        models2.push_back(model);
    }
  }

  return models2; // return only models found in several paths
}

// double ModelScenario::getModelPathProbability(const ModelPath& hn) const
// {
//   auto models= hn.getModels();

//   double fprob = 1;
//   for (auto& model:models)
//   {
//     double x = 0;
//     for (const auto fn:hn.getPathNode(model))
//       x += model->getNProbability(static_cast<size_t>(fn));

//     fprob *= x;
//   }

//   return fprob;
// }
