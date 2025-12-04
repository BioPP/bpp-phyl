// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <numeric>

#include "ModelScenario.h"

using namespace bpp;
using namespace std;

bool ModelScenario::complete()
{
  ModelPath nhn;
  for (const auto& mp : vModelPaths_)
  {
    nhn += *mp;
  }

  auto rest = make_shared<ModelPath> ();
  auto models = nhn.getModels();
  for (const auto& model:models)
  {
    Vuint v((unsigned int)model->getNumberOfModels());
    iota(v.begin(), v.end(), 0);
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

void ModelScenario::changeModel(shared_ptr<MixedTransitionModelInterface> m1,
    shared_ptr<MixedTransitionModelInterface> m2)
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

  shared_ptr<MixedTransitionModelInterface> pfSM(vModelPaths_[0]->getLeadModel());

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

string ModelScenario::toString() const
{
  string output;

  for (const auto& mp :vModelPaths_)
  {
    output += "<" + mp->toString() + ">";
  }

  return output;
}

vector<shared_ptr<MixedTransitionModelInterface>> ModelScenario::getModels() const
{
  vector<shared_ptr<MixedTransitionModelInterface>> models, models2;

  for (const auto& mp : vModelPaths_)
  {
    auto vmodel = mp->getModels();
    for (auto& model:vmodel)
    {
      if (find(models.begin(), models.end(), model) == models.end())
        models.push_back(model);
      else if (find(models2.begin(), models2.end(), model) == models2.end())
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
