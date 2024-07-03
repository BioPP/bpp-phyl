// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "RateAcrossSitesSubstitutionProcess.h"

using namespace bpp;
using namespace std;

RateAcrossSitesSubstitutionProcess::RateAcrossSitesSubstitutionProcess(
    std::shared_ptr<BranchModelInterface> model,
    std::shared_ptr<DiscreteDistributionInterface> rdist,
    std::shared_ptr<const PhyloTree> tree,
    std::shared_ptr<FrequencySetInterface> rootFrequencies) :
  AbstractParameterAliasable(""),
  AbstractAutonomousSubstitutionProcess(tree, rootFrequencies, model ? model->getNamespace() : ""),
  model_(model),
  rDist_(rdist)
{
  if (!model)
    throw Exception("RateAcrossSitesSubstitutionProcess. A model instance must be provided.");
  if (!rdist)
    throw Exception("RateAcrossSitesSubstitutionProcess. A rate distribution instance must be provided.");

  addParameters_(model->getIndependentParameters()); // Substitution model

  addParameters_(rdist->getIndependentParameters()); // Rate
                                                     // distribution
}

RateAcrossSitesSubstitutionProcess::RateAcrossSitesSubstitutionProcess(
    std::shared_ptr<BranchModelInterface> model,
    std::shared_ptr<DiscreteDistributionInterface> rdist,
    std::shared_ptr<ParametrizablePhyloTree> tree,
    std::shared_ptr<FrequencySetInterface> rootFrequencies) :
  AbstractParameterAliasable(""),
  AbstractAutonomousSubstitutionProcess(tree, rootFrequencies, model ? model->getNamespace() : ""),
  model_(model),
  rDist_(rdist)
{
  if (!model)
    throw Exception("RateAcrossSitesSubstitutionProcess. A model instance must be provided.");
  if (!rdist)
    throw Exception("RateAcrossSitesSubstitutionProcess. A rate distribution instance must be provided.");

  addParameters_(model->getIndependentParameters()); // Substitution model

  addParameters_(rdist->getIndependentParameters()); // Rate
                                                     // distribution
}

RateAcrossSitesSubstitutionProcess::RateAcrossSitesSubstitutionProcess(const RateAcrossSitesSubstitutionProcess& rassp) :
  AbstractParameterAliasable(rassp),
  AbstractAutonomousSubstitutionProcess(rassp),
  model_(rassp.model_->clone()),
  rDist_(rassp.rDist_->clone())
{
  if (modelScenario_)
    modelScenario_->changeModel(std::dynamic_pointer_cast<MixedTransitionModelInterface>(rassp.model_), std::dynamic_pointer_cast<MixedTransitionModelInterface>(model_));
}


RateAcrossSitesSubstitutionProcess& RateAcrossSitesSubstitutionProcess::operator=(const RateAcrossSitesSubstitutionProcess& rassp)
{
  AbstractParameterAliasable::operator=(rassp);
  AbstractAutonomousSubstitutionProcess::operator=(rassp);
  model_.reset(rassp.model_->clone());
  rDist_.reset(rassp.rDist_->clone());

  if (modelScenario_)
    modelScenario_->changeModel(std::dynamic_pointer_cast<MixedTransitionModelInterface>(rassp.model_), std::dynamic_pointer_cast<MixedTransitionModelInterface>(model_));

  return *this;
}

void RateAcrossSitesSubstitutionProcess::setModelScenario(std::shared_ptr<ModelScenario> modelpath)
{
  auto vmod = modelpath->getModels();

  if (vmod.size() == 0) // as if no scenario
    return;

  if (vmod.size() != 1)
    throw Exception("RateAcrossSitesSubstitutionProcess::setModelScenario: model path must have exactly one model.");

  if (vmod[0] != model_)
    throw Exception("RateAcrossSitesSubstitutionProcess::setModelScenario: models are different " + vmod[0]->getNModel(0)->getName() + " != " + model_->getName());

  modelScenario_ = modelpath;
}

void RateAcrossSitesSubstitutionProcess::fireParameterChanged(const ParameterList& pl)
{
  // Update rate distribution:
  rDist_->matchParametersValues(pl);
  model_->matchParametersValues(pl);

  // Transition probabilities have changed and need to be recomputed:
  AbstractAutonomousSubstitutionProcess::fireParameterChanged(pl);
}
