// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SimpleSubstitutionProcess.h"

using namespace bpp;
using namespace std;

SimpleSubstitutionProcess::SimpleSubstitutionProcess(
    std::shared_ptr<BranchModelInterface> model,
    std::shared_ptr<ParametrizablePhyloTree> tree,
    std::shared_ptr<FrequencySetInterface> rootFrequencies) :
  AbstractParameterAliasable(""),
  AbstractAutonomousSubstitutionProcess(tree, rootFrequencies, model ? model->getNamespace() : ""),
  model_(model)
{
  if (!model)
    throw Exception("SimpleSubstitutionProcess. A model instance must be provided.");

  // Add parameters:
  addParameters_(model->getIndependentParameters()); // Substitution model
}

SimpleSubstitutionProcess::SimpleSubstitutionProcess(
    std::shared_ptr<BranchModelInterface> model,
    std::shared_ptr<const PhyloTree> tree,
    std::shared_ptr<FrequencySetInterface> rootFrequencies) :
  AbstractParameterAliasable(""),
  AbstractAutonomousSubstitutionProcess(tree, rootFrequencies, model ? model->getNamespace() : ""),
  model_(model)
{
  if (!model)
    throw Exception("SimpleSubstitutionProcess. A model instance must be provided.");

  // Add parameters:
  addParameters_(model->getIndependentParameters()); // Substitution model
}

SimpleSubstitutionProcess::SimpleSubstitutionProcess(const SimpleSubstitutionProcess& ssp) :
  AbstractParameterAliasable(ssp),
  AbstractAutonomousSubstitutionProcess(ssp),
  model_(ssp.model_->clone())
{
  if (modelScenario_)
    modelScenario_->changeModel(dynamic_pointer_cast<MixedTransitionModelInterface>(ssp.model_), dynamic_pointer_cast<MixedTransitionModelInterface>(model_));
}

SimpleSubstitutionProcess& SimpleSubstitutionProcess::operator=(const SimpleSubstitutionProcess& ssp)
{
  AbstractParameterAliasable::operator=(ssp);
  AbstractAutonomousSubstitutionProcess::operator=(ssp);
  model_.reset(ssp.model_->clone());

  if (modelScenario_)
    modelScenario_->changeModel(dynamic_pointer_cast<MixedTransitionModelInterface>(ssp.model_), dynamic_pointer_cast<MixedTransitionModelInterface>(model_));

  return *this;
}

void SimpleSubstitutionProcess::setModelScenario(std::shared_ptr<ModelScenario> modelpath)
{
  auto vmod = modelpath->getModels();

  if (vmod.size() == 0) // as if no scenario
    return;

  if (vmod.size() != 1)
    throw Exception("SimpleSubstitutionProcess::setModelScenario: model path must have exactly one model.");

  auto mixed_ = dynamic_pointer_cast<MixedTransitionModelInterface>(model_);
  if (!mixed_)
    throw Exception("SimpleSubstitutionProcess::setModelScenario: model must be mixed.");

  if (vmod[0] != mixed_)
    throw Exception("SimpleSubstitutionProcess::setModelScenario: models are different " + vmod[0]->getNModel(0)->getName() + " != " + mixed_->getNModel(0)->getName());

  modelScenario_ = modelpath;
}

void SimpleSubstitutionProcess::fireParameterChanged(const ParameterList& pl)
{
  // Transition probabilities have changed and need to be recomputed:
  AbstractAutonomousSubstitutionProcess::fireParameterChanged(pl);

  model_->matchParametersValues(pl);
}
