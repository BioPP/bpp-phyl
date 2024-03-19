// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>

#include "../FrequencySet/ProteinFrequencySet.h"
#include "../MixtureOfSubstitutionModels.h"
#include "LLG08_EHO.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

LLG08_EHO::LLG08_EHO(
    shared_ptr<const ProteicAlphabet> alpha) :
  AbstractParameterAliasable("LLG08_EHO."),
  AbstractWrappedModel("LLG08_EHO."),
  AbstractWrappedTransitionModel("LLG08_EHO."),
  AbstractTotallyWrappedTransitionModel("LLG08_EHO."),
  AbstractBiblioTransitionModel("LLG08_EHO."),
  AbstractBiblioMixedTransitionModel("LLG08_EHO.")
{
  // build the submodel

  vector<unique_ptr<TransitionModelInterface>> vpSM;
  vpSM.push_back(make_unique<LLG08_EHO::EmbeddedModel>(alpha, "Extended"));
  vpSM.push_back(make_unique<LLG08_EHO::EmbeddedModel>(alpha, "Helix"));
  vpSM.push_back(make_unique<LLG08_EHO::EmbeddedModel>(alpha, "Other"));

  Vdouble vrate, vproba;

  for (auto& vi : vpSM)
  {
    vproba.push_back(dynamic_cast<LLG08_EHO::EmbeddedModel&>(*vi).getProportion());
    vrate.push_back(vi->getRate());
  }

  mixedModelPtr_.reset(new MixtureOfSubstitutionModels(alpha, vpSM, vproba, vrate));

  string name, st;
  ParameterList pl = mixedModelPtr_->getParameters();
  for (size_t i = 0; i < pl.size(); ++i)
  {
    name = pl[i].getName();
    lParPmodel_.addParameter(Parameter(pl[i]));
    st = mixedModelPtr_->getParameterNameWithoutNamespace(name);
    mapParNamesFromPmodel_[name] = st;
    addParameter_(new Parameter("LLG08_EHO." + st,
          mixedModelPtr_->getParameterValue(st),
          mixedModelPtr_->parameter(st).hasConstraint() ? shared_ptr<ConstraintInterface>(mixedModelPtr_->parameter(st).getConstraint()->clone()) : nullptr));
  }

  updateMatrices_();
}

/**************** sub model classes ***********************/

LLG08_EHO::EmbeddedModel::EmbeddedModel(
    shared_ptr<const ProteicAlphabet> alpha,
    string name) :
  AbstractParameterAliasable(name),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), name),
  proportion_(1),
  name_(name)
{
#include "__LLG08_EHOExchangeabilityCode"
#include "__LLG08_EHOFrequenciesCode"
#include "__LLG08_EHORatesProps"
  updateMatrices_();
}
