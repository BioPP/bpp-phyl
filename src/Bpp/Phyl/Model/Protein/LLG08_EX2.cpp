// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>

#include "../FrequencySet/ProteinFrequencySet.h"
#include "../MixtureOfSubstitutionModels.h"
#include "LLG08_EX2.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

LLG08_EX2::LLG08_EX2(shared_ptr<const ProteicAlphabet> alpha) :
  AbstractParameterAliasable("LLG08_EX2."),
  AbstractWrappedModel("LLG08_EX2."),
  AbstractWrappedTransitionModel("LLG08_EX2."),
  AbstractTotallyWrappedTransitionModel("LLG08_EX2."),
  AbstractBiblioTransitionModel("LLG08_EX2."),
  AbstractBiblioMixedTransitionModel("LLG08_EX2.")
{
  // build the submodel

  vector<unique_ptr<TransitionModelInterface>> vpSM;
  vpSM.push_back(make_unique<LLG08_EX2::EmbeddedModel>(alpha, "Buried"));
  vpSM.push_back(make_unique<LLG08_EX2::EmbeddedModel>(alpha, "Exposed"));

  Vdouble vrate, vproba;

  for (auto& vi : vpSM)
  {
    vproba.push_back(dynamic_cast<LLG08_EX2::EmbeddedModel&>(*vi).getProportion());
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
    addParameter_(new Parameter("LLG08_EX2." + st,
          mixedModelPtr_->getParameterValue(st),
          mixedModelPtr_->parameter(st).hasConstraint() ? std::shared_ptr<ConstraintInterface>(mixedModelPtr_->parameter(st).getConstraint()->clone()) : 0));
  }

  updateMatrices_();
}

/**************** sub model classes *******************/

LLG08_EX2::EmbeddedModel::EmbeddedModel(
    shared_ptr<const ProteicAlphabet> alpha,
    string name) :
  AbstractParameterAliasable(name),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), name),
  proportion_(1),
  name_(name)
{
#include "__LLG08_EX2ExchangeabilityCode"
#include "__LLG08_EX2FrequenciesCode"
#include "__LLG08_EX2RatesProps"
  updateMatrices_();
}
