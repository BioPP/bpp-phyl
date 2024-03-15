// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>

#include "../FrequencySet/ProteinFrequencySet.h"
#include "../MixtureOfSubstitutionModels.h"
#include "LG10_EX_EHO.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

LG10_EX_EHO::LG10_EX_EHO(shared_ptr<const ProteicAlphabet> alpha) :
  AbstractParameterAliasable("LG10_EX_EHO."),
  AbstractWrappedModel("LG10_EX_EHO."),
  AbstractWrappedTransitionModel("LG10_EX_EHO."),
  AbstractTotallyWrappedTransitionModel("LG10_EX_EHO."),
  AbstractBiblioTransitionModel("LG10_EX_EHO."),
  AbstractBiblioMixedTransitionModel("LG10_EX_EHO.")
{
  // build the submodel

  vector<unique_ptr<TransitionModelInterface>> vpSM;
  vpSM.push_back(make_unique<LG10_EX_EHO::EmbeddedModel>(alpha, "BUR_EXT"));
  vpSM.push_back(make_unique<LG10_EX_EHO::EmbeddedModel>(alpha, "BUR_HEL"));
  vpSM.push_back(make_unique<LG10_EX_EHO::EmbeddedModel>(alpha, "BUR_OTH"));
  vpSM.push_back(make_unique<LG10_EX_EHO::EmbeddedModel>(alpha, "EXP_EXT"));
  vpSM.push_back(make_unique<LG10_EX_EHO::EmbeddedModel>(alpha, "EXP_HEL"));
  vpSM.push_back(make_unique<LG10_EX_EHO::EmbeddedModel>(alpha, "EXP_OTH"));

  Vdouble vrate, vproba;

  for (auto& vi : vpSM)
  {
    vproba.push_back((dynamic_cast<LG10_EX_EHO::EmbeddedModel&>(*vi)).getProportion());
    vrate.push_back(vi->getRate());
  }

  mixedModelPtr_.reset(new MixtureOfSubstitutionModels(alpha, vpSM, vproba, vrate));

  string name, st;
  ParameterList pl = mixedModelPtr_->getParameters();
  for (unsigned int i = 0; i < pl.size(); ++i)
  {
    name = pl[i].getName();
    lParPmodel_.addParameter(Parameter(pl[i]));
    st = mixedModelPtr_->getParameterNameWithoutNamespace(name);
    mapParNamesFromPmodel_[name] = st;
    addParameter_(new Parameter("LG10_EX_EHO." + st,
                                mixedModelPtr_->getParameterValue(st),
                                mixedModelPtr_->parameter(st).hasConstraint() ? std::shared_ptr<ConstraintInterface>(mixedModelPtr_->parameter(st).getConstraint()->clone()) : 0));
  }

  updateMatrices_();
}

/**************** sub model classes ********************/

LG10_EX_EHO::EmbeddedModel::EmbeddedModel(
    shared_ptr<const ProteicAlphabet> alpha,
    string name) :
  AbstractParameterAliasable(name),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), name),
  proportion_(1),
  name_(name)
{
#include "__LG10_EX_EHOExchangeabilityCode"
#include "__LG10_EX_EHOFrequenciesCode"
#include "__LG10_EX_EHORatesProps"
  updateMatrices_();
}
