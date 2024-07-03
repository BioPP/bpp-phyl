// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>

#include "../FrequencySet/ProteinFrequencySet.h"
#include "../MixtureOfSubstitutionModels.h"
#include "LGL08_CAT.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

LGL08_CAT::LGL08_CAT(
    std::shared_ptr<const ProteicAlphabet> alpha,
    unsigned int nbCat) :
  AbstractParameterAliasable("LGL08_CAT."),
  AbstractWrappedModel("LGL08_CAT."),
  AbstractWrappedTransitionModel("LGL08_CAT."),
  AbstractTotallyWrappedTransitionModel("LGL08_CAT."),
  AbstractBiblioTransitionModel("LGL08_CAT."),
  AbstractBiblioMixedTransitionModel("LGL08_CAT.")
{
  // build the submodel

  vector<unique_ptr<TransitionModelInterface>> vpSM;
  for (unsigned int i = 1; i < nbCat + 1; ++i)
  {
    vpSM.push_back(make_unique<LGL08_CAT::EmbeddedModel>(alpha, "C" + TextTools::toString(i), nbCat));
  }

  Vdouble vrate, vproba;

  for (auto& vi : vpSM)
  {
    vproba.push_back((dynamic_cast<LGL08_CAT::EmbeddedModel&>(*vi)).getProportion());
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
    addParameter_(new Parameter("LGL08_CAT." + st,
          mixedModelPtr_->getParameterValue(st),
          mixedModelPtr_->parameter(st).hasConstraint() ? shared_ptr<ConstraintInterface>(mixedModelPtr_->parameter(st).getConstraint()->clone()) : 0));
  }

  updateMatrices_();
}

/**************** sub model classes */ // ////////

LGL08_CAT::EmbeddedModel::EmbeddedModel(
    std::shared_ptr<const ProteicAlphabet> alpha,
    string name,
    unsigned int nbCat) :
  AbstractParameterAliasable(name),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), name),
  proportion_(1),
  name_(name)
{
  // Exchangeabilities:
  for (unsigned int i = 0; i < 20; ++i)
  {
    for (unsigned int j = 0; j < 20; ++j)
    {
      if (i == j)
        exchangeability_(i, i) = -19.;
      else
        exchangeability_(i, j) = 1.;
    }
  }

  // Equilibrium frequencies, rates and proportions:
  if (nbCat == 10)
  {
#include "__CATC10FrequenciesCode"
#include "__CATC10RatesProps"
  }
  else if (nbCat == 20)
  {
#include "__CATC20FrequenciesCode"
#include "__CATC20RatesProps"
  }
  else if (nbCat == 30)
  {
#include "__CATC30FrequenciesCode"
#include "__CATC30RatesProps"
  }
  else if (nbCat == 40)
  {
#include "__CATC40FrequenciesCode"
#include "__CATC40RatesProps"
  }
  else if (nbCat == 50)
  {
#include "__CATC50FrequenciesCode"
#include "__CATC50RatesProps"
  }
  else if (nbCat == 60)
  {
#include "__CATC60FrequenciesCode"
#include "__CATC60RatesProps"
  }
  else
    throw Exception("LGL08_CAT.cpp: incorrect number of profiles. This number has to be 10, 20, 30, 40, 50 or 60.");

  updateMatrices_();
}
