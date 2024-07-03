// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/NumConstants.h>

#include "../FrequencySet/CodonFrequencySet.h"
#include "../Nucleotide/GTR.h"
#include "KCM.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

KCM::KCM(
    std::shared_ptr<const GeneticCode> gc,
    bool oneModel) :
  AbstractParameterAliasable("KCM" + string(oneModel ? "7" : "19") + "."),
  AbstractWrappedModel("KCM" + string(oneModel ? "7" : "19") + "."),
  AbstractWrappedTransitionModel("KCM" + string(oneModel ? "7" : "19") + "."),
  AbstractTotallyWrappedTransitionModel("KCM" + string(oneModel ? "7" : "19") + "."),
  AbstractBiblioTransitionModel("KCM" + string(oneModel ? "7" : "19") + "."),
  AbstractWrappedSubstitutionModel("KCM" + string(oneModel ? "7" : "19") + "."),
  AbstractTotallyWrappedSubstitutionModel("KCM" + string(oneModel ? "7" : "19") + "."),
  AbstractBiblioSubstitutionModel("KCM" + string(oneModel ? "7" : "19") + "."),
  pmodel_(),
  oneModel_(oneModel)
{
  shared_ptr<const NucleicAlphabet> nalph = gc->codonAlphabet().getNucleicAlphabet();

  if (oneModel)
    pmodel_.reset(new KroneckerCodonDistanceSubstitutionModel(
          gc,
          make_unique<GTR>(nalph)));
  else
    pmodel_.reset(new KroneckerCodonDistanceSubstitutionModel(
          gc,
          make_unique<GTR>(nalph),
          make_unique<GTR>(nalph),
          make_unique<GTR>(nalph)));

  string name = "KCM" + string(oneModel ? "7" : "19") + ".";

  pmodel_->setNamespace(name);

  addParameters_(pmodel_->getParameters());

  getParameter_("beta").setName(name + "omega"),

  lParPmodel_.addParameters(pmodel_->getParameters());

  vector<std::string> v = lParPmodel_.getParameterNames();

  for (auto& vi : v)
  {
    mapParNamesFromPmodel_[vi] = getParameterNameWithoutNamespace(vi);
  }

  mapParNamesFromPmodel_[name + "beta"] = "omega";

  updateMatrices_();
}


KCM::KCM(const KCM& kcm) :
  AbstractParameterAliasable(kcm),
  AbstractWrappedModel(kcm),
  AbstractWrappedTransitionModel(kcm),
  AbstractTotallyWrappedTransitionModel(kcm),
  AbstractBiblioTransitionModel(kcm),
  AbstractWrappedSubstitutionModel(kcm),
  AbstractTotallyWrappedSubstitutionModel(kcm),
  AbstractBiblioSubstitutionModel(kcm),
  pmodel_(new KroneckerCodonDistanceSubstitutionModel(*kcm.pmodel_)),
  oneModel_(kcm.oneModel_)
{}

KCM& KCM::operator=(const KCM& kcm)
{
  AbstractBiblioSubstitutionModel::operator=(kcm);

  oneModel_ = kcm.oneModel_;
  pmodel_.reset(new KroneckerCodonDistanceSubstitutionModel(*kcm.pmodel_));
  return *this;
}
