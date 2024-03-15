// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../FrequencySet/CodonFrequencySet.h"
#include "../Nucleotide/K80.h"
#include "GY94.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

GY94::GY94(
    shared_ptr<const GeneticCode> gc,
    unique_ptr<CodonFrequencySetInterface> codonFreqs) :
  AbstractParameterAliasable("GY94."),
  AbstractWrappedModel("GY94."),
  AbstractWrappedTransitionModel("GY94."),
  AbstractTotallyWrappedTransitionModel("GY94."),
  AbstractBiblioTransitionModel("GY94."),
  AbstractWrappedSubstitutionModel("GY94."),
  AbstractTotallyWrappedSubstitutionModel("GY94."),
  AbstractBiblioSubstitutionModel("GY94."),
  gacd_(new GranthamAAChemicalDistance()),
  pmodel_(new CodonDistanceFrequenciesSubstitutionModel(
      gc,
      make_unique<K80>(gc->codonAlphabet().getNucleicAlphabet()), 
      move(codonFreqs),
      gacd_))
{
  addParameter_(new Parameter("GY94.kappa", 1, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("GY94.V", 10000, Parameter::R_PLUS_STAR));

  pmodel_->setNamespace("GY94.");
  addParameters_(pmodel_->frequencySet().getParameters());

  lParPmodel_.addParameters(pmodel_->getParameters());

  vector<std::string> v = pmodel_->frequencySet().getParameters().getParameterNames();
  for (unsigned int i = 0; i < v.size(); i++)
  {
    mapParNamesFromPmodel_[v[i]] = getParameterNameWithoutNamespace(v[i]);
  }

  mapParNamesFromPmodel_["GY94.123_K80.kappa"] = "kappa";
  mapParNamesFromPmodel_["GY94.alpha"] = "V";

  updateMatrices_();
}

GY94::GY94(const GY94& gy94) :
  AbstractParameterAliasable(gy94),
  AbstractWrappedModel(gy94),
  AbstractWrappedTransitionModel(gy94),
  AbstractTotallyWrappedTransitionModel(gy94),
  AbstractBiblioTransitionModel(gy94),
  AbstractWrappedSubstitutionModel(gy94),
  AbstractTotallyWrappedSubstitutionModel(gy94),
  AbstractBiblioSubstitutionModel(gy94),
  gacd_(),
  pmodel_(new CodonDistanceFrequenciesSubstitutionModel(*gy94.pmodel_))
{}

GY94& GY94::operator=(const GY94& gy94)
{
  AbstractBiblioSubstitutionModel::operator=(gy94);
  pmodel_.reset(new CodonDistanceFrequenciesSubstitutionModel(*gy94.pmodel_));
  return *this;
}

GY94::~GY94() {}
