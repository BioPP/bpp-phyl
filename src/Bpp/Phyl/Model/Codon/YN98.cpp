// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/NumConstants.h>

#include "../FrequencySet/CodonFrequencySet.h"
#include "../Nucleotide/K80.h"
#include "YN98.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

YN98::YN98(
    shared_ptr<const GeneticCode> gc,
    unique_ptr<CodonFrequencySetInterface> codonFreqs) :
  AbstractParameterAliasable("YN98."),
  AbstractWrappedModel("YN98."),
  AbstractWrappedTransitionModel("YN98."),
  AbstractTotallyWrappedTransitionModel("YN98."),
  AbstractBiblioTransitionModel("YN98."),
  AbstractWrappedSubstitutionModel("YN98."),
  AbstractTotallyWrappedSubstitutionModel("YN98."),
  AbstractBiblioSubstitutionModel("YN98."),
  pmodel_(new CodonDistanceFrequenciesSubstitutionModel(gc, make_unique<K80>(gc->codonAlphabet().getNucleicAlphabet()), std::move(codonFreqs)))
{
  computeFrequencies(false);

  addParameter_(new Parameter("YN98.kappa", 1, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("YN98.omega", 1, make_shared<IntervalConstraint>(0.001, 999, true, true)));

  pmodel_->setNamespace("YN98.");
  addParameters_(pmodel_->codonFrequencySet().getParameters());

  lParPmodel_.addParameters(pmodel_->getParameters());

  vector<std::string> v = pmodel_->codonFrequencySet().getParameters().getParameterNames();

  for (auto& vi : v)
  {
    mapParNamesFromPmodel_[vi] = getParameterNameWithoutNamespace(vi);
  }
  mapParNamesFromPmodel_["YN98.123_K80.kappa"] = "kappa";
  mapParNamesFromPmodel_["YN98.beta"] = "omega";

  updateMatrices_();
}


YN98::YN98(const YN98& yn98) :
  AbstractParameterAliasable(yn98),
  AbstractWrappedModel(yn98),
  AbstractWrappedTransitionModel(yn98),
  AbstractTotallyWrappedTransitionModel(yn98),
  AbstractBiblioTransitionModel(yn98),
  AbstractWrappedSubstitutionModel(yn98),
  AbstractTotallyWrappedSubstitutionModel(yn98),
  AbstractBiblioSubstitutionModel(yn98),
  pmodel_(yn98.pmodel_->clone())
{}

YN98& YN98::operator=(const YN98& yn98)
{
  AbstractBiblioSubstitutionModel::operator=(yn98);
  pmodel_.reset(yn98.pmodel_->clone());
  return *this;
}
