// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <string>

#include "../FrequencySet/CodonFrequencySet.h"
#include "../Nucleotide/GTR.h"
#include "MG94xREV.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

MG94xREV::MG94xREV(
    std::shared_ptr<const GeneticCode> gc,
    std::unique_ptr<CodonFrequencySetInterface> codonFreqs) :
  AbstractParameterAliasable("MG94xREV."),
  AbstractWrappedModel("MG94xREV."),
  AbstractWrappedTransitionModel("MG94xREV."),
  AbstractTotallyWrappedTransitionModel("MG94xREV."),
  AbstractBiblioTransitionModel("MG94xREV."),
  AbstractWrappedSubstitutionModel("MG94xREV."),
  AbstractTotallyWrappedSubstitutionModel("MG94xREV."),
  AbstractBiblioSubstitutionModel("MG94xREV."),
  pmodel_(new CodonDistancePhaseFrequenciesSubstitutionModel(
        gc,
        make_unique<GTR>(gc->codonAlphabet().getNucleicAlphabet()),
        std::move(codonFreqs)))
{
  computeFrequencies(false);

  addParameter_(new Parameter("MG94xREV.a", 1, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("MG94xREV.b", 1, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("MG94xREV.c", 1, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("MG94xREV.d", 1, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("MG94xREV.e", 1, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("MG94xREV.omega", 1, std::make_shared<IntervalConstraint>(0.002, 999, true, true)));

  pmodel_->setNamespace("MG94xREV.");
  addParameters_(pmodel_->frequencySet().getParameters());

  // get only non frequency parameters from GTR
  auto& pl = pmodel_->getParameters();
  for (auto& name : pl.getParameterNames())
  {
    auto param=pl.parameter(name);
    if (param.getName().substr(0,5)!="theta")
      lParPmodel_.addParameter(param);
  }
  
  vector<std::string> v = pmodel_->frequencySet().getParameters().getParameterNames();
  for (unsigned int i = 0; i < v.size(); ++i)
  {
    mapParNamesFromPmodel_[v[i]] = getParameterNameWithoutNamespace(v[i]);
  }

  mapParNamesFromPmodel_["MG94xREV.beta"] = "omega";
  mapParNamesFromPmodel_["MG94xREV.123_GTR.a"] = "a";
  mapParNamesFromPmodel_["MG94xREV.123_GTR.b"] = "b";
  mapParNamesFromPmodel_["MG94xREV.123_GTR.c"] = "c";
  mapParNamesFromPmodel_["MG94xREV.123_GTR.d"] = "d";
  mapParNamesFromPmodel_["MG94xREV.123_GTR.e"] = "e";

  updateMatrices_();
}

MG94xREV::MG94xREV(const MG94xREV& mg94) :
  AbstractParameterAliasable(mg94),
  AbstractWrappedModel(mg94),
  AbstractWrappedTransitionModel(mg94),
  AbstractTotallyWrappedTransitionModel(mg94),
  AbstractBiblioTransitionModel(mg94),
  AbstractWrappedSubstitutionModel(mg94),
  AbstractTotallyWrappedSubstitutionModel(mg94),
  AbstractBiblioSubstitutionModel(mg94),
  pmodel_(new CodonDistancePhaseFrequenciesSubstitutionModel(*mg94.pmodel_))
{}

MG94xREV& MG94xREV::operator=(const MG94xREV& mg94)
{
  AbstractBiblioSubstitutionModel::operator=(mg94);
  pmodel_.reset(new CodonDistancePhaseFrequenciesSubstitutionModel(*mg94.pmodel_));
  return *this;
}

MG94xREV::~MG94xREV() {}
