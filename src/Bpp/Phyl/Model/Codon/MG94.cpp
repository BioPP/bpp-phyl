// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../FrequencySet/CodonFrequencySet.h"
#include "../Nucleotide/K80.h"
#include "MG94.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

MG94::MG94(
    shared_ptr<const GeneticCode> gc,
    unique_ptr<CodonFrequencySetInterface> codonFreqs) :
  AbstractParameterAliasable("MG94."),
  AbstractWrappedModel("MG94."),
  AbstractWrappedTransitionModel("MG94."),
  AbstractTotallyWrappedTransitionModel("MG94."),
  AbstractBiblioTransitionModel("MG94."),
  AbstractWrappedSubstitutionModel("MG94."),
  AbstractTotallyWrappedSubstitutionModel("MG94."),
  AbstractBiblioSubstitutionModel("MG94."),
  pmodel_(new CodonDistancePhaseFrequenciesSubstitutionModel(
        gc,
        make_unique<K80>(gc->codonAlphabet().getNucleicAlphabet()),
        std::move(codonFreqs)))
{
  addParameter_(new Parameter("MG94.rho", 1, std::make_shared<IntervalConstraint>(0.002, 999, true, true)));

  pmodel_->setNamespace("MG94.");
  addParameters_(pmodel_->frequencySet().getParameters());

  lParPmodel_.addParameters(pmodel_->getParameters());

  vector<std::string> v = pmodel_->frequencySet().getParameters().getParameterNames();
  for (unsigned int i = 0; i < v.size(); ++i)
  {
    mapParNamesFromPmodel_[v[i]] = getParameterNameWithoutNamespace(v[i]);
  }

  mapParNamesFromPmodel_["MG94.beta"] = "rho";

  updateMatrices_();
}

MG94::MG94(const MG94& mg94) :
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

MG94& MG94::operator=(const MG94& mg94)
{
  AbstractBiblioSubstitutionModel::operator=(mg94);
  pmodel_.reset(new CodonDistancePhaseFrequenciesSubstitutionModel(*mg94.pmodel_));
  return *this;
}

MG94::~MG94() {}
