// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>

#include "../MixtureOfASubstitutionModel.h"
#include "DFP07.h"
#include "DFPDistanceFrequenciesSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

DFP07::DFP07(
    shared_ptr<const GeneticCode> gCode,
    unique_ptr<ProteinSubstitutionModelInterface> pAAmodel,
    unique_ptr<CodonFrequencySetInterface> codonFreqs) :
  AbstractParameterAliasable("DFP07."),
  AbstractWrappedModel("DFP07."),
  AbstractWrappedTransitionModel("DFP07."),
  AbstractTotallyWrappedTransitionModel("DFP07."),
  AbstractBiblioTransitionModel("DFP07."),
  AbstractBiblioMixedTransitionModel("DFP07."),
  mixedSubModelPtr_(nullptr),
  synfrom_(),
  synto_()
{
  // build the submodel

  auto codmodel = make_unique<DFPDistanceFrequenciesSubstitutionModel>(gCode, std::move(codonFreqs));
  auto submodel = make_unique<CodonSameAARateSubstitutionModel>(std::move(pAAmodel), std::move(codmodel),  unique_ptr<CodonFrequencySetInterface>(nullptr), gCode);

  map<string, unique_ptr<DiscreteDistributionInterface>> mpdd;
  vector<double> v1, v2;
  v1.push_back(0.5); v1.push_back(1);
  v2.push_back(0.5); v2.push_back(0.5);
  mpdd["DFPDistFreq.beta"] = make_unique<SimpleDiscreteDistribution>(v1, v2);

  mixedModelPtr_.reset(new MixtureOfASubstitutionModel(gCode->getSourceAlphabet(), std::move(submodel), mpdd));
  mixedSubModelPtr_ = dynamic_cast<const MixtureOfASubstitutionModel*>(&mixedModel());

  vector<int> supportedChars = mixedSubModelPtr_->getAlphabetStates();
  // map the parameters

  lParPmodel_.addParameters(mixedModelPtr_->getParameters());

  vector<std::string> v = mixedModelPtr_->getNModel(0)->getParameters().getParameterNames();

  for (auto vi : v)
  {
    if (vi != "SameAARate.DFPDistFreq.beta")
      mapParNamesFromPmodel_[vi] = vi.substr(23);
  }

  mapParNamesFromPmodel_["SameAARate.DFPDistFreq.beta_Simple.V1"] = "omega";
  mapParNamesFromPmodel_["SameAARate.DFPDistFreq.beta_Simple.theta1"] = "p0";

  // specific parameters

  string st;
  for (auto it : mapParNamesFromPmodel_)
  {
    st = mixedModelPtr_->getParameterNameWithoutNamespace(it.first);
    if (st != "DFPDistFreq.beta_Simple.V1")
    {
      addParameter_(new Parameter("DFP07." + it.second, mixedModelPtr_->getParameterValue(st),
                                  mixedModelPtr_->parameter(st).hasConstraint() ? std::shared_ptr<ConstraintInterface>(mixedModelPtr_->parameter(st).getConstraint()->clone()) : 0));
    }
  }

  addParameter_(new Parameter("DFP07.omega", 0.5, std::make_shared<IntervalConstraint>(0.002, 1, true, false, 0.002)));

  // look for synonymous codons
  for (synfrom_ = 1; synfrom_ < supportedChars.size(); ++synfrom_)
  {
    for (synto_ = 0; synto_ < synfrom_; ++synto_)
    {
      if (gCode->areSynonymous(supportedChars[synfrom_], supportedChars[synto_])
          && (mixedSubModelPtr_->subNModel(0).Qij(synfrom_, synto_) != 0)
          && (mixedSubModelPtr_->subNModel(1).Qij(synfrom_, synto_) != 0))
        break;
    }
    if (synto_ < synfrom_)
      break;
  }

  if (synto_ == supportedChars.size())
    throw Exception("Impossible to find synonymous codons");

  // update matrice

  computeFrequencies(false);
  updateMatrices_();
}

void DFP07::updateMatrices_()
{
  AbstractBiblioTransitionModel::updateMatrices_();

  // homogeneization of the synonymous substitution rates

  Vdouble vd;

  vd.push_back(1. / mixedSubModelPtr_->subNModel(0).Qij(synfrom_, synto_));
  vd.push_back(1. / mixedSubModelPtr_->subNModel(1).Qij(synfrom_, synto_));

  mixedModelPtr_->setVRates(vd);
}
