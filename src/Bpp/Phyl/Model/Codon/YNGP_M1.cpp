// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>

#include "../MixtureOfASubstitutionModel.h"
#include "YN98.h"
#include "YNGP_M1.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

YNGP_M1::YNGP_M1(
    shared_ptr<const GeneticCode> gc,
    unique_ptr<CodonFrequencySetInterface> codonFreqs) :
  AbstractParameterAliasable("YNGP_M1."),
  AbstractWrappedModel("YNGP_M1."),
  AbstractWrappedTransitionModel("YNGP_M1."),
  AbstractTotallyWrappedTransitionModel("YNGP_M1."),
  AbstractBiblioTransitionModel("YNGP_M1."),
  YNGP_M("YNGP_M1.")
{
  // build the submodel

  vector<double> v1, v2;
  v1.push_back(0.5); v1.push_back(1);
  v2.push_back(0.5); v2.push_back(0.5);

  auto psdd = make_unique<SimpleDiscreteDistribution>(v1, v2);

  map<string, unique_ptr<DiscreteDistributionInterface>> mpdd;
  mpdd["omega"] = move(psdd);

  auto yn98 = make_unique<YN98>(gc, move(codonFreqs));

  mixedModelPtr_.reset(new MixtureOfASubstitutionModel(gc->getCodonAlphabet(), move(yn98), mpdd));
  mixedSubModelPtr_ = dynamic_cast<const MixtureOfASubstitutionModel*>(&mixedModel());

  vector<int> supportedChars = mixedModelPtr_->getAlphabetStates();
  // map the parameters

  lParPmodel_.addParameters(mixedModelPtr_->getParameters());

  vector<std::string> v = dynamic_cast<const YN98&>(mixedModelPtr_->nModel(0)).frequencySet().getParameters().getParameterNames();

  for (auto& vi : v)
  {
    mapParNamesFromPmodel_[vi] = vi.substr(5);
  }

  mapParNamesFromPmodel_["YN98.kappa"] = "kappa";
  mapParNamesFromPmodel_["YN98.omega_Simple.V1"] = "omega";
  mapParNamesFromPmodel_["YN98.omega_Simple.theta1"] = "p0";

  // specific parameters

  string st;
  for (auto it : mapParNamesFromPmodel_)
  {
    st = mixedModelPtr_->getParameterNameWithoutNamespace(it.first);
    if (st != "omega_Simple.V1")
    {
      addParameter_(new Parameter("YNGP_M1." + it.second, mixedModelPtr_->getParameterValue(st),
            mixedModelPtr_->parameter(st).hasConstraint() ? std::shared_ptr<ConstraintInterface>(mixedModelPtr_->parameter(st).getConstraint()->clone()) : 0));
    }
  }

  addParameter_(new Parameter("YNGP_M1.omega", 0.5, std::make_shared<IntervalConstraint>(0.002, 1, true, false, 0.002)));

  // look for synonymous codons
  for (synfrom_ = 1; synfrom_ < supportedChars.size(); ++synfrom_)
  {
    for (synto_ = 0; synto_ < synfrom_; ++synto_)
    {
      if (gc->areSynonymous(supportedChars[synfrom_], supportedChars[synto_])
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

void YNGP_M1::updateMatrices_()
{
  AbstractBiblioTransitionModel::updateMatrices_();

  // homogeneization of the synonymous substitution rates

  Vdouble vd;

  vd.push_back(1. / mixedSubModelPtr_->subNModel(0).Qij(synfrom_, synto_));
  vd.push_back(1. / mixedSubModelPtr_->subNModel(1).Qij(synfrom_, synto_));

  mixedModelPtr_->setVRates(vd);
}
