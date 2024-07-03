// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Prob/BetaDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/MixtureOfDiscreteDistributions.h>
#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>
#include <Bpp/Text/TextTools.h>

#include "../MixtureOfASubstitutionModel.h"
#include "YN98.h"
#include "YNGP_M8.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

YNGP_M8::YNGP_M8(
    std::shared_ptr<const GeneticCode> gc,
    std::unique_ptr<CodonFrequencySetInterface> codonFreqs,
    unsigned int nclass,
    bool neutral) :
  AbstractParameterAliasable(neutral ? "YNGP_M8a." : "YNGP_M8."),
  AbstractWrappedModel(neutral ? "YNGP_M8a." : "YNGP_M8."),
  AbstractWrappedTransitionModel(neutral ? "YNGP_M8a." : "YNGP_M8."),
  AbstractTotallyWrappedTransitionModel(neutral ? "YNGP_M8a." : "YNGP_M8."),
  AbstractBiblioTransitionModel(neutral ? "YNGP_M8a." : "YNGP_M8."),
  YNGP_M(neutral ? "YNGP_M8a." : "YNGP_M8."),
  neutral_(neutral)
{
  if (nclass <= 0)
    throw Exception("Bad number of classes for model " + getName() + ": " + TextTools::toString(nclass));

  // build the submodel

  auto pbdd = make_unique<BetaDiscreteDistribution>(nclass, 2, 2);

  vector<double> val = { neutral_ ? 1. : 2. };
  vector<double> prob = { 1. };
  auto psdd = make_unique<SimpleDiscreteDistribution>(val, std::move(prob));

  vector<unique_ptr<DiscreteDistributionInterface>> v_distr;
  v_distr.push_back(std::move(pbdd));
  v_distr.push_back(std::move(psdd));
  prob.clear();
  prob.push_back(0.5);
  prob.push_back(0.5);

  // Distribution on omega = mixture (Beta + Simple of size 1)
  auto pmodd = make_unique<MixtureOfDiscreteDistributions>(v_distr, prob);

  map<string, unique_ptr<DiscreteDistributionInterface>> mpdd;
  mpdd["omega"] = std::move(pmodd);

  auto yn98 = make_unique<YN98>(gc, std::move(codonFreqs));

  mixedModelPtr_.reset(new MixtureOfASubstitutionModel(gc->getSourceAlphabet(), std::move(yn98), mpdd));
  mixedSubModelPtr_ = dynamic_cast<const MixtureOfASubstitutionModel*>(&mixedModel());

  vector<int> supportedChars = mixedSubModelPtr_->getAlphabetStates();

  // mapping the parameters

  ParameterList pl = mixedModelPtr_->getParameters();
  for (size_t i = 0; i < pl.size(); ++i)
  {
    if (neutral_ && pl[i].getName() == "YN98.omega_Mixture.2_Simple.V1")
      continue;
    lParPmodel_.addParameter(Parameter(pl[i]));
  }

  vector<std::string> v = dynamic_cast<const YN98&>(mixedModelPtr_->nModel(0)).frequencySet().getParameters().getParameterNames();

  for (auto& vi : v)
  {
    mapParNamesFromPmodel_[vi] = vi.substr(5);
  }

  mapParNamesFromPmodel_["YN98.kappa"] = "kappa";
  mapParNamesFromPmodel_["YN98.omega_Mixture.theta1"] = "p0";
  mapParNamesFromPmodel_["YN98.omega_Mixture.1_Beta.alpha"] = "p";
  mapParNamesFromPmodel_["YN98.omega_Mixture.1_Beta.beta"] = "q";
  if (!neutral_)
    mapParNamesFromPmodel_["YN98.omega_Mixture.2_Simple.V1"] = "omegas";

  // specific parameters

  string st;
  for (auto& it : mapParNamesFromPmodel_)
  {
    st = mixedModelPtr_->getParameterNameWithoutNamespace(it.first);
    if (it.second != "omegas")
      addParameter_(new Parameter(getName() + "." + it.second, mixedModelPtr_->getParameterValue(st),
            mixedModelPtr_->parameter(st).hasConstraint() ? std::shared_ptr<ConstraintInterface>(mixedModelPtr_->parameter(st).getConstraint()->clone()) : 0));
  }

  if (!neutral_)
    addParameter_(new Parameter("YNGP_M8.omegas", 2., make_shared<IntervalConstraint>(1, 1, false)));

  // look for synonymous codons
  for (synfrom_ = 1; synfrom_ < supportedChars.size(); synfrom_++)
  {
    for (synto_ = 0; synto_ < synfrom_; synto_++)
    {
      if ((gc->areSynonymous(supportedChars[synfrom_], supportedChars[synto_]))
          && (mixedSubModelPtr_->subNModel(0).Qij(synfrom_, synto_) != 0)
          && (mixedSubModelPtr_->subNModel(1).Qij(synfrom_, synto_) != 0))
        break;
    }
    if (synto_ < synfrom_)
      break;
  }

  if (synto_ == gc->getSourceAlphabet()->getSize())
    throw Exception("Impossible to find synonymous codons");

  // update Matrices
  computeFrequencies(false);
  updateMatrices_();
}

void YNGP_M8::updateMatrices_()
{
  AbstractBiblioTransitionModel::updateMatrices_();

  // homogeneization of the synonymous substittion rates

  Vdouble vd;

  for (unsigned int i = 0; i < mixedModelPtr_->getNumberOfModels(); ++i)
  {
    vd.push_back(1 / mixedSubModelPtr_->subNModel(i).Qij(synfrom_, synto_));
  }

  mixedModelPtr_->setVRates(vd);
}
