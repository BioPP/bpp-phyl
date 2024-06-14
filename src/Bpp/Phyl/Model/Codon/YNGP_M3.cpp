// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>

#include "../MixtureOfASubstitutionModel.h"
#include "YN98.h"
#include "YNGP_M3.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

YNGP_M3::YNGP_M3(
    shared_ptr<const GeneticCode> gc,
    unique_ptr<CodonFrequencySetInterface> codonFreqs,
    unsigned int nbOmega) :
  AbstractParameterAliasable("YNGP_M3."),
  AbstractWrappedModel("YNGP_M3."),
  AbstractWrappedTransitionModel("YNGP_M3."),
  AbstractTotallyWrappedTransitionModel("YNGP_M3."),
  AbstractBiblioTransitionModel("YNGP_M3."),
  YNGP_M("YNGP_M3.")
{
  if (nbOmega < 1)
    throw Exception("At least one omega is necessary in the YNGP_M3 model");

  // build the submodel

  vector<double> v1, v2;
  v1.push_back(0.5);
  for (unsigned int i = 1; i < nbOmega; ++i)
  {
    v1.push_back(0.5 + 0.5 * i);
  }

  for (unsigned int i = 0; i < nbOmega; ++i)
  {
    v2.push_back(1. / nbOmega);
  }

  auto psdd = make_unique<SimpleDiscreteDistribution>(v1, v2, 0.002);

  map<string, unique_ptr<DiscreteDistributionInterface>> mpdd;
  mpdd["omega"] = std::move(psdd);

  auto yn98 = make_unique<YN98>(gc, std::move(codonFreqs));

  mixedModelPtr_.reset(new MixtureOfASubstitutionModel(gc->getSourceAlphabet(), std::move(yn98), mpdd));
  mixedSubModelPtr_ = dynamic_cast<const MixtureOfASubstitutionModel*>(&mixedModel());

  vector<int> supportedChars = mixedSubModelPtr_->getAlphabetStates();

  // mapping the parameters

  ParameterList pl = mixedModelPtr_->getParameters();
  for (size_t i = 0; i < pl.size(); ++i)
  {
    lParPmodel_.addParameter(Parameter(pl[i]));
  }

  vector<std::string> v = dynamic_cast<const YN98&>(mixedModelPtr_->nModel(0)).frequencySet().getParameters().getParameterNames();
  for (auto& vi : v)
  {
    mapParNamesFromPmodel_[vi] = vi.substr(5);
  }

  mapParNamesFromPmodel_["YN98.kappa"] = "kappa";

  for (size_t i = 1; i < nbOmega; ++i)
  {
    mapParNamesFromPmodel_["YN98.omega_Simple.theta" + TextTools::toString(i)] = "theta" + TextTools::toString(i);
  }


  mapParNamesFromPmodel_["YN98.omega_Simple.V1"] = "omega0";
  for (size_t i = 1; i < nbOmega; ++i)
  {
    mapParNamesFromPmodel_["YN98.omega_Simple.V" + TextTools::toString(i + 1)] = "delta" + TextTools::toString(i);
  }

  // specific parameters

  string st;
  for (auto& it : mapParNamesFromPmodel_)
  {
    st = mixedModelPtr_->getParameterNameWithoutNamespace(it.first);
    if (it.second.substr(0, 5) != "delta")
      addParameter_(new Parameter("YNGP_M3." + it.second, mixedModelPtr_->getParameterValue(st),
            mixedModelPtr_->parameter(st).hasConstraint() ? shared_ptr<ConstraintInterface>(mixedModelPtr_->parameter(st).getConstraint()->clone()) : 0));
  }

  for (size_t i = 1; i < nbOmega; ++i)
  {
    addParameter_(new Parameter("YNGP_M3.delta" + TextTools::toString(i), 0.5, make_shared<IntervalConstraint>(0.002, 999, true, true, 0.002)));
  }

  // look for synonymous codons
  for (synfrom_ = 1; synfrom_ < supportedChars.size(); synfrom_++)
  {
    for (synto_ = 0; synto_ < synfrom_; synto_++)
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

  // update Matrices
  computeFrequencies(false);
  updateMatrices_();
}

void YNGP_M3::updateMatrices_()
{
  for (unsigned int i = 0; i < lParPmodel_.size(); ++i)
  {
    if (mapParNamesFromPmodel_.find(lParPmodel_[i].getName()) != mapParNamesFromPmodel_.end())
    {
      const string& np = lParPmodel_[i].getName();
      if (np.size() > 19 && np[18] == 'V')
      {
        size_t ind = TextTools::to<size_t>(np.substr(19));
        double x = getParameterValue("omega0");
        for (unsigned j = 1; j < ind; j++)
        {
          x += getParameterValue("delta" + TextTools::toString(j));
        }

        const auto parConst = lParPmodel_[i].getConstraint();
        lParPmodel_[i].setValue(parConst ? (parConst->isCorrect(x) ? x : parConst->getAcceptedLimit(x)) : x);
      }
      else
        lParPmodel_[i].setValue(parameter(getParameterNameWithoutNamespace(mapParNamesFromPmodel_[np])).getValue());
    }
  }

  mixedModelPtr_->matchParametersValues(lParPmodel_);

  // homogeneization of the synonymous substitution rates

  Vdouble vd;

  for (unsigned int i = 0; i < mixedModelPtr_->getNumberOfModels(); ++i)
  {
    vd.push_back(1 / mixedSubModelPtr_->subNModel(i).Qij(synfrom_, synto_));
  }

  mixedModelPtr_->setVRates(vd);
}
