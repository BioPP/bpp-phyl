//
// File: DFP07.cpp
// Authors:
//   Laurent Gueguen
// Created: 2010-05-08 00:00:00
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

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
