//
// File: YNGP_M2.cpp
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
#include "YN98.h"
#include "YNGP_M2.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

YNGP_M2::YNGP_M2(
    shared_ptr<const GeneticCode> gc,
    unique_ptr<CodonFrequencySetInterface> codonFreqs) :
  AbstractParameterAliasable("YNGP_M2."),
  AbstractWrappedModel("YNGP_M2."),
  AbstractWrappedTransitionModel("YNGP_M2."),
  AbstractTotallyWrappedTransitionModel("YNGP_M2."),
  AbstractBiblioTransitionModel("YNGP_M2."),
  YNGP_M("YNGP_M2.")
{
  // build the submodel

  vector<double> v1, v2;
  v1.push_back(0.5);
  v1.push_back(1);
  v1.push_back(2);
  v2.push_back(0.333333);
  v2.push_back(0.333333);
  v2.push_back(0.333334);

  auto psdd = make_unique<SimpleDiscreteDistribution>(v1, v2);

  map<string, unique_ptr<DiscreteDistribution>> mpdd;
  mpdd["omega"] = move(psdd);

  auto yn98 = make_unique<YN98>(gc, move(codonFreqs));

  mixedModelPtr_.reset(new MixtureOfASubstitutionModel(gc->getSourceAlphabet(), move(yn98), mpdd));
  mixedSubModelPtr_ = dynamic_cast<const MixtureOfASubstitutionModel*>(&mixedModel());

  vector<int> supportedChars = yn98->getAlphabetStates();

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
  mapParNamesFromPmodel_["YN98.omega_Simple.V1"] = "omega0";
  mapParNamesFromPmodel_["YN98.omega_Simple.theta1"] = "theta1";
  mapParNamesFromPmodel_["YN98.omega_Simple.V3"] = "omega2";
  mapParNamesFromPmodel_["YN98.omega_Simple.theta2"] = "theta2";

  // specific parameters

  string st;
  for (auto it : mapParNamesFromPmodel_)
  {
    st = mixedModelPtr_->getParameterNameWithoutNamespace(it.first);
    if (it.second.substr(0, 5) != "omega")
      addParameter_(new Parameter("YNGP_M2." + it.second, mixedModelPtr_->getParameterValue(st),
                                  mixedModelPtr_->getParameter(st).hasConstraint() ? std::shared_ptr<Constraint>(mixedModelPtr_->getParameter(st).getConstraint()->clone()) : 0));
  }

  addParameter_(new Parameter("YNGP_M2.omega0", 0.5, std::make_shared<IntervalConstraint>(0.002, 1, true, false)));

  addParameter_(new Parameter("YNGP_M2.omega2", 2, std::make_shared<IntervalConstraint>(1, 999, false, false, 0.002)));

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

  // update Matrices

  computeFrequencies(false);
  updateMatrices_();
}

void YNGP_M2::updateMatrices_()
{
  AbstractBiblioTransitionModel::updateMatrices_();

  // homogeneization of the synonymous substittion rates

  Vdouble vd;

  vd.push_back(1 / mixedSubModelPtr_->subNModel(0).Qij(synfrom_, synto_));
  vd.push_back(1 / mixedSubModelPtr_->subNModel(1).Qij(synfrom_, synto_));
  vd.push_back(1 / mixedSubModelPtr_->subNModel(2).Qij(synfrom_, synto_));

  mixedModelPtr_->setVRates(vd);
}
