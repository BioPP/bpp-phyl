//
// File: YNGP_M8.cpp
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
    shared_ptr<const GeneticCode> gc,
    unique_ptr<CodonFrequencySetInterface> codonFreqs,
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
  auto psdd = make_unique<SimpleDiscreteDistribution>(val, move(prob));

  vector<unique_ptr<DiscreteDistribution>> v_distr;
  v_distr.push_back(move(pbdd)); 
  v_distr.push_back(move(psdd));
  prob.clear(); 
  prob.push_back(0.5);
  prob.push_back(0.5);

  // Distribution on omega = mixture (Beta + Simple of size 1)
  auto pmodd = make_unique<MixtureOfDiscreteDistributions>(v_distr, prob);

  map<string, unique_ptr<DiscreteDistribution>> mpdd;
  mpdd["omega"] = move(pmodd);

  auto yn98 = make_unique<YN98>(gc, move(codonFreqs));

  mixedModelPtr_.reset(new MixtureOfASubstitutionModel(gc->getSourceAlphabet(), move(yn98), mpdd));
  mixedSubModelPtr_ = dynamic_cast<const MixtureOfASubstitutionModel*>(&mixedModel());

  vector<int> supportedChars = mixedSubModelPtr_->getAlphabetStates();

  // mapping the parameters

  ParameterList pl = mixedModelPtr_->getParameters();
  for (size_t i = 0; i < pl.size(); ++i)
  {
    if (neutral_ && pl[i].getName()=="YN98.omega_Mixture.2_Simple.V1")
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
      addParameter_(new Parameter(getName()+"." + it.second, mixedModelPtr_->getParameterValue(st),
                                  mixedModelPtr_->getParameter(st).hasConstraint() ? std::shared_ptr<ConstraintInterface>(mixedModelPtr_->getParameter(st).getConstraint()->clone()) : 0));
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

