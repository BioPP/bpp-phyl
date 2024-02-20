//
// File: YNGP_M9.cpp
// Authors:
//   Laurent Gueguen
// Created: lundi 15 dÃÂ©cembre 2014, ÃÂ  15h 21
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
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/MixtureOfDiscreteDistributions.h>
#include <Bpp/Text/TextTools.h>

#include "../MixtureOfASubstitutionModel.h"
#include "YN98.h"
#include "YNGP_M9.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

YNGP_M9::YNGP_M9(
    shared_ptr<const GeneticCode> gc,
    unique_ptr<CodonFrequencySetInterface> codonFreqs,
    unsigned int nbBeta,
    unsigned int nbGamma) :
  AbstractParameterAliasable("YNGP_M9."),
  AbstractWrappedModel("YNGP_M9."),
  AbstractWrappedTransitionModel("YNGP_M9."),
  AbstractTotallyWrappedTransitionModel("YNGP_M9."),
  AbstractBiblioTransitionModel("YNGP_M9."),
  YNGP_M("YNGP_M9."),
  nBeta_(nbBeta),
  nGamma_(nbGamma)
{
  if (nbBeta <= 0)
    throw Exception("Bad number of classes for beta distribution of model YNGP_M9: " + TextTools::toString(nbBeta));
  if (nbGamma <= 0)
    throw Exception("Bad number of classes for gamma distribution of model YNGP_M9: " + TextTools::toString(nbGamma));

  // build the submodel

  auto pbdd = make_unique<BetaDiscreteDistribution>(nbBeta, 2, 2);
  auto pgdd = make_unique<GammaDiscreteDistribution>(nbGamma, 1, 1);

  vector<unique_ptr<DiscreteDistributionInterface>> v_distr;
  v_distr.push_back(move(pbdd));
  v_distr.push_back(move(pgdd));
  vector<double> prob;
  prob.push_back(0.5);
  prob.push_back(0.5);

  auto pmodd = make_unique<MixtureOfDiscreteDistributions>(v_distr, prob);

  map<string, unique_ptr<DiscreteDistributionInterface>> mpdd;
  mpdd["omega"] = move(pmodd);

  auto yn98 = make_unique<YN98>(gc, move(codonFreqs));
  mixedModelPtr_.reset(new MixtureOfASubstitutionModel(gc->getSourceAlphabet(), move(yn98), mpdd));
  mixedSubModelPtr_ = dynamic_cast<const MixtureOfASubstitutionModel*>(&mixedModel());

  vector<int> supportedChars = mixedModelPtr_->getAlphabetStates();

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
  mapParNamesFromPmodel_["YN98.omega_Mixture.theta1"] = "p0";
  mapParNamesFromPmodel_["YN98.omega_Mixture.1_Beta.alpha"] = "p";
  mapParNamesFromPmodel_["YN98.omega_Mixture.1_Beta.beta"] = "q";
  mapParNamesFromPmodel_["YN98.omega_Mixture.2_Gamma.alpha"] = "alpha";
  mapParNamesFromPmodel_["YN98.omega_Mixture.2_Gamma.beta"] = "beta";

  // specific parameters

  string st;
  for (auto& it : mapParNamesFromPmodel_)
  {
    st = mixedModelPtr_->getParameterNameWithoutNamespace(it.first);
    addParameter_(new Parameter("YNGP_M9." + it.second, mixedModelPtr_->getParameterValue(st),
                                mixedModelPtr_->parameter(st).hasConstraint() ? shared_ptr<ConstraintInterface>(mixedModelPtr_->parameter(st).getConstraint()->clone()) : 0));
  }

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

void YNGP_M9::updateMatrices_()
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

