//
// File: LG10_EX_EHO.cpp
// Authors:
//   Mathieu Groussin
// Created: jeudi 21 octobre 2010, ÃÂ  14h 35
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

#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>

#include "../FrequencySet/ProteinFrequencySet.h"
#include "../MixtureOfSubstitutionModels.h"
#include "LG10_EX_EHO.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

LG10_EX_EHO::LG10_EX_EHO(shared_ptr<const ProteicAlphabet> alpha) :
  AbstractParameterAliasable("LG10_EX_EHO."),
  AbstractWrappedModel("LG10_EX_EHO."),
  AbstractWrappedTransitionModel("LG10_EX_EHO."),
  AbstractTotallyWrappedTransitionModel("LG10_EX_EHO."),
  AbstractBiblioTransitionModel("LG10_EX_EHO."),
  AbstractBiblioMixedTransitionModel("LG10_EX_EHO.")
{
  // build the submodel

  vector<unique_ptr<TransitionModelInterface>> vpSM;
  vpSM.push_back(make_unique<LG10_EX_EHO::EmbeddedModel>(alpha, "BUR_EXT"));
  vpSM.push_back(make_unique<LG10_EX_EHO::EmbeddedModel>(alpha, "BUR_HEL"));
  vpSM.push_back(make_unique<LG10_EX_EHO::EmbeddedModel>(alpha, "BUR_OTH"));
  vpSM.push_back(make_unique<LG10_EX_EHO::EmbeddedModel>(alpha, "EXP_EXT"));
  vpSM.push_back(make_unique<LG10_EX_EHO::EmbeddedModel>(alpha, "EXP_HEL"));
  vpSM.push_back(make_unique<LG10_EX_EHO::EmbeddedModel>(alpha, "EXP_OTH"));

  Vdouble vrate, vproba;

  for (auto& vi : vpSM)
  {
    vproba.push_back((dynamic_cast<LG10_EX_EHO::EmbeddedModel&>(*vi)).getProportion());
    vrate.push_back(vi->getRate());
  }

  mixedModelPtr_.reset(new MixtureOfSubstitutionModels(alpha, vpSM, vproba, vrate));

  string name, st;
  ParameterList pl = mixedModelPtr_->getParameters();
  for (unsigned int i = 0; i < pl.size(); ++i)
  {
    name = pl[i].getName();
    lParPmodel_.addParameter(Parameter(pl[i]));
    st = mixedModelPtr_->getParameterNameWithoutNamespace(name);
    mapParNamesFromPmodel_[name] = st;
    addParameter_(new Parameter("LG10_EX_EHO." + st,
                                mixedModelPtr_->getParameterValue(st),
                                mixedModelPtr_->getParameter(st).hasConstraint() ? std::shared_ptr<ConstraintInterface>(mixedModelPtr_->getParameter(st).getConstraint()->clone()) : 0));
  }

  updateMatrices_();
}

/**************** sub model classes ********************/

LG10_EX_EHO::EmbeddedModel::EmbeddedModel(
    shared_ptr<const ProteicAlphabet> alpha,
    string name) :
  AbstractParameterAliasable(name),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), name),
  proportion_(1),
  name_(name)
{
#include "__LG10_EX_EHOExchangeabilityCode"
#include "__LG10_EX_EHOFrequenciesCode"
#include "__LG10_EX_EHORatesProps"
  updateMatrices_();
}
