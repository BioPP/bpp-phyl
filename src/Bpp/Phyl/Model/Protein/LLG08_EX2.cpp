//
// File: LLG08_EX2.cpp
// Authors:
//   Laurent Gueguen
// Created: mardi 12 octobre 2010, ÃÂ  09h 43
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
#include "LLG08_EX2.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

LLG08_EX2::LLG08_EX2(shared_ptr<const ProteicAlphabet> alpha) :
  AbstractParameterAliasable("LLG08_EX2."),
  AbstractWrappedModel("LLG08_EX2."),
  AbstractWrappedTransitionModel("LLG08_EX2."),
  AbstractTotallyWrappedTransitionModel("LLG08_EX2."),
  AbstractBiblioTransitionModel("LLG08_EX2."),
  AbstractBiblioMixedTransitionModel("LLG08_EX2.")
{
  // build the submodel

  vector<unique_ptr<TransitionModelInterface>> vpSM;
  vpSM.push_back(make_unique<LLG08_EX2::EmbeddedModel>(alpha, "Buried"));
  vpSM.push_back(make_unique<LLG08_EX2::EmbeddedModel>(alpha, "Exposed"));

  Vdouble vrate, vproba;

  for (auto& vi : vpSM)
  {
    vproba.push_back(dynamic_cast<LLG08_EX2::EmbeddedModel&>(*vi).getProportion());
    vrate.push_back(vi->getRate());
  }

  mixedModelPtr_.reset(new MixtureOfSubstitutionModels(alpha, vpSM, vproba, vrate));

  string name, st;
  ParameterList pl = mixedModelPtr_->getParameters();
  for (size_t i = 0; i < pl.size(); ++i)
  {
    name = pl[i].getName();
    lParPmodel_.addParameter(Parameter(pl[i]));
    st = mixedModelPtr_->getParameterNameWithoutNamespace(name);
    mapParNamesFromPmodel_[name] = st;
    addParameter_(new Parameter("LLG08_EX2." + st,
                                mixedModelPtr_->getParameterValue(st),
                                mixedModelPtr_->getParameter(st).hasConstraint() ? std::shared_ptr<ConstraintInterface>(mixedModelPtr_->getParameter(st).getConstraint()->clone()) : 0));
  }

  updateMatrices_();
}

/**************** sub model classes *******************/

LLG08_EX2::EmbeddedModel::EmbeddedModel(
    shared_ptr<const ProteicAlphabet> alpha,
    string name) :
  AbstractParameterAliasable(name),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), name),
  proportion_(1),
  name_(name)
{
#include "__LLG08_EX2ExchangeabilityCode"
#include "__LLG08_EX2FrequenciesCode"
#include "__LLG08_EX2RatesProps"
  updateMatrices_();
}
