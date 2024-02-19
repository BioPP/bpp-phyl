//
// File: LGL08_CAT.cpp
// Authors:
//   Mathieu Groussin
// Created: 2012-12-11 00:00:00
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
#include "LGL08_CAT.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

LGL08_CAT::LGL08_CAT(
    shared_ptr<const ProteicAlphabet> alpha,
    unsigned int nbCat) :
  AbstractParameterAliasable("LGL08_CAT."),
  AbstractWrappedModel("LGL08_CAT."),
  AbstractWrappedTransitionModel("LGL08_CAT."),
  AbstractTotallyWrappedTransitionModel("LGL08_CAT."),
  AbstractBiblioTransitionModel("LGL08_CAT."),
  AbstractBiblioMixedTransitionModel("LGL08_CAT.")
{
  // build the submodel

  vector<unique_ptr<TransitionModelInterface>> vpSM;
  for (unsigned int i = 1; i < nbCat + 1; ++i)
  {
    vpSM.push_back(make_unique<LGL08_CAT::EmbeddedModel>(alpha, "C" + TextTools::toString(i), nbCat));
  }

  Vdouble vrate, vproba;

  for (auto& vi : vpSM)
  {
    vproba.push_back((dynamic_cast<LGL08_CAT::EmbeddedModel&>(*vi)).getProportion());
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
    addParameter_(new Parameter("LGL08_CAT." + st,
                                mixedModelPtr_->getParameterValue(st),
                                mixedModelPtr_->parameter(st).hasConstraint() ? shared_ptr<ConstraintInterface>(mixedModelPtr_->parameter(st).getConstraint()->clone()) : 0));
  }

  updateMatrices_();
}

/**************** sub model classes */ // ////////

LGL08_CAT::EmbeddedModel::EmbeddedModel(
    shared_ptr<const ProteicAlphabet> alpha,
    string name, 
    unsigned int nbCat) :
  AbstractParameterAliasable(name),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), name),
  proportion_(1),
  name_(name)
{
  // Exchangeabilities:
  for (unsigned int i = 0; i < 20; ++i)
  {
    for (unsigned int j = 0; j < 20; ++j)
    {
      if (i == j)
        exchangeability_(i, i) = -19.;
      else
        exchangeability_(i, j) = 1.;
    }
  }

  // Equilibrium frequencies, rates and proportions:
  if (nbCat == 10)
  {
#include "__CATC10FrequenciesCode"
#include "__CATC10RatesProps"
  }
  else if (nbCat == 20)
  {
#include "__CATC20FrequenciesCode"
#include "__CATC20RatesProps"
  }
  else if (nbCat == 30)
  {
#include "__CATC30FrequenciesCode"
#include "__CATC30RatesProps"
  }
  else if (nbCat == 40)
  {
#include "__CATC40FrequenciesCode"
#include "__CATC40RatesProps"
  }
  else if (nbCat == 50)
  {
#include "__CATC50FrequenciesCode"
#include "__CATC50RatesProps"
  }
  else if (nbCat == 60)
  {
#include "__CATC60FrequenciesCode"
#include "__CATC60RatesProps"
  }
  else
    throw Exception("LGL08_CAT.cpp: incorrect number of profiles. This number has to be 10, 20, 30, 40, 50 or 60.");

  updateMatrices_();
}

