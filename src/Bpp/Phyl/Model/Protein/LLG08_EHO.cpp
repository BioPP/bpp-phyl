//
// File: LLG08_EHO.cpp
// Created by:  Laurent Gueguen
// Created on: jeudi 7 octobre 2010, à 21h 47
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)
   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#include "LLG08_EHO.h"
#include "../FrequenciesSet/ProteinFrequenciesSet.h"

#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>

using namespace bpp;

using namespace std;

/******************************************************************************/

LLG08_EHO::LLG08_EHO(const ProteicAlphabet* alpha) :
  AbstractBiblioMixedSubstitutionModel("LLG08_EHO."),
  pmixmodel_(0)
{
  // build the submodel

  vector<SubstitutionModel*> vpSM;
  vpSM.push_back(new LLG08_EHO::EmbeddedModel(alpha, "Extended"));
  vpSM.push_back(new LLG08_EHO::EmbeddedModel(alpha, "Helix"));
  vpSM.push_back(new LLG08_EHO::EmbeddedModel(alpha, "Other"));

  Vdouble vrate, vproba;

  for (unsigned int i = 0; i < 3; i++)
  {
    vproba.push_back((dynamic_cast<LLG08_EHO::EmbeddedModel*>(vpSM[i]))->getProportion());
    vrate.push_back(vpSM[i]->getRate());
  }

  pmixmodel_.reset(new MixtureOfSubstitutionModels(alpha, vpSM, vproba, vrate));

  string name, st;
  ParameterList pl = pmixmodel_->getParameters();
  for (size_t i = 0; i < pl.size(); i++)
  {
    name = pl[i].getName();
    lParPmodel_.addParameter(Parameter(pl[i]));
    st = pmixmodel_->getParameterNameWithoutNamespace(name);
    mapParNamesFromPmodel_[name] = st;
    addParameter_(new Parameter("LLG08_EHO." + st,
                                pmixmodel_->getParameterValue(st),
                                pmixmodel_->getParameter(st).hasConstraint() ? pmixmodel_->getParameter(st).getConstraint()->clone() : 0, true));
  }

  updateMatrices();
}

LLG08_EHO::LLG08_EHO(const LLG08_EHO& mod2) : AbstractBiblioMixedSubstitutionModel(mod2),
  pmixmodel_(new MixtureOfSubstitutionModels(*mod2.pmixmodel_))
{}

LLG08_EHO& LLG08_EHO::operator=(const LLG08_EHO& mod2)
{
  AbstractBiblioMixedSubstitutionModel::operator=(mod2);

  pmixmodel_.reset(new MixtureOfSubstitutionModels(*mod2.pmixmodel_));
  return *this;
}

LLG08_EHO::~LLG08_EHO() {}

/**************** sub model classes */ // ////////

LLG08_EHO::EmbeddedModel::EmbeddedModel(const ProteicAlphabet* alpha, string name) :
  AbstractParameterAliasable(name),
  AbstractReversibleProteinSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), name),
  proportion_(1),
  name_(name)
{
#include "__LLG08_EHOExchangeabilityCode"
#include "__LLG08_EHOFrequenciesCode"
#include "__LLG08_EHORatesProps"
  updateMatrices();
}


