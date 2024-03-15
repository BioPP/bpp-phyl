// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/NumConstants.h>

#include "AbstractCodonAARateSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractCodonAARateSubstitutionModel::AbstractCodonAARateSubstitutionModel(
  shared_ptr<ProteinSubstitutionModelInterface> pmodel,
  shared_ptr<const GeneticCode> pgencode,
  const string& prefix,
  bool paramSynRate) :
  AbstractParameterAliasable(prefix),
  pAAmodel_(pmodel),
  pgencode_(pgencode),
  beta_(19),
  gamma_(1),
  stateMap_(new CanonicalStateMap(pgencode->getSourceAlphabet(), false))
{
  if (paramSynRate)
    addParameter_(new Parameter(prefix + "gamma", 1, std::make_shared<IntervalConstraint>(NumConstants::SMALL(), 999, true, true)));

  addParameter_(new Parameter(prefix + "beta", 1, Parameter::R_PLUS_STAR));

  pAAmodel_->enableEigenDecomposition(false);

  pAAmodel_->setNamespace(prefix + pAAmodel_->getNamespace());
  addParameters_(pAAmodel_->getParameters());
}

void AbstractCodonAARateSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  pAAmodel_->matchParametersValues(parameters);

  if (hasParameter("gamma"))
    gamma_ = getParameterValue("gamma");

  beta_ = getParameterValue("beta");
}

double AbstractCodonAARateSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  int si(stateMap_->getAlphabetStateAsInt(i)), sj(stateMap_->getAlphabetStateAsInt(j));

  return pgencode_->areSynonymous(si, sj) ? gamma_ :
         beta_* pAAmodel_->Qij(pAAmodel_->getModelStates(pgencode_->translate(si))[0],
                               pAAmodel_->getModelStates(pgencode_->translate(sj))[0]);
}
