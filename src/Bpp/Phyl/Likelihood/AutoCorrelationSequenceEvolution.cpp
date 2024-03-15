// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AutoCorrelationSequenceEvolution.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

AutoCorrelationSequenceEvolution::AutoCorrelationSequenceEvolution(
  std::shared_ptr<SubstitutionProcessCollection> processColl,
  std::vector<size_t>& nProc) :
  MultiProcessSequenceEvolution(processColl, nProc),
  hmmAlph_(),
  autoCorrTransMat_()
{
  hmmAlph_ = make_unique<HmmProcessAlphabet>(processColl_, nProc);

  autoCorrTransMat_ = make_unique<AutoCorrelationTransitionMatrix>(hmmAlph_, "AutoCorr.");

  // initialize parameters:
  addParameters_(hmmAlph_->getParameters());
  addParameters_(autoCorrTransMat_->getParameters());
}


void AutoCorrelationSequenceEvolution::setNamespace(const std::string& nameSpace)
{
  deleteParameters_(hmmAlph_->getParameters().getParameterNames());
  deleteParameters_(autoCorrTransMat_->getParameters().getParameterNames());

  hmmAlph_->setNamespace(nameSpace);
  autoCorrTransMat_->setNamespace(nameSpace);

  addParameters_(hmmAlph_->getParameters());
  addParameters_(autoCorrTransMat_->getParameters());
}


void AutoCorrelationSequenceEvolution::fireParameterChanged(const ParameterList& parameters)
{
  MultiProcessSequenceEvolution::fireParameterChanged(parameters);

  hmmAlph_->matchParametersValues(parameters);
  autoCorrTransMat_->matchParametersValues(parameters);
}

ParameterList AutoCorrelationSequenceEvolution::getNonDerivableParameters() const
{
  ParameterList pl = MultiProcessSequenceEvolution::getNonDerivableParameters();

  pl.includeParameters(autoCorrTransMat_->getParameters());

  return pl;
}
