// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "HmmSequenceEvolution.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

HmmSequenceEvolution::HmmSequenceEvolution(
    shared_ptr<SubstitutionProcessCollection> processColl,
    vector<size_t>& nProc) :
  MultiProcessSequenceEvolution(processColl, nProc),
  hmmAlph_(),
  hmmTransMat_()
{
  hmmAlph_ = make_shared<HmmProcessAlphabet>(processColl_, nProc);

  hmmTransMat_ = make_shared<FullHmmTransitionMatrix>(hmmAlph_, "HMM.");

  // initialize parameters:
  addParameters_(hmmAlph_->getParameters());
  addParameters_(hmmTransMat_->getParameters());
}


void HmmSequenceEvolution::setNamespace(const std::string& nameSpace)
{
  deleteParameters_(hmmAlph_->getParameters().getParameterNames());
  deleteParameters_(hmmTransMat_->getParameters().getParameterNames());

  hmmAlph_->setNamespace(nameSpace);
  hmmTransMat_->setNamespace(nameSpace);

  addParameters_(hmmAlph_->getParameters());
  addParameters_(hmmTransMat_->getParameters());
}

void HmmSequenceEvolution::fireParameterChanged(const ParameterList& parameters)
{
  MultiProcessSequenceEvolution::fireParameterChanged(parameters);
  hmmAlph_->matchParametersValues(parameters);
  hmmTransMat_->matchParametersValues(parameters);
}
