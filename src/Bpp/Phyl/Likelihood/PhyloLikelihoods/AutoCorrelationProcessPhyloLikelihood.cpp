// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Hmm/AutoCorrelationTransitionMatrix.h>

#include "AutoCorrelationProcessPhyloLikelihood.h"
#include "HmmLikelihood_DF.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

AutoCorrelationProcessPhyloLikelihood::AutoCorrelationProcessPhyloLikelihood(
  shared_ptr<const AlignmentDataInterface> data,
  shared_ptr<AutoCorrelationSequenceEvolution> processSeqEvol,
  shared_ptr<CollectionNodes> collNodes,
  size_t nSeqEvol,
  size_t nData) :
  AbstractPhyloLikelihood(collNodes->context()),
  AbstractAlignedPhyloLikelihood(collNodes->context(), data->getNumberOfSites()),
  AbstractSingleDataPhyloLikelihood(
      collNodes->context(),
      data->getNumberOfSites(),
      (processSeqEvol->getSubstitutionProcessNumbers().size() != 0)
          ? processSeqEvol->substitutionProcess(processSeqEvol->getSubstitutionProcessNumbers()[0]).getNumberOfStates()
	  : 0,
      nData),
  AbstractSequencePhyloLikelihood(collNodes->context(), processSeqEvol, nData),
  AbstractParametrizable(""),
  MultiProcessSequencePhyloLikelihood(data, processSeqEvol, collNodes, nSeqEvol, nData),
  Hpep_(),
  hmm_()
{
  auto alphyl = make_shared<HmmPhyloAlphabet>(*this);

  Hpep_ = make_shared<HmmPhyloEmissionProbabilities>(alphyl);

  hmm_ = shared_ptr<HmmLikelihood_DF>(new HmmLikelihood_DF(context(), processSeqEvol->getHmmProcessAlphabet(), processSeqEvol->getHmmTransitionMatrix(), Hpep_));

  addParameters_(hmm_->getParameters());
}

void AutoCorrelationProcessPhyloLikelihood::setNamespace(const std::string& nameSpace)
{
  MultiProcessSequencePhyloLikelihood::setNamespace(nameSpace);
  hmm_->setNamespace(nameSpace);
}


void AutoCorrelationProcessPhyloLikelihood::fireParameterChanged(const ParameterList& parameters)
{
  MultiProcessSequencePhyloLikelihood::fireParameterChanged(parameters);

  hmm_->matchParametersValues(parameters);
}
