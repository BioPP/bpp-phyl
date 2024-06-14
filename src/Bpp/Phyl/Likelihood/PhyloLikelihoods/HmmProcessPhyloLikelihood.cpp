// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Hmm/FullHmmTransitionMatrix.h>

#include "HmmLikelihood_DF.h"
#include "HmmProcessPhyloLikelihood.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

HmmProcessPhyloLikelihood::HmmProcessPhyloLikelihood(
    shared_ptr<const AlignmentDataInterface> data,
    shared_ptr<HmmSequenceEvolution> processSeqEvol,
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
  MultiProcessSequencePhyloLikelihood(data, processSeqEvol, collNodes, nData),
  Hpep_(),
  hmm_()
{
  auto alphyl = make_shared<HmmPhyloAlphabet>(*this);

  Hpep_ = make_shared<HmmPhyloEmissionProbabilities>(alphyl);

  hmm_ = make_shared<HmmLikelihood_DF>(context(), processSeqEvol->getHmmProcessAlphabet(), processSeqEvol->getHmmTransitionMatrix(), Hpep_);

  addParameters_(hmm_->getParameters());
}

void HmmProcessPhyloLikelihood::setNamespace(const std::string& nameSpace)
{
  MultiProcessSequencePhyloLikelihood::setNamespace(nameSpace);
  hmm_->setNamespace(nameSpace);
}

void HmmProcessPhyloLikelihood::fireParameterChanged(const ParameterList& parameters)
{
  MultiProcessSequencePhyloLikelihood::fireParameterChanged(parameters);
  hmm_->matchParametersValues(parameters);
}
