// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Hmm/AutoCorrelationTransitionMatrix.h>

#include "AlignedPhyloLikelihoodAutoCorrelation.h"
#include "HmmLikelihood_DF.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

AlignedPhyloLikelihoodAutoCorrelation::AlignedPhyloLikelihoodAutoCorrelation(
    Context& context,
    std::shared_ptr<PhyloLikelihoodContainer> pC,
    const std::vector<size_t>& nPhylo,
    bool inCollection) :
  AbstractPhyloLikelihood(context),
  AbstractParametrizable(""),
  AbstractPhyloLikelihoodSet(context, pC, {}, inCollection),
  AbstractAlignedPhyloLikelihood(context, 0),
  AbstractAlignedPhyloLikelihoodSet(context, pC, nPhylo, inCollection),
  hma_(),
  htm_(),
  hpep_(),
  hmm_()
{
  hma_ = make_shared<HmmPhyloAlphabet>(*this);

  htm_ = make_shared<AutoCorrelationTransitionMatrix>(hma_, "AutoCorr.");

  hpep_ = make_shared<HmmPhyloEmissionProbabilities>(hma_);

  hmm_ = shared_ptr<HmmLikelihood_DF>(new HmmLikelihood_DF(this->context(), hma_, htm_, hpep_));

  addParameters_(htm_->getParameters());
}

void AlignedPhyloLikelihoodAutoCorrelation::setNamespace(const std::string& nameSpace)
{
  AbstractAlignedPhyloLikelihoodSet::setNamespace(nameSpace);
  hmm_->setNamespace(nameSpace);
}

void AlignedPhyloLikelihoodAutoCorrelation::fireParameterChanged(const ParameterList& parameters)
{
  AbstractAlignedPhyloLikelihoodSet::fireParameterChanged(parameters);
  hmm_->matchParametersValues(parameters);
  hma_->matchParametersValues(parameters);
  htm_->matchParametersValues(parameters);
  hpep_->matchParametersValues(parameters);
}
