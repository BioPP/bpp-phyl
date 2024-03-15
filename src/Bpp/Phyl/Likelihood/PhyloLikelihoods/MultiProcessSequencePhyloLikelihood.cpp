// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../DataFlow/CollectionNodes.h"
#include "MultiProcessSequencePhyloLikelihood.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

MultiProcessSequencePhyloLikelihood::MultiProcessSequencePhyloLikelihood(
  shared_ptr<const AlignmentDataInterface> data,
  shared_ptr<MultiProcessSequenceEvolution> processSeqEvol,
  shared_ptr<CollectionNodes> collNodes,
  size_t nSeqEvol,
  size_t nData) :
  AbstractPhyloLikelihood(collNodes->context()),
  AbstractAlignedPhyloLikelihood(collNodes->context(), data->getNumberOfSites()),
  AbstractSequencePhyloLikelihood(collNodes->context(), processSeqEvol, nSeqEvol, nData),
  AbstractParametrizableSequencePhyloLikelihood(collNodes->context(), processSeqEvol, nSeqEvol),
  mSeqEvol_(processSeqEvol),
  vLikCal_()
{
  resetParameters_(); // Do not keep the original parameters to get ConfiguredParameters
  // initialize parameters:

  const vector<size_t>& nProc = processSeqEvol->getSubstitutionProcessNumbers();

  for (auto n:nProc)
  {
    auto liksing = std::make_shared<LikelihoodCalculationSingleProcess>(collNodes,
                                                                        data,
                                                                        n);
    liksing->makeLikelihoods();
    vLikCal_.push_back(liksing);
    shareParameters_(liksing->getParameters());
  }
}

/******************************************************************************/

void MultiProcessSequencePhyloLikelihood::setData(shared_ptr<const AlignmentDataInterface> sites, size_t nData)
{
  AbstractSequencePhyloLikelihood::setData(sites, nData);

  for (auto& lik : vLikCal_)
  {
    lik->setData(sites);
  }
}

/******************************************************************************/

VVdouble MultiProcessSequencePhyloLikelihood::getLikelihoodPerSitePerProcess() const
{
  VVdouble l(getNumberOfSites());
  for (size_t i = 0; i < l.size(); ++i)
  {
    Vdouble* l_i = &l[i];
    l_i->resize(getNumberOfSubstitutionProcess());
    for (size_t c = 0; c < l_i->size(); ++c)
    {
      (*l_i)[c] = convert(getLikelihoodForASiteForAProcess(i, c));
    }
  }
  return l;
}
