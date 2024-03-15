// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Seq/Container/SiteContainerTools.h>

#include "PartitionProcessPhyloLikelihood.h"
#include "PhyloLikelihoodContainer.h"
#include "SingleProcessPhyloLikelihood.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

PartitionProcessPhyloLikelihood::PartitionProcessPhyloLikelihood(
    Context& context,
    shared_ptr<PartitionSequenceEvolution> processSeqEvol,
    size_t nSeqEvol) :
  AbstractPhyloLikelihood(context),
  AbstractAlignedPhyloLikelihood(context, 0),
  AbstractSingleDataPhyloLikelihood(
	context, processSeqEvol->getNumberOfSites(),
       	(processSeqEvol->getSubstitutionProcessNumbers().size() != 0) ? processSeqEvol->substitutionProcess(processSeqEvol->getSubstitutionProcessNumbers()[0]).getNumberOfStates() : 0, 0),
  AbstractParametrizable(""),
  AbstractSequencePhyloLikelihood(context, processSeqEvol, nSeqEvol),
  AbstractPhyloLikelihoodSet(context, std::make_shared<PhyloLikelihoodContainer>(context, processSeqEvol->getCollection())),
  mSeqEvol_(processSeqEvol),
  vProcPos_(),
  mData_(),
  likCal_(new AlignedLikelihoodCalculation(context))
{
  // make new shared Parameters
  auto collNodes = pPhyloCont_->getCollectionNodes();

  //

  map<size_t, vector<size_t> >& mProcPos = processSeqEvol->mapOfProcessSites();

  vProcPos_.resize(processSeqEvol->getNumberOfSites());

  // Build the PhyloLikelihoodContainer
  auto pC = getPhyloContainer();

  for (const auto& it : mProcPos)
  {
    auto nProcess = it.first;

    auto l = std::make_shared<LikelihoodCalculationSingleProcess>(collNodes, nProcess);

    auto nPL = make_shared<SingleProcessPhyloLikelihood>(context, l, nProcess);

    pC->addPhyloLikelihood(nProcess, nPL);

    for (size_t i = 0; i < it.second.size(); ++i)
    {
      vProcPos_[it.second[i]].nProc = nProcess;
      vProcPos_[it.second[i]].pos = i;
    }

    addPhyloLikelihood(nProcess);
  }
}


/******************************************************************************/

PartitionProcessPhyloLikelihood::PartitionProcessPhyloLikelihood(
    Context& context,
    shared_ptr<const AlignmentDataInterface> data,
    shared_ptr<PartitionSequenceEvolution> processSeqEvol,
    size_t nSeqEvol,
    size_t nData) :
  AbstractPhyloLikelihood(context),
  AbstractAlignedPhyloLikelihood(context, data->getNumberOfSites()),
  AbstractSingleDataPhyloLikelihood(
	context, data->getNumberOfSites(),
       	(processSeqEvol->getSubstitutionProcessNumbers().size() != 0) ? processSeqEvol->substitutionProcess(processSeqEvol->getSubstitutionProcessNumbers()[0]).getNumberOfStates() : 0, nData),
  AbstractParametrizable(""),
  AbstractSequencePhyloLikelihood(context, processSeqEvol, nSeqEvol, nData),
  AbstractPhyloLikelihoodSet(context, std::make_shared<PhyloLikelihoodContainer>(context, processSeqEvol->getCollection())),
  mSeqEvol_(processSeqEvol),
  vProcPos_(),
  mData_(),
  likCal_(new AlignedLikelihoodCalculation(context))
{
  if (data->getNumberOfSites() != processSeqEvol->getNumberOfSites())
    throw BadIntegerException("PartitionProcessPhyloLikelihood::PartitionProcessPhyloLikelihood, data and sequence process lengths do not match.", static_cast<int>(data->getNumberOfSites()));

  // get new shared Parameters
  auto collNodes = pPhyloCont_->getCollectionNodes();


  map<size_t, vector<size_t> >& mProcPos = processSeqEvol->mapOfProcessSites();

  vProcPos_.resize(processSeqEvol->getNumberOfSites());

  // Build the PhyloLikelihoodContainer
  auto pC = getPhyloContainer();

  for (const auto& it:mProcPos)
  {
    auto nProcess = it.first;

    auto l = std::make_shared<LikelihoodCalculationSingleProcess>(collNodes, nProcess);

    auto nPL = make_shared<SingleProcessPhyloLikelihood>(context, l, nProcess, nData);

    pC->addPhyloLikelihood(nProcess, nPL);

    for (size_t i = 0; i < it.second.size(); ++i)
    {
      vProcPos_[it.second[i]].nProc = nProcess;
      vProcPos_[it.second[i]].pos = i;
    }

    addPhyloLikelihood(nProcess);
  }

  // set Data (will build calculations)
  setData(data, nData);
}

/******************************************************************************/

PartitionProcessPhyloLikelihood::PartitionProcessPhyloLikelihood(
    shared_ptr<const AlignmentDataInterface> data,
    shared_ptr<PartitionSequenceEvolution> processSeqEvol,
    std::shared_ptr<CollectionNodes> collNodes,
    size_t nSeqEvol,
    size_t nData) :
  AbstractPhyloLikelihood(collNodes->context()),
  AbstractAlignedPhyloLikelihood(collNodes->context(), data->getNumberOfSites()),
  AbstractSingleDataPhyloLikelihood(
	collNodes->context(), data->getNumberOfSites(),
       	(processSeqEvol->getSubstitutionProcessNumbers().size() != 0) ? processSeqEvol->substitutionProcess(processSeqEvol->getSubstitutionProcessNumbers()[0]).getNumberOfStates() : 0, nData),
  AbstractParametrizable(""),
  AbstractSequencePhyloLikelihood(collNodes->context(), processSeqEvol, nSeqEvol, nData),
  AbstractPhyloLikelihoodSet(collNodes->context(), std::make_shared<PhyloLikelihoodContainer>(collNodes->context(), collNodes)),
  mSeqEvol_(processSeqEvol),
  vProcPos_(),
  mData_(),
  likCal_(make_shared<AlignedLikelihoodCalculation>(collNodes->context()))
{
  if (data->getNumberOfSites() != processSeqEvol->getNumberOfSites())
    throw BadIntegerException("PartitionProcessPhyloLikelihood::PartitionProcessPhyloLikelihood, data and sequence process lengths do not match.", static_cast<int>(data->getNumberOfSites()));

  map<size_t, vector<size_t> >& mProcPos = processSeqEvol->mapOfProcessSites();

  vProcPos_.resize(processSeqEvol->getNumberOfSites());

  // Build the PhyloLikelihoodContainer
  auto pC = getPhyloContainer();

  for (const auto& it : mProcPos)
  {
    auto nProcess = it.first;

    auto l = std::make_shared<LikelihoodCalculationSingleProcess>(collNodes, nProcess);

    auto nPL = make_shared<SingleProcessPhyloLikelihood>(this->context(), l, nProcess, nData);

    pC->addPhyloLikelihood(nProcess, nPL);

    for (size_t i = 0; i < it.second.size(); ++i)
    {
      vProcPos_[it.second[i]].nProc = nProcess;
      vProcPos_[it.second[i]].pos = i;
    }

    addPhyloLikelihood(nProcess);
  }

  // set Data (will build calculations)
  setData(data, nData);
}

/******************************************************************************/

void PartitionProcessPhyloLikelihood::setData(
    shared_ptr<const AlignmentDataInterface> data,
    size_t nData)
{
  if (data->getNumberOfSites() != mSeqEvol_->getNumberOfSites())
    throw BadIntegerException("PartitionProcessPhyloLikelihood::PartitionProcessPhyloLikelihood, data and sequence process lengths do not match.", static_cast<int>(data->getNumberOfSites()));

  AbstractSequencePhyloLikelihood::setData(data, nData);

  const auto& mProcPos = mSeqEvol_->mapOfProcessSites();

  for (const auto& it : mProcPos)
  {
    shared_ptr<AlignmentDataInterface> st = SiteContainerTools::getSelectedSites(*data, it.second);
    getPhyloContainer()->setData(st, it.first);
    mData_[it.first] = st;
  }

  //  Then build calculations
  makeLikCal_();
}

/******************************************************************************/

void PartitionProcessPhyloLikelihood::makeLikCal_()
{
  NodeRefVec vLik;

  map<size_t, size_t> mProcInd; // map of the index of the process numbers in the list of vLikCal_

  // get the RowVectors of site likelihoods
  for (const auto& lik : vLikCal_)
  {
    vLik.push_back(dynamic_pointer_cast<AlignedLikelihoodCalculation>(lik)->getSiteLikelihoods(false));
  }

  // build the reverse map of Phylo indexes
  std::map<size_t, size_t> nProc2iProc;
  for (size_t i = 0; i < nPhylo_.size(); ++i)
  {
    nProc2iProc[nPhylo_[i]] = i;
  }

  // and the matching  between sequence & partition
  /*
   * Matrix of matching positions : site X (index of T in the vector of Ts, position for corresponding T)
   *
   */

  MatchingType matching;

  matching.resize(static_cast<Eigen::Index>(getNumberOfSites()), 2);

  for (size_t i = 0; i < getNumberOfSites(); ++i)
  {
    matching(Eigen::Index(i), 0) = nProc2iProc[vProcPos_[i].nProc];
    matching(Eigen::Index(i), 1) = vProcPos_[i].pos;
  }

  auto matchingDF = NumericConstant<MatchingType>::create(context(), matching);

  vLik.push_back(matchingDF);

  auto sL = CWiseMatching<RowLik, ReductionOf<RowLik> >::create(context(), std::move(vLik), RowVectorDimension (getNumberOfSites()));

  likCal_->setSiteLikelihoods(sL);

  auto lik = SumOfLogarithms<RowLik>::create(this->context(), {sL}, RowVectorDimension (getNumberOfSites()));

  likCal_->setLikelihoodNode(lik);

  // using bpp::DotOptions;
  // writeGraphToDot(
  //   "partition.dot", {lik.get()});//, DotOptions::DetailedNodeInfo | DotOp
}

