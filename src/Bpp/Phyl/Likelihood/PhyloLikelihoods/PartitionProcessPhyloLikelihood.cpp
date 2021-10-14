//
// File: PartitionProcessPhyloLikelihood.cpp
// Authors:
//   Laurent GuÃÂ©guen
// Created: samedi 16 mai 2015, ÃÂ  13h 54
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

#include <Bpp/Seq/Container/SiteContainerTools.h>

#include "PartitionProcessPhyloLikelihood.h"
#include "PhyloLikelihoodContainer.h"
#include "SingleProcessPhyloLikelihood.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

PartitionProcessPhyloLikelihood::PartitionProcessPhyloLikelihood(
  Context& context,
  PartitionSequenceEvolution& processSeqEvol,
  size_t nSeqEvol,
  bool verbose,
  bool patterns) :
  AbstractPhyloLikelihood(context),
  AbstractAlignedPhyloLikelihood(context, 0),
  SequencePhyloLikelihood(context, processSeqEvol, nSeqEvol),
  SetOfAbstractPhyloLikelihood(context, std::make_shared<PhyloLikelihoodContainer>(context, processSeqEvol.getCollection())),
  mSeqEvol_(processSeqEvol),
  vProcPos_(),
  mData_(),
  likCal_(new AlignedLikelihoodCalculation(context))
{
  // make new shared Parameters
  auto collNodes = pPhyloCont_->getCollectionNodes();

  //

  map<size_t, vector<size_t> >& mProcPos = processSeqEvol.getMapOfProcessSites();

  vProcPos_.resize(processSeqEvol.getNumberOfSites());

  // Build the PhyloLikelihoodContainer
  auto pC = getPhyloContainer();

  for (const auto& it:mProcPos)
  {
    auto nProcess = it.first;

    auto l = std::make_shared<LikelihoodCalculationSingleProcess>(*collNodes, nProcess);

    auto nPL = make_shared<SingleProcessPhyloLikelihood>(context, l, nProcess);

    pC->sharePhyloLikelihood(nProcess, nPL);

    for (size_t i = 0; i < it.second.size(); i++)
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
  const AlignedValuesContainer& data,
  PartitionSequenceEvolution& processSeqEvol,
  size_t nSeqEvol,
  size_t nData,
  bool verbose,
  bool patterns) :
  AbstractPhyloLikelihood(context),
  AbstractAlignedPhyloLikelihood(context, data.getNumberOfSites()),
  SequencePhyloLikelihood(context, processSeqEvol, nSeqEvol, nData),
  SetOfAbstractPhyloLikelihood(context, std::make_shared<PhyloLikelihoodContainer>(context, processSeqEvol.getCollection())),
  mSeqEvol_(processSeqEvol),
  vProcPos_(),
  mData_(),
  likCal_(new AlignedLikelihoodCalculation(context))
{
  if (data.getNumberOfSites() != processSeqEvol.getNumberOfSites())
    throw BadIntegerException("PartitionProcessPhyloLikelihood::PartitionProcessPhyloLikelihood, data and sequence process lengths do not match.", (int)data.getNumberOfSites());

  // get new shared Parameters
  auto collNodes = pPhyloCont_->getCollectionNodes();


  map<size_t, vector<size_t> >& mProcPos = processSeqEvol.getMapOfProcessSites();

  vProcPos_.resize(processSeqEvol.getNumberOfSites());

  // Build the PhyloLikelihoodContainer
  auto pC = getPhyloContainer();

  for (const auto& it:mProcPos)
  {
    auto nProcess = it.first;

    auto l = std::make_shared<LikelihoodCalculationSingleProcess>(*collNodes, nProcess);

    auto nPL = make_shared<SingleProcessPhyloLikelihood>(context, l, nProcess, nData);

    pC->sharePhyloLikelihood(nProcess, nPL);

    for (size_t i = 0; i < it.second.size(); i++)
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
  const AlignedValuesContainer& data,
  PartitionSequenceEvolution& processSeqEvol,
  std::shared_ptr<CollectionNodes> collNodes,
  size_t nSeqEvol,
  size_t nData,
  bool verbose,
  bool patterns) :
  AbstractPhyloLikelihood(collNodes->getContext()),
  AbstractAlignedPhyloLikelihood(collNodes->getContext(), data.getNumberOfSites()),
  SequencePhyloLikelihood(collNodes->getContext(), processSeqEvol, nSeqEvol, nData),
  SetOfAbstractPhyloLikelihood(collNodes->getContext(), std::make_shared<PhyloLikelihoodContainer>(collNodes->getContext(), collNodes)),
  mSeqEvol_(processSeqEvol),
  vProcPos_(),
  mData_(),
  likCal_(new AlignedLikelihoodCalculation(collNodes->getContext()))
{
  if (data.getNumberOfSites() != processSeqEvol.getNumberOfSites())
    throw BadIntegerException("PartitionProcessPhyloLikelihood::PartitionProcessPhyloLikelihood, data and sequence process lengths do not match.", (int)data.getNumberOfSites());

  map<size_t, vector<size_t> >& mProcPos = processSeqEvol.getMapOfProcessSites();

  vProcPos_.resize(processSeqEvol.getNumberOfSites());

  // Build the PhyloLikelihoodContainer
  auto pC = getPhyloContainer();

  for (const auto& it:mProcPos)
  {
    auto nProcess = it.first;

    auto l = std::make_shared<LikelihoodCalculationSingleProcess>(*collNodes, nProcess);

    auto nPL = make_shared<SingleProcessPhyloLikelihood>(getContext(), l, nProcess, nData);

    pC->sharePhyloLikelihood(nProcess, nPL);

    for (size_t i = 0; i < it.second.size(); i++)
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

void PartitionProcessPhyloLikelihood::setData(const AlignedValuesContainer& data, size_t nData)
{
  if (data.getNumberOfSites() != mSeqEvol_.getNumberOfSites())
    throw BadIntegerException("PartitionProcessPhyloLikelihood::PartitionProcessPhyloLikelihood, data and sequence process lengths do not match.", (int)data.getNumberOfSites());

  SequencePhyloLikelihood::setData(data, nData);

  const auto& mProcPos = mSeqEvol_.getMapOfProcessSites();

  for (const auto& it:mProcPos)
  {
    auto st = shared_ptr<AlignedValuesContainer>(SiteContainerTools::getSelectedSites(data, it.second));
    mData_[it.first] = st;
    getPhyloContainer()->setData(*st, it.first);
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
  for (const auto& lik:vLikCal_)
  {
    vLik.push_back(dynamic_cast<LikelihoodCalculationSingleProcess*>(lik.get())->getSiteLikelihoods(false));
  }

  // build the reverse map of Phylo indexes
  std::map<size_t, size_t> nProc2iProc;
  for (size_t i = 0; i < nPhylo_.size(); i++)
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

  for (size_t i = 0; i < getNumberOfSites(); i++)
  {
    matching(Eigen::Index(i), 0) = nProc2iProc[vProcPos_[i].nProc];
    matching(Eigen::Index(i), 1) = vProcPos_[i].pos;
  }

  auto matchingDF = NumericConstant<MatchingType>::create(getContext(), matching);

  vLik.push_back(matchingDF);

  auto sL = CWiseMatching<RowLik, ReductionOf<RowLik> >::create(getContext(), std::move(vLik), RowVectorDimension (getNumberOfSites()));

  likCal_->setSiteLikelihoods(sL);

  auto lik = SumOfLogarithms<RowLik>::create (getContext(), {sL}, RowVectorDimension (getNumberOfSites()));

  likCal_->setLikelihoodNode(lik);

  // using bpp::DotOptions;
  // writeGraphToDot(
  //   "partition.dot", {lik.get()});//, DotOptions::DetailedNodeInfo | DotOp
}
