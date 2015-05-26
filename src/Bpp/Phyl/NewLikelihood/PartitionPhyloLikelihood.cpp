//
// File: PartitionPhyloLikelihood.cpp
// Created by: Laurent Guéguen
// Created on: samedi 16 mai 2015, à 13h 54
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

#include "PartitionPhyloLikelihood.h"
#include "SingleProcessPhyloLikelihood.h"
#include "SingleRecursiveTreeLikelihoodCalculation.h"
#include "DoubleRecursiveTreeLikelihoodCalculation.h"

#include <Bpp/Seq/Container/SiteContainerTools.h>

using namespace std;
using namespace bpp;

/******************************************************************************/

PartitionPhyloLikelihood::PartitionPhyloLikelihood(
  PartitionSequenceEvolution& processSeqEvol,
  char recursivity,
  size_t nSeqEvol,
  bool verbose,
  bool patterns) :
  SequencePhyloLikelihood(processSeqEvol, nSeqEvol),
  SumOfDataPhyloLikelihood(),
  mSeqEvol_(processSeqEvol),
  vProcPos_()
{
  SubstitutionProcessCollection& processColl = processSeqEvol.getCollection();

  map<size_t, vector<size_t> >& mProcPos=processSeqEvol.getMapOfProcessSites();

  vProcPos_.resize(processSeqEvol.getNumberOfSites());

  map<size_t, SingleProcessPhyloLikelihood*> mSP;
  
  for (std::map<size_t, std::vector<size_t> >::iterator it=mProcPos.begin(); it!=mProcPos.end(); it++)
  {
    TreeLikelihoodCalculation* rt;

    if (recursivity=='S')
      rt = new SingleRecursiveTreeLikelihoodCalculation(
        &processColl.getSubstitutionProcess(it->first),
        it==mProcPos.begin(), patterns);
    else if (recursivity=='D')
      rt = new DoubleRecursiveTreeLikelihoodCalculation(        
        &processColl.getSubstitutionProcess(it->first),
        it==mProcPos.begin());
    else throw(Exception("PartitionPhyloLikelihood::PartitionPhyloLikelihood: unknown recursivity : " + recursivity));

    addSingleDataPhylolikelihood(it->first, new SingleProcessPhyloLikelihood(&processColl.getSubstitutionProcess(it->first), rt, it->first));
    
    for (size_t i = 0; i<it->second.size();i++)
    {
      vProcPos_[it->second[i]].nProc=it->first;
      vProcPos_[it->second[i]].pos=i;      
    }
  }
}


/******************************************************************************/

PartitionPhyloLikelihood::PartitionPhyloLikelihood(
  const SiteContainer& data,
  PartitionSequenceEvolution& processSeqEvol,
  char recursivity,
  size_t nSeqEvol,
  size_t nData,
  bool verbose,
  bool patterns) :
  SequencePhyloLikelihood(processSeqEvol, nSeqEvol),
  SumOfDataPhyloLikelihood(),
  mSeqEvol_(processSeqEvol),
  vProcPos_()
{
  if (data.getNumberOfSites()!=processSeqEvol.getNumberOfSites())
    throw BadIntegerException("PartitionPhyloLikelihood::PartitionPhyloLikelihood, data and sequence process lengths do not match.", (int)data.getNumberOfSites());

  SubstitutionProcessCollection& processColl = processSeqEvol.getCollection();

  map<size_t, vector<size_t> >& mProcPos=processSeqEvol.getMapOfProcessSites();

  vProcPos_.resize(processSeqEvol.getNumberOfSites());

  map<size_t, SingleProcessPhyloLikelihood*> mSP;
  
  for (std::map<size_t, std::vector<size_t> >::iterator it=mProcPos.begin(); it!=mProcPos.end(); it++)
  {
    TreeLikelihoodCalculation* rt;
    
    if (recursivity=='S')
      rt = new SingleRecursiveTreeLikelihoodCalculation(
        &processColl.getSubstitutionProcess(it->first), it==mProcPos.begin(), patterns);
    else if (recursivity=='D')
      rt = new DoubleRecursiveTreeLikelihoodCalculation(
        &processColl.getSubstitutionProcess(it->first), it==mProcPos.begin());
    else throw(Exception("PartitionPhyloLikelihood::PartitionPhyloLikelihood: unknown recursivity : " + recursivity));

    addSingleDataPhylolikelihood(it->first, new SingleProcessPhyloLikelihood(&processColl.getSubstitutionProcess(it->first), rt, it->first));
    
    
    for (size_t i = 0; i<it->second.size();i++)
    {
      vProcPos_[it->second[i]].nProc=it->first;
      vProcPos_[it->second[i]].pos=i;      
    }
  }

  
  setData(data, nData);
}

/******************************************************************************/

void PartitionPhyloLikelihood::setData(const SiteContainer& data, size_t nData)
{
  if (data.getNumberOfSites()!=mSeqEvol_.getNumberOfSites())
    throw BadIntegerException("PartitionPhyloLikelihood::PartitionPhyloLikelihood, data and sequence process lengths do not match.", (int)data.getNumberOfSites());

  SequencePhyloLikelihood::setData(data, nData);

  const map<size_t, vector<size_t> >& mProcPos=mSeqEvol_.getMapOfProcessSites();

  for (std::map<size_t, std::vector<size_t> >::const_iterator it=mProcPos.begin(); it!=mProcPos.end(); it++)
  {
    SiteContainer* st=SiteContainerTools::getSelectedSites(data, it->second);

    SumOfDataPhyloLikelihood::setData(*st, it->first);
  }
   
 }

/******************************************************************************/

double PartitionPhyloLikelihood::getLikelihoodForASite(size_t site) const
{
  return getSingleDataPhylolikelihood(vProcPos_[site].nProc)->getLikelihoodForASite(vProcPos_[site].pos);
}

/******************************************************************************/

double PartitionPhyloLikelihood::getDLogLikelihoodForASite(size_t site) const
{
  return getSingleDataPhylolikelihood(vProcPos_[site].nProc)->getDLogLikelihoodForASite(vProcPos_[site].pos);
}

/******************************************************************************/

double PartitionPhyloLikelihood::getD2LogLikelihoodForASite(size_t site) const
{
  return getSingleDataPhylolikelihood(vProcPos_[site].nProc)->getD2LogLikelihoodForASite(vProcPos_[site].pos);
}
