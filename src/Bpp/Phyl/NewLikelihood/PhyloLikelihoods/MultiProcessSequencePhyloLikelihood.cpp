//
// File: MultiProcessSequencePhyloLikelihood.cpp
// Created by: Laurent Guéguen
// Created on: vendredi 12 juillet 2013, à 00h 32
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

#include "MultiProcessSequencePhyloLikelihood.h"

#include "../DataFlow/CollectionNodes.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

MultiProcessSequencePhyloLikelihood::MultiProcessSequencePhyloLikelihood(
  Context& context,
  const AlignedValuesContainer& data,
  MultiProcessSequenceEvolution& processSeqEvol,
  size_t nSeqEvol,
  size_t nData,
  bool verbose,
  bool patterns) :
  AbstractPhyloLikelihood(context),
  AbstractAlignedPhyloLikelihood(context, data.getNumberOfSites()),
  AbstractSequencePhyloLikelihood(context, processSeqEvol, nSeqEvol, nData),
  mSeqEvol_(processSeqEvol),
  vLikCal_()
{
  resetParameters_(); // Do not keep the original parameters to get ConfiguredParameters
 // initialize parameters:

  const SubstitutionProcessCollection& processColl = processSeqEvol.getCollection();

  CollectionNodes cN(context, processColl);

  const vector<size_t>& nProc = processSeqEvol.getSubstitutionProcessNumbers();
  
  for (auto n:nProc)
  {
    auto liksing = std::make_shared<LikelihoodCalculationSingleProcess>(cN,
                                                                        data,
                                                                        n);
    liksing->makeLikelihoods();
    vLikCal_.push_back(liksing);
    shareParameters_(liksing->getParameters());
  }

}

/******************************************************************************/

void MultiProcessSequencePhyloLikelihood::setData(const AlignedValuesContainer& sites, size_t nData)
{
  AbstractSequencePhyloLikelihood::setData(sites, nData);
  
  for (auto& lik : vLikCal_)
    lik->setData(sites);
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
        (*l_i)[c] = getLikelihoodForASiteForAProcess(i, c);
      }
  }
  return l;
}

