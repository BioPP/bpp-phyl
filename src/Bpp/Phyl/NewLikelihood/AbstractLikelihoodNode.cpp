//
// File: AbstractLikelihoodNode.cpp
// Created by: Laurent Guéguen
// Created on: mercredi 1 juillet 2015, à 16h 22
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

#include "AbstractLikelihoodNode.h"
#include "ComputingNode.h"
#include "../PatternTools.h"

using namespace bpp;
using namespace std;


void AbstractLikelihoodNode::getPosteriorProbabilitiesForEachState(VVdouble& vPP) const
{
  size_t nSites=nodeLikelihoods_.size();
  size_t nStates=nodeLikelihoods_[0].size();

  vPP.resize(nSites);
  
  const VVdouble& larray = getLikelihoodArray(ComputingNode::D0);

  for (size_t i = 0; i < nSites; i++)
  {
    const Vdouble * larray_i = & larray[i];
    Vdouble * vPP_i = & vPP[i];
    vPP_i->resize(nStates);
    
    double likelihood;
    if (usesLog())
    {
      likelihood=VectorTools::logSumExp(* larray_i);

      for(size_t s = 0; s < nStates; s++)
        (* vPP_i)[s] = exp((* larray_i)[s] - likelihood);
    }
    else
    {
      likelihood=VectorTools::sum((* larray_i));
      
      for(size_t s = 0; s < nStates; s++)
        (* vPP_i)[s] = (* larray_i)[s] / likelihood;
    }
  }
}

void AbstractLikelihoodNode::setUseLog(bool useLog)
{
  if (useLog==usesLog())
    return;

  if (isUp2date(ComputingNode::D0))
  {
    size_t nSites=nodeLikelihoods_.size();
    size_t nStates=nodeLikelihoods_[0].size();

    if (useLog)
    {
      for (size_t i = 0; i < nSites; i++)
      {
        Vdouble* nodeLikelihoods_i_ = &(nodeLikelihoods_[i]);
        for(size_t s = 0; s < nStates; s++)
          (*nodeLikelihoods_i_)[s]=log((*nodeLikelihoods_i_)[s]);
      }
    }
    for (size_t i = 0; i < nSites; i++)
    {
      Vdouble* nodeLikelihoods_i_ = &(nodeLikelihoods_[i]);
      for(size_t s = 0; s < nStates; s++)
        (*nodeLikelihoods_i_)[s]=exp((*nodeLikelihoods_i_)[s]);
    }
  }
  
  usesLog_=useLog;
}


void AbstractLikelihoodNode::setUseLogDownward(bool useLog)
{
  setUseLog(useLog);

  size_t nS=getNumberOfSons();
  for (size_t i=0; i<nS; i++)
    static_cast<AbstractLikelihoodNode*>(getSon(i))->setUseLogDownward(useLog);
}
    
