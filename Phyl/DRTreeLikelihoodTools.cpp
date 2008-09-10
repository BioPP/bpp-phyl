//
// File: DRTreeLikelihoodTools.cpp
// Created by: Julien Dutheil
// Created on: Mon Janv 17 09:56 2005
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#include "DRTreeLikelihoodTools.h"
#include <NumCalc/VectorTools.h>

using namespace bpp;

//-----------------------------------------------------------------------------------------

VVVdouble DRTreeLikelihoodTools::getPosteriorProbabilitiesForEachStateForEachRate(
  const DRTreeLikelihood & drl,
  int nodeId)
{
  unsigned int nSites   = drl.getLikelihoodData()->getNumberOfDistinctSites();
  unsigned int nClasses = drl.getNumberOfClasses();
  unsigned int nStates  = drl.getNumberOfStates();
  VVVdouble postProb(nSites);
  
  const DiscreteDistribution * rDist = drl.getRateDistribution();
  Vdouble rcProbs = rDist->getProbabilities();
  if(drl.getTree()->isLeaf(nodeId))
  {
    VVdouble larray = drl.getLikelihoodData()->getLeafLikelihoods(nodeId);
    for(unsigned int i = 0; i < nSites; i++)
    {
      VVdouble * postProb_i = & postProb[i];
      postProb_i->resize(nClasses);
      Vdouble * larray_i = & larray[i];
      for(unsigned int c = 0; c < nClasses; c++)
      {
        Vdouble * postProb_i_c = & (* postProb_i)[c];
        postProb_i_c->resize(nStates);
        double * rcProb = & rcProbs[c];
        for(unsigned int x = 0; x < nStates; x++)
        {
          (* postProb_i_c)[x] = (* larray_i)[x] * (* rcProb);
        }
      }
    }
  }
  else
  {
    VVVdouble larray;
    drl.computeLikelihoodAtNode(nodeId, larray);
    
    Vdouble likelihoods(nSites, 0);
    for(unsigned int i = 0; i < nSites; i++)
    {
      VVdouble * larray_i = & larray[i];
      for(unsigned int c = 0; c < nClasses; c++)
      {
        Vdouble * larray_i_c = & (* larray_i)[c];
        for(unsigned int s = 0; s < nStates; s++)
        {
          likelihoods[i] += (* larray_i_c)[s];
        }
      }
    }
    
    for(unsigned int i = 0; i < nSites; i++)
    {
      VVdouble * postProb_i = & postProb[i];
      postProb_i->resize(nClasses);
      VVdouble * larray_i = & larray[i];
      double likelihood = likelihoods[i];
      for(unsigned int c = 0; c < nClasses; c++)
      {
        Vdouble * postProb_i_c = & (* postProb_i)[c];
        postProb_i_c->resize(nStates);
        Vdouble * larray_i_c = & (* larray_i)[c];
        for(unsigned int x = 0; x < nStates; x++)
        {
          (* postProb_i_c)[x] = (* larray_i_c)[x] / likelihood;
        }
      }
    }
  }
  return postProb;
}

//-----------------------------------------------------------------------------------------

