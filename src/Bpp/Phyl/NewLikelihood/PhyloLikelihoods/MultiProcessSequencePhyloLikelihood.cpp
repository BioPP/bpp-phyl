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

#include "../RecursiveLikelihoodTreeCalculation.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

MultiProcessSequencePhyloLikelihood::MultiProcessSequencePhyloLikelihood(
  const SiteContainer& data,
  MultiProcessSequenceEvolution& processSeqEvol,
  size_t nSeqEvol,
  size_t nData,
  bool verbose,
  bool patterns) :
  AbstractPhyloLikelihood(),
  AbstractAlignedPhyloLikelihood(data.getNumberOfSites()),
  AbstractSequencePhyloLikelihood(processSeqEvol, nSeqEvol, nData),
  mSeqEvol_(processSeqEvol),
  vpTreelik_()
{
  // initialize parameters:

  const SubstitutionProcessCollection& processColl = processSeqEvol.getCollection();

  const vector<size_t>& nProc = processSeqEvol.getSubstitutionProcessNumbers();
  
  for (size_t i = 0; i < nProc.size(); i++)
  {
    RecursiveLikelihoodTreeCalculation* rt = new RecursiveLikelihoodTreeCalculation(
      &processColl.getSubstitutionProcess(nProc[i]), i == 0, patterns);
    vpTreelik_.push_back(rt);
  }

  setData(data, nData);
}

/******************************************************************************/

void MultiProcessSequencePhyloLikelihood::setData(const SiteContainer& sites, size_t nData)
{
  AbstractSequencePhyloLikelihood::setData(sites, nData);
  
  for (size_t i = 0; i < vpTreelik_.size(); i++)
  {
    vpTreelik_[i]->setData(sites);
  }
  updateLikelihood();
  computeLikelihood();
}

/******************************************************************************/

VVdouble MultiProcessSequencePhyloLikelihood::getLikelihoodForEachSiteForEachProcess() const
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


/******************************************************************************/

void MultiProcessSequencePhyloLikelihood::computeDLogLikelihoodForAProcess(std::string& variable, size_t p) const
{
  size_t i=0;
  
  try {
    i=(size_t)atoi(variable.substr(variable.rfind('_')+1).c_str());
    if (p+1==i){
      vpTreelik_[p]->computeTreeDLogLikelihood(variable);
      return;
    }
  }
  catch (exception& e){}
  
  vector<string> valias= mSeqEvol_.getCollection().getAlias(variable);
  
  for (size_t v=0; v<valias.size();v++)
  {
    try {
      i=(size_t)atoi(valias[v].substr(valias[v].rfind('_')+1).c_str());
      if (p+1==i){
        vpTreelik_[p]->computeTreeDLogLikelihood(valias[v]);
        return;
      }
    }
    catch (exception& e){
      continue;
    }
  }
  
  vpTreelik_[p]->computeTreeDLogLikelihood("");

}


/************************************************************/

void MultiProcessSequencePhyloLikelihood::computeD2LogLikelihoodForAProcess(std::string& variable, size_t p) const
{
  size_t i=0;
  
  try {
    i=(size_t)atoi(variable.substr(variable.rfind('_')+1).c_str());
    if (p+1==i){
      vpTreelik_[p]->computeTreeD2LogLikelihood(variable);
      return;
    }
  }
  catch (exception& e){}
  
  vector<string> valias= mSeqEvol_.getCollection().getAlias(variable);
  
  for (size_t v=0; v<valias.size();v++)
  {
    try {
      i=(size_t)atoi(valias[v].substr(valias[v].rfind('_')+1).c_str());
      if (p+1==i){
        vpTreelik_[p]->computeTreeD2LogLikelihood(valias[v]);
        return;
      }
    }
    catch (exception& e){
      continue;
    }
  }
  
  vpTreelik_[p]->computeTreeD2LogLikelihood("");

}
