//
// File: MultiProcessPhyloLikelihood.cpp
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

#include "MultiProcessPhyloLikelihood.h"

#include "SingleRecursiveTreeLikelihoodCalculation.h"
#include "DoubleRecursiveTreeLikelihoodCalculation.h"

using namespace std;
using namespace bpp;
using namespace newlik;

/******************************************************************************/

MultiProcessPhyloLikelihood::MultiProcessPhyloLikelihood(
  SubstitutionProcessCollection* processColl,
  vector<size_t> nProc,
  char recursivity,
  bool verbose,
  bool patterns) :
  AbstractSingleDataPhyloLikelihood(""),
  processColl_(processColl),
  nProc_(nProc),
  recursivity_('S'),
  computeFirstOrderDerivatives_(true),
  computeSecondOrderDerivatives_(true),
  initialized_(false),
  verbose_(verbose),
  nbSites_(0),
  nbStates_(0),
  minusLogLik_(0),
  vpTreelik_()
{
  // initialize parameters:

  for (size_t i=0; i<nProc_.size(); i++)
    includeParameters_(processColl_->getSubstitutionProcessParameters(nProc_[i]));
  
  if (recursivity_=='S')
    for (size_t i = 0; i < nProc_.size(); i++)
    {
      SingleRecursiveTreeLikelihoodCalculation* rt = new SingleRecursiveTreeLikelihoodCalculation(
        processColl_->getSubstitutionProcess(nProc_[i]), patterns, i == 0);
      vpTreelik_.push_back(rt);
    }
  else if (recursivity_=='D')
    for (size_t i = 0; i < nProc_.size(); i++)
    {
      DoubleRecursiveTreeLikelihoodCalculation* rt = new DoubleRecursiveTreeLikelihoodCalculation(
        processColl_->getSubstitutionProcess(nProc_[i]), i == 0);
      vpTreelik_.push_back(rt);
    }
  else throw(Exception("MultiProcessPhyloLikelihood::MultiProcessPhyloLikelihood: unknown recursivity : " + recursivity_));
  
}

/******************************************************************************/

MultiProcessPhyloLikelihood::MultiProcessPhyloLikelihood(
  const SiteContainer& data,
  SubstitutionProcessCollection* processColl,
  vector<size_t> nProc,
  char recursivity,
  size_t nData,
  bool verbose,
  bool patterns) :
  AbstractSingleDataPhyloLikelihood("", nData),
  processColl_(processColl),
  nProc_(nProc),
  recursivity_('S'),
  computeFirstOrderDerivatives_(true),
  computeSecondOrderDerivatives_(true),
  initialized_(false),
  verbose_(verbose),
  nbSites_(0),
  nbStates_(0),
  minusLogLik_(0),
  vpTreelik_()
{
  // initialize parameters:
  for (size_t i=0; i<nProc_.size(); i++)
    includeParameters_(processColl_->getSubstitutionProcessParameters(nProc_[i]));

  if (recursivity_=='S')
    for (size_t i = 0; i < nProc_.size(); i++)
    {
      SingleRecursiveTreeLikelihoodCalculation* rt = new SingleRecursiveTreeLikelihoodCalculation(data, processColl_->getSubstitutionProcess(nProc_[i]), patterns, i == 0);
      vpTreelik_.push_back(rt);
    }
  else if (recursivity_=='D')
    for (size_t i = 0; i < nProc_.size(); i++)
    {
      DoubleRecursiveTreeLikelihoodCalculation* rt = new DoubleRecursiveTreeLikelihoodCalculation(data, processColl_->getSubstitutionProcess(nProc_[i]), i == 0);
      vpTreelik_.push_back(rt);
    }
  else throw(Exception("MultiProcessPhyloLikelihood::MultiProcessPhyloLikelihood: unknown recursivity : " + recursivity_));

  setData(data, nData);
}

/******************************************************************************/

void MultiProcessPhyloLikelihood::setData(const SiteContainer& sites, size_t nData)
{
  for (size_t i = 0; i < vpTreelik_.size(); i++)
  {
    vpTreelik_[i]->setData(sites);
  }

  setNData(nData);
  
  nbSites_ = sites.getNumberOfSites();
  nbStates_ = sites.getAlphabet()->getSize();

  initialized_ = true;
}

/******************************************************************************/

void MultiProcessPhyloLikelihood::fireParameterChanged(const ParameterList& parameters)
{
  processColl_->matchParametersValues(parameters);

  if (parameters.size() > 0)
    for (size_t i = 0; i < vpTreelik_.size(); i++)
      vpTreelik_[i]->computeTreeLikelihood();
}

/******************************************************************************/

ParameterList MultiProcessPhyloLikelihood::getSubstitutionProcessParameters() const
{
  ParameterList pl;
  
  for (size_t i=0; i<nProc_.size(); i++)
    pl.includeParameters(processColl_->getSubstitutionProcessParameters(nProc_[i]));

  return pl;
}

ParameterList MultiProcessPhyloLikelihood::getSubstitutionModelParameters() const
{ 
  ParameterList pl;
  
  for (size_t i=0; i<nProc_.size(); i++)
    pl.includeParameters(processColl_->getSubstitutionProcess(nProc_[i])->getSubstitutionModelParameters());

  return pl;
}

ParameterList MultiProcessPhyloLikelihood::getRateDistributionParameters() const
{
  ParameterList pl;
  
  for (size_t i=0; i<nProc_.size(); i++)
    pl.includeParameters(processColl_->getSubstitutionProcess(nProc_[i])->getRateDistributionParameters());

  return pl;
}

ParameterList MultiProcessPhyloLikelihood::getRootFrequenciesParameters() const
{
  ParameterList pl;
  
  for (size_t i=0; i<nProc_.size(); i++)
    pl.includeParameters(processColl_->getSubstitutionProcess(nProc_[i])->getRootFrequenciesParameters());

  return pl;
}

ParameterList MultiProcessPhyloLikelihood::getBranchLengthParameters() const
{ 
  ParameterList pl;
  
  for (size_t i=0; i<nProc_.size(); i++)
    pl.includeParameters(processColl_->getSubstitutionProcess(nProc_[i])->getBranchLengthParameters());

  return pl;
}

ParameterList MultiProcessPhyloLikelihood::getNonDerivableParameters() const
{ 
  ParameterList pl;
  
  for (size_t i=0; i<nProc_.size(); i++)
    pl.includeParameters(processColl_->getSubstitutionProcess(nProc_[i])->getNonDerivableParameters());

  return pl;
}

/******************************************************************************/

Vdouble MultiProcessPhyloLikelihood::getLikelihoodForEachSite() const
{
  Vdouble l(getNumberOfSites());
  for (unsigned int i = 0; i < l.size(); ++i)
  {
    l[i] = getLikelihoodForASite(i);
  }
  return l;
}

/******************************************************************************/

VVdouble MultiProcessPhyloLikelihood::getLikelihoodForEachSiteForEachProcess() const
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


void MultiProcessPhyloLikelihood::setParameters(const ParameterList& parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

double MultiProcessPhyloLikelihood::getValue() const
throw (Exception)
{
  if (!isInitialized())
    throw Exception("MultiProcessPhyloLikelihood::getValue(). Instance is not initialized.");
  return minusLogLik_;
}

/******************************************************************************/

double MultiProcessPhyloLikelihood::getFirstOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("MultiProcessPhyloLikelihood::getFirstOrderDerivative().", variable);
  if (!processColl_->hasBranchLengthParameter(variable))
  {
    throw Exception("Derivatives are only implemented for branch length parameters.");
  }
  computeDLogLikelihood_(variable);
  return -getDLogLikelihood();
}

/******************************************************************************/

double MultiProcessPhyloLikelihood::getSecondOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("MultiProcessPhyloLikelihood::getSecondOrderDerivative().", variable);
  if (!processColl_->hasBranchLengthParameter(variable))
  {
    throw Exception("Derivatives are only implemented for branch length parameters.");
  }
  computeD2LogLikelihood_(variable);
  return -getD2LogLikelihood();
}

/******************************************************************************/

void MultiProcessPhyloLikelihood::computeDLogLikelihoodForAProcess(std::string& variable, size_t p) const
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
  
  vector<string> valias= processColl_->getAlias(variable);
  
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

void MultiProcessPhyloLikelihood::computeD2LogLikelihoodForAProcess(std::string& variable, size_t p) const
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
  
  vector<string> valias= processColl_->getAlias(variable);
  
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
