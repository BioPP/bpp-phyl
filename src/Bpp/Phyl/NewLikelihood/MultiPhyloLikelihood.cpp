//
// File: MultiPhyloLikelihood.cpp
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

#include "MultiPhyloLikelihood.h"

#include "SingleRecursiveTreeLikelihoodCalculation.h"
#include "DoubleRecursiveTreeLikelihoodCalculation.h"

using namespace std;
using namespace bpp;
using namespace newlik;

/******************************************************************************/

MultiPhyloLikelihood::MultiPhyloLikelihood(
  SubstitutionProcessCollection* processColl,
  char recursivity,
  bool verbose,
  bool patterns) :
  AbstractParametrizable(""),
  data_(0),
  processColl_(processColl),
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
  
  addParameters_(processColl_->getIndependentParameters());

  size_t nbP = processColl_->getNumberOfSubstitutionProcess();

  if (recursivity_=='S')
    for (size_t i = 0; i < nbP; i++)
    {
      SingleRecursiveTreeLikelihoodCalculation* rt = new SingleRecursiveTreeLikelihoodCalculation(
        processColl_->getSubstitutionProcess(i), patterns, i == 0);
      vpTreelik_.push_back(rt);
    }
  else if (recursivity_=='D')
    for (size_t i = 0; i < nbP; i++)
    {
      DoubleRecursiveTreeLikelihoodCalculation* rt = new DoubleRecursiveTreeLikelihoodCalculation(
        processColl_->getSubstitutionProcess(i), i == 0);
      vpTreelik_.push_back(rt);
    }
  else throw(Exception("MultiPhyloLikelihood::MultiPhyloLikelihood: unknown recursivity : " + recursivity_));
  
}

/******************************************************************************/

MultiPhyloLikelihood::MultiPhyloLikelihood(
  const SiteContainer& data,
  SubstitutionProcessCollection* processColl,
  char recursivity,
  bool verbose,
  bool patterns) :
  AbstractParametrizable(""),
  data_(&data),
  processColl_(processColl),
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
  addParameters_(processColl_->getIndependentParameters());

  size_t nbP = processColl_->getNumberOfSubstitutionProcess();

  if (recursivity_=='S')
    for (size_t i = 0; i < nbP; i++)
    {
      SingleRecursiveTreeLikelihoodCalculation* rt = new SingleRecursiveTreeLikelihoodCalculation(*data_->clone(),
        processColl_->getSubstitutionProcess(i), patterns, i == 0);
      vpTreelik_.push_back(rt);
    }
  else if (recursivity_=='D')
    for (size_t i = 0; i < nbP; i++)
    {
      DoubleRecursiveTreeLikelihoodCalculation* rt = new DoubleRecursiveTreeLikelihoodCalculation(*data_->clone(),
        processColl_->getSubstitutionProcess(i), i == 0);
      vpTreelik_.push_back(rt);
    }
  else throw(Exception("MultiPhyloLikelihood::MultiPhyloLikelihood: unknown recursivity : " + recursivity_));
  
  setData(data);
}

/******************************************************************************/

void MultiPhyloLikelihood::setData(const SiteContainer& sites)
{
  for (size_t i = 0; i < vpTreelik_.size(); i++)
  {
    vpTreelik_[i]->setData(sites);
  }

  nbSites_ = sites.getNumberOfSites();
  nbStates_ = sites.getAlphabet()->getSize();

  initialized_ = true;
}

/******************************************************************************/

void MultiPhyloLikelihood::fireParameterChanged(const ParameterList& parameters)
{
  processColl_->matchParametersValues(parameters);

  if (parameters.size() > 0)
    for (size_t i = 0; i < vpTreelik_.size(); i++)
      vpTreelik_[i]->computeTreeLikelihood();
}

/******************************************************************************/

ParameterList MultiPhyloLikelihood::getSubstitutionProcessParameters() const
{
  return processColl_->getSubstitutionProcessParameters();
}

ParameterList MultiPhyloLikelihood::getSubstitutionModelParameters() const
{
  return processColl_->getSubstitutionModelParameters();
}

ParameterList MultiPhyloLikelihood::getRateDistributionParameters() const
{
  return processColl_->getRateDistributionParameters();
}

ParameterList MultiPhyloLikelihood::getRootFrequenciesParameters() const
{
  return processColl_->getRootFrequenciesParameters();
}

ParameterList MultiPhyloLikelihood::getBranchLengthsParameters() const
{
  return processColl_->getBranchLengthsParameters();
}

std::vector<const TreeTemplate<Node>* > MultiPhyloLikelihood::getTrees() const
{
  std::vector<const TreeTemplate<Node>* > vT;
  size_t i=1;
  while (processColl_->hasTreeNumber(i)){
    vT.push_back(&(processColl_->getTree(i)->getTree()));
    i++;
  }
  return vT;
}

/******************************************************************************/

Vdouble MultiPhyloLikelihood::getLikelihoodForEachSite() const
{
  Vdouble l(getNumberOfSites());
  for (unsigned int i = 0; i < l.size(); ++i)
  {
    l[i] = getLikelihoodForASite(i);
  }
  return l;
}

/******************************************************************************/

VVdouble MultiPhyloLikelihood::getLikelihoodForEachSiteForEachProcess() const
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


void MultiPhyloLikelihood::setParameters(const ParameterList& parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

double MultiPhyloLikelihood::getValue() const
throw (Exception)
{
  if (!isInitialized())
    throw Exception("MultiPhyloLikelihood::getValue(). Instance is not initialized.");
  return minusLogLik_;
}

/******************************************************************************/

double MultiPhyloLikelihood::getFirstOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("MultiPhyloLikelihood::getFirstOrderDerivative().", variable);
  if (!processColl_->hasBranchLengthsParameter(variable))
  {
    throw Exception("Derivatives are only implemented for branch length parameters.");
  }
  computeDLogLikelihood_(variable);
  return -getDLogLikelihood();
}

/******************************************************************************/

double MultiPhyloLikelihood::getSecondOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("MultiPhyloLikelihood::getSecondOrderDerivative().", variable);
  if (!processColl_->hasBranchLengthsParameter(variable))
  {
    throw Exception("Derivatives are only implemented for branch length parameters.");
  }
  computeD2LogLikelihood_(variable);
  return -getD2LogLikelihood();
}

/******************************************************************************/

void MultiPhyloLikelihood::computeDLogLikelihoodForAProcess(std::string& variable, size_t p) const
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

void MultiPhyloLikelihood::computeD2LogLikelihoodForAProcess(std::string& variable, size_t p) const
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
