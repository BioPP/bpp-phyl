//
// File: SinglePhyloLikelihood.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 17 17:57:21 2003
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

#include "SinglePhyloLikelihood.h"

using namespace std;
using namespace bpp;
using namespace newlik;

/******************************************************************************/

SinglePhyloLikelihood::SinglePhyloLikelihood(
  SubstitutionProcess* process,
  TreeLikelihoodCalculation* tlComp,
  bool verbose) :
  AbstractParametrizable(""),
  tlComp_(tlComp),
  process_(process),
  computeFirstOrderDerivatives_(true),
  computeSecondOrderDerivatives_(true),
  minusLogLik_(1),
  verbose_(verbose)
{
  if (tlComp->getSubstitutionProcess() != process_.get())
    throw Exception("SinglePhyloLikelihood::SinglePhyloLikelihood Error :  given process must be the same as the one of TreeLikelihoodCalculation");
  
  // initialize parameters:
  addParameters_(process_->getIndependentParameters());

  tlComp_->computeTreeLikelihood();
  minusLogLik_ = - tlComp_->getLogLikelihood();
}

/******************************************************************************/

void SinglePhyloLikelihood::setParameters(const ParameterList& parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

void SinglePhyloLikelihood::fireParameterChanged(const ParameterList& params)
{
  process_->matchParametersValues(params);

  if (params.size() > 0)
  {
    tlComp_->computeTreeLikelihood();
    minusLogLik_ = - tlComp_->getLogLikelihood();
  }
}

/******************************************************************************/

Vdouble SinglePhyloLikelihood::getLikelihoodForEachSite() const
{
  Vdouble l(getNumberOfSites());
  for (size_t i = 0; i < l.size(); ++i)
  {
    l[i] = getLikelihoodForASite(i);
  }
  return l;
}

/******************************************************************************/

VVdouble SinglePhyloLikelihood::getLikelihoodForEachSiteForEachState() const
{
  VVdouble l(getNumberOfSites());
  for (size_t i = 0; i < l.size(); ++i)
  {
    Vdouble* l_i = &l[i];
    l_i->resize(getNumberOfStates());
    for (size_t x = 0; x < l_i->size(); ++x)
    {
      (*l_i)[x] = tlComp_->getLikelihoodForASiteForAState(i, static_cast<int>(x));
    }
  }
  return l;
}

/******************************************************************************/

VVdouble SinglePhyloLikelihood::getLikelihoodForEachSiteForEachClass() const
{
  VVdouble l(getNumberOfSites());
  for (size_t i = 0; i < l.size(); ++i)
  {
    Vdouble* l_i = &l[i];
    l_i->resize(getNumberOfClasses());
    for (size_t c = 0; c < l_i->size(); ++c)
    {
      (*l_i)[c] = tlComp_->getLikelihoodForASiteForAClass(i, c);
    }
  }
  return l;
}

/******************************************************************************/

VVVdouble SinglePhyloLikelihood::getLikelihoodForEachSiteForEachClassForEachState() const
{
  VVVdouble l(getNumberOfSites());
  for (size_t i = 0; i < l.size(); ++i)
  {
    VVdouble* l_i = &l[i];
    l_i->resize(getNumberOfClasses());
    for (size_t c = 0; c < l_i->size(); ++c)
    {
      Vdouble* l_ic = &(*l_i)[c];
      l_ic->resize(getNumberOfStates());
      for (size_t x = 0; x < l_ic->size(); ++x)
      {
        (*l_ic)[x] = tlComp_->getLikelihoodForASiteForAClassForAState(i, c, static_cast<int>(x));
      }
    }
  }
  return l;
}

/******************************************************************************/

ParameterList SinglePhyloLikelihood::getNonDerivableParameters() const
{
  ParameterList pl=getSubstitutionModelParameters();
  pl.addParameters(getRootFrequenciesParameters());
  pl.addParameters(getRateDistributionParameters());

  return pl;
}

/******************************************************************************/

VVdouble SinglePhyloLikelihood::getPosteriorProbabilitiesOfEachClass() const
{
  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  VVdouble pb = getLikelihoodForEachSiteForEachClass();
  Vdouble l = getLikelihoodForEachSite();
  for (size_t i = 0; i < nbSites; ++i)
  {
    for (size_t j = 0; j < nbClasses; ++j)
    {
      pb[i][j] = pb[i][j] * process_->getProbabilityForModel(j) / l[i];
    }
  }
  return pb;
}

/******************************************************************************/

vector<size_t> SinglePhyloLikelihood::getClassWithMaxPostProbOfEachSite() const
{
  size_t nbSites = getNumberOfSites();
  VVdouble l = getLikelihoodForEachSiteForEachClass();
  vector<size_t> classes(nbSites);
  for (size_t i = 0; i < nbSites; ++i)
  {
    classes[i] = VectorTools::whichMax<double>(l[i]);
  }
  return classes;
}

/******************************************************************************/

double SinglePhyloLikelihood::getValue() const
throw (Exception)
{
  if (!isInitialized())
    throw Exception("SinglePhyloLikelihood::getValue(). Instance is not initialized.");
  return minusLogLik_;
}

/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/

double SinglePhyloLikelihood::getFirstOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("SinglePhyloLikelihood::getFirstOrderDerivative().", variable);
  if (process_->hasTransitionProbabilitiesParameter(variable))
  {
    throw Exception("Derivatives are only implemented for branch length parameters.");
  }

  tlComp_->computeTreeDLikelihood(variable);
  return - tlComp_->getDLogLikelihood();
}

/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/

  double SinglePhyloLikelihood::getSecondOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("SinglePhyloLikelihood::getSecondOrderDerivative().", variable);
  if (process_->hasTransitionProbabilitiesParameter(variable))
  {
    throw Exception("Derivatives are only implemented for branch length parameters.");
  }

  tlComp_->computeTreeD2Likelihood(variable);
  return - tlComp_->getD2LogLikelihood();
}

/******************************************************************************/

