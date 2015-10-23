//
// File: SingleProcessPhyloLikelihood.cpp
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

#include "SingleProcessPhyloLikelihood.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

SingleProcessPhyloLikelihood::SingleProcessPhyloLikelihood(
  SubstitutionProcess* process,
  LikelihoodTreeCalculation* tlComp,
  size_t nProc,
  size_t nData) :
  AbstractPhyloLikelihood(),
  AbstractAlignedPhyloLikelihood(tlComp->getNumberOfSites()),
  AbstractSingleDataPhyloLikelihood(tlComp->getNumberOfSites(), process->getNumberOfStates(), nData),
  AbstractParametrizable(""),
  tlComp_(tlComp),
  process_(process),
  nProc_(nProc)
{
  if (tlComp->getSubstitutionProcess() != process_)
    throw Exception("SingleProcessPhyloLikelihood::SingleProcessPhyloLikelihood Error :  given process must be the same as the one of LikelihoodTreeCalculation");

  // initialize INDEPENDENT parameters:

  addParameters_(process_->getBranchLengthParameters(true));
  addParameters_(process_->getSubstitutionModelParameters(true));
  addParameters_(process_->getRateDistributionParameters(true));
  addParameters_(process_->getRootFrequenciesParameters(true));
}

/******************************************************************************/

void SingleProcessPhyloLikelihood::fireParameterChanged(const ParameterList& params)
{
  update();
  process_->matchParametersValues(params);
}

/******************************************************************************/

VVdouble SingleProcessPhyloLikelihood::getLikelihoodForEachSiteForEachState() const
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

VVdouble SingleProcessPhyloLikelihood::getLikelihoodForEachSiteForEachClass() const
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

VVVdouble SingleProcessPhyloLikelihood::getLikelihoodForEachSiteForEachClassForEachState() const
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

ParameterList SingleProcessPhyloLikelihood::getNonDerivableParameters() const
{
  // patch, to be fixed properly later
  return getParameters();

  ParameterList pl = getSubstitutionModelParameters();
  pl.addParameters(getRootFrequenciesParameters());
  pl.addParameters(getRateDistributionParameters());

  return pl;
}

/******************************************************************************/

VVdouble SingleProcessPhyloLikelihood::getPosteriorProbabilitiesOfEachClass() const
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

vector<size_t> SingleProcessPhyloLikelihood::getClassWithMaxPostProbOfEachSite() const
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

Vdouble SingleProcessPhyloLikelihood::getPosteriorRateOfEachSite() const
{
  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  VVdouble pb = getLikelihoodForEachSiteForEachClass();
  Vdouble l  = getLikelihoodForEachSite();
  Vdouble rates(nbSites, 0.);
  for (size_t i = 0; i < nbSites; i++)
  {
    for (size_t j = 0; j < nbClasses; j++)
    {
      rates[i] += (pb[i][j] / l[i]) * process_->getProbabilityForModel(j) *  process_->getProbabilityForModel(j);
    }
  }
  return rates;
}

/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/

void SingleProcessPhyloLikelihood::computeDLogLikelihood_(const string& variable) const
{
  // patch, to be fixed properly later
  throw Exception("SingleProcessPhyloLikelihood::1sd Derivative is not implemented for " + variable + " parameter.");
  
 tlComp_->computeTreeDLogLikelihood(variable);
}

/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/

void SingleProcessPhyloLikelihood::computeD2LogLikelihood_(const string& variable) const
{
  // patch, to be fixed properly later
  throw Exception("SingleProcessPhyloLikelihood::2nd Derivative is not implemented for " + variable + " parameter.");

  tlComp_->computeTreeD2LogLikelihood(variable);
}

/******************************************************************************/

