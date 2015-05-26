//
// File: OneProcessSequencePhyloLikelihood.cpp
// Created by: Julien Dutheil
// Created on: mardi 28 avril 2015, à 13h 11
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

#include "OneProcessSequencePhyloLikelihood.h"
#include "SingleRecursiveTreeLikelihoodCalculation.h"
#include "DoubleRecursiveTreeLikelihoodCalculation.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

OneProcessSequencePhyloLikelihood::OneProcessSequencePhyloLikelihood(
  OneProcessSequenceEvolution& evol,
  char recursivity,
  size_t nSeqEvol,
  bool verbose,
  bool patterns) :
  AbstractSequencePhyloLikelihood(evol, nSeqEvol),
  mSeqEvol_(evol),
  tlComp_()
{
  SubstitutionProcess& sp=evol.getSubstitutionProcess();

  if (recursivity=='S')
    tlComp_ = std::auto_ptr<TreeLikelihoodCalculation>(new SingleRecursiveTreeLikelihoodCalculation(&sp, true, patterns));
  else
    if (recursivity=='D')
      tlComp_ = std::auto_ptr<TreeLikelihoodCalculation>(new DoubleRecursiveTreeLikelihoodCalculation(&sp, true));
    else throw(Exception("OneProcessPhyloLikelihood::OneProcessPhyloLikelihood: unknown recursivity : " + recursivity));

  tlComp_->resetToCompute();
}

/******************************************************************************/

OneProcessSequencePhyloLikelihood::OneProcessSequencePhyloLikelihood(
  const SiteContainer& data,
  OneProcessSequenceEvolution& evol,
  char recursivity,
  size_t nSeqEvol,
  size_t nData,
  bool verbose,
  bool patterns) :
  AbstractSequencePhyloLikelihood(evol, nSeqEvol),
  mSeqEvol_(evol),
  tlComp_()
{
  SubstitutionProcess& sp=evol.getSubstitutionProcess();

  if (recursivity=='S')
    tlComp_ = std::auto_ptr<TreeLikelihoodCalculation>(new SingleRecursiveTreeLikelihoodCalculation(&sp, true, patterns));
  else
    if (recursivity=='D')
      tlComp_ = std::auto_ptr<TreeLikelihoodCalculation>(new DoubleRecursiveTreeLikelihoodCalculation(&sp, true));
    else throw(Exception("OneProcessPhyloLikelihood::OneProcessPhyloLikelihood: unknown recursivity : " + recursivity));

  tlComp_->resetToCompute();

  setData(data, nData);
}



/******************************************************************************/

void OneProcessSequencePhyloLikelihood::fireParameterChanged(const ParameterList& params)
{
  AbstractSequencePhyloLikelihood::fireParameterChanged(params);
  
  tlComp_->resetToCompute();
}

/******************************************************************************/

VVdouble OneProcessSequencePhyloLikelihood::getLikelihoodForEachSiteForEachState() const
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

VVdouble OneProcessSequencePhyloLikelihood::getLikelihoodForEachSiteForEachClass() const
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

VVVdouble OneProcessSequencePhyloLikelihood::getLikelihoodForEachSiteForEachClassForEachState() const
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

VVdouble OneProcessSequencePhyloLikelihood::getPosteriorProbabilitiesOfEachClass() const
{
  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  VVdouble pb = getLikelihoodForEachSiteForEachClass();
  Vdouble l = getLikelihoodForEachSite();

  for (size_t i = 0; i < nbSites; ++i)
  {
    for (size_t j = 0; j < nbClasses; ++j)
    {
      pb[i][j] = pb[i][j] * getSubstitutionProcess().getProbabilityForModel(j) / l[i];
    }
  }
  return pb;
}

/******************************************************************************/

vector<size_t> OneProcessSequencePhyloLikelihood::getClassWithMaxPostProbOfEachSite() const
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

Vdouble OneProcessSequencePhyloLikelihood::getPosteriorRateOfEachSite() const
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
      rates[i] += (pb[i][j] / l[i]) * getSubstitutionProcess().getProbabilityForModel(j) *  getSubstitutionProcess().getProbabilityForModel(j);
    }
  }
  return rates;
}

/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/

double OneProcessSequencePhyloLikelihood::getFirstOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    return 0;

  if (!hasDerivableParameter(variable))
  {
    throw Exception("Non derivable parameter: " + variable);
  }

  tlComp_->computeTreeDLogLikelihood(variable);
  return - tlComp_->getDLogLikelihood();
}

/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/

double OneProcessSequencePhyloLikelihood::getSecondOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    return 0;
  
  if (!hasDerivableParameter(variable))
  {
    throw Exception("Non derivable parameter: " + variable);
  }

  tlComp_->computeTreeD2LogLikelihood(variable);
  return - tlComp_->getD2LogLikelihood();
}

/******************************************************************************/

