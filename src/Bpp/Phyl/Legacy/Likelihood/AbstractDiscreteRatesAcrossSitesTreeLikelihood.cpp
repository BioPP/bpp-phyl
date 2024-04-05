// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/VectorTools.h>

#include "AbstractDiscreteRatesAcrossSitesTreeLikelihood.h"

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

AbstractDiscreteRatesAcrossSitesTreeLikelihood::AbstractDiscreteRatesAcrossSitesTreeLikelihood(
    shared_ptr<DiscreteDistributionInterface> rDist,
    bool verbose) :
  rateDistribution_(rDist)
{
  AbstractTreeLikelihood::enableDerivatives(true);
}

/******************************************************************************/

ParameterList AbstractDiscreteRatesAcrossSitesTreeLikelihood::getRateDistributionParameters() const
{
  if (!initialized_)
    throw Exception("AbstractDiscreteRatesAcrossSitesTreeLikelihood::getRateDistributionParameters(). Object is not initialized.");
  return rateDistribution_->getParameters().getCommonParametersWith(getParameters());
}

/******************************************************************************/

ParameterList AbstractDiscreteRatesAcrossSitesTreeLikelihood::getDerivableParameters() const
{
  if (!initialized_)
    throw Exception("AbstractDiscreteRatesAcrossSitesTreeLikelihood::getDerivableParameters(). Object is not initialized.");
  return getBranchLengthsParameters();
}

/******************************************************************************/

ParameterList AbstractDiscreteRatesAcrossSitesTreeLikelihood::getNonDerivableParameters() const
{
  if (!initialized_)
    throw Exception("AbstractDiscreteRatesAcrossSitesTreeLikelihood::getNonDerivableParameters(). Object is not initialized.");
  ParameterList tmp = getSubstitutionModelParameters();
  tmp.addParameters(getRateDistributionParameters());
  return tmp;
}

/******************************************************************************/

VVdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getLikelihoodPerSitePerRateClass() const
{
  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  VVdouble l(nbSites);
  for (size_t i = 0; i < nbSites; i++)
  {
    l[i].resize(nbClasses);
    for (size_t j = 0; j < nbClasses; j++)
    {
      l[i][j] = getLikelihoodForASiteForARateClass(i, j);
    }
  }
  return l;
}

/******************************************************************************/

double AbstractDiscreteRatesAcrossSitesTreeLikelihood::getLikelihoodForASiteForAState(size_t site, int state) const
{
  size_t nbClasses = getNumberOfClasses();
  double l = 0;
  for (size_t i = 0; i < nbClasses; i++)
  {
    l += getLikelihoodForASiteForARateClassForAState(site, i, state) * rateDistribution_->getProbability(i);
  }
  return l;
}

/******************************************************************************/

double AbstractDiscreteRatesAcrossSitesTreeLikelihood::getLogLikelihoodForASiteForAState(size_t site, int state) const
{
  size_t nbClasses = getNumberOfClasses();
  double l = 0;
  for (size_t i = 0; i < nbClasses; i++)
  {
    l += getLikelihoodForASiteForARateClassForAState(site, i, state) * rateDistribution_->getProbability(i);
  }
  // if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
  return log(l);
}

/******************************************************************************/

VVdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getLogLikelihoodPerSitePerRateClass() const
{
  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  VVdouble l(nbSites);
  for (size_t i = 0; i < nbSites; i++)
  {
    l[i] = Vdouble(nbClasses);
    for (size_t j = 0; j < nbClasses; j++)
    {
      l[i][j] = getLogLikelihoodForASiteForARateClass(i, j);
    }
  }
  return l;
}

/******************************************************************************/

VVVdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getLikelihoodPerSitePerRateClassPerState() const
{
  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  size_t nbStates  = getNumberOfStates();
  VVVdouble l(nbSites);
  for (size_t i = 0; i < nbSites; i++)
  {
    l[i].resize(nbClasses);
    for (size_t j = 0; j < nbClasses; j++)
    {
      l[i][j].resize(nbStates);
      for (size_t x = 0; x < nbStates; x++)
      {
        l[i][j][x] = getLikelihoodForASiteForARateClassForAState(i, j, static_cast<int>(x));
      }
    }
  }
  return l;
}

/******************************************************************************/

VVVdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getLogLikelihoodPerSitePerRateClassPerState() const
{
  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  size_t nbStates  = getNumberOfStates();
  VVVdouble l(nbSites);
  for (size_t i = 0; i < nbSites; i++)
  {
    l[i].resize(nbClasses);
    for (size_t j = 0; j < nbClasses; j++)
    {
      l[i][j].resize(nbStates);
      for (size_t x = 0; x < nbStates; x++)
      {
        l[i][j][x] = getLogLikelihoodForASiteForARateClassForAState(i, j, static_cast<int>(x));
      }
    }
  }
  return l;
}

/*******************************************************************************/

VVdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getPosteriorProbabilitiesPerRate() const
{
  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  VVdouble pb = getLikelihoodPerSitePerRateClass();
  Vdouble l  = getLikelihoodPerSite();
  for (size_t i = 0; i < nbSites; i++)
  {
    for (size_t j = 0; j < nbClasses; j++)
    {
      pb[i][j] = pb[i][j] * rateDistribution_->getProbability(j) / l[i];
    }
  }
  return pb;
}

/******************************************************************************/

Vdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getPosteriorRatePerSite() const
{
  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  VVdouble lr = getLikelihoodPerSitePerRateClass();
  Vdouble l  = getLikelihoodPerSite();
  Vdouble rates(nbSites, 0.);
  for (size_t i = 0; i < nbSites; i++)
  {
    for (size_t j = 0; j < nbClasses; j++)
    {
      rates[i] += (lr[i][j] / l[i]) * rateDistribution_->getProbability(j) *  rateDistribution_->getCategory(j);
    }
  }
  return rates;
}

/******************************************************************************/

vector<size_t> AbstractDiscreteRatesAcrossSitesTreeLikelihood::getRateClassWithMaxPostProbPerSite() const
{
  size_t nbSites = getNumberOfSites();
  VVdouble l = getLikelihoodPerSitePerRateClass();
  vector<size_t> classes(nbSites);
  for (size_t i = 0; i < nbSites; i++)
  {
    classes[i] = VectorTools::whichMax<double>(l[i]);
  }
  return classes;
}

/******************************************************************************/

Vdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getRateWithMaxPostProbPerSite() const
{
  size_t nbSites = getNumberOfSites();
  VVdouble l = getLikelihoodPerSitePerRateClass();
  Vdouble rates(nbSites);
  for (size_t i = 0; i < nbSites; i++)
  {
    rates[i] = rateDistribution_->getCategory(VectorTools::whichMax<double>(l[i]));
  }
  return rates;
}

/******************************************************************************/

void AbstractDiscreteRatesAcrossSitesTreeLikelihood::resetLikelihoodArray(
    VVVdouble& likelihoodArray)
{
  size_t nbSites   = likelihoodArray.size();
  size_t nbClasses = likelihoodArray[0].size();
  size_t nbStates  = likelihoodArray[0][0].size();
  for (size_t i = 0; i < nbSites; i++)
  {
    for (size_t c = 0; c < nbClasses; c++)
    {
      for (size_t s = 0; s < nbStates; s++)
      {
        likelihoodArray[i][c][s] = 1.;
      }
    }
  }
}

/******************************************************************************/

void AbstractDiscreteRatesAcrossSitesTreeLikelihood::displayLikelihoodArray(
    const VVVdouble& likelihoodArray)
{
  size_t nbSites   = likelihoodArray.size();
  size_t nbClasses = likelihoodArray[0].size();
  size_t nbStates  = likelihoodArray[0][0].size();
  for (size_t i = 0; i < nbSites; i++)
  {
    cout << "Site " << i << ":" << endl;
    for (size_t c = 0; c < nbClasses; c++)
    {
      cout << "Rate class " << c;
      for (size_t s = 0; s < nbStates; s++)
      {
        cout << "\t" << likelihoodArray[i][c][s];
      }
      cout << endl;
    }
    cout << endl;
  }
}

/******************************************************************************/

VVdouble AbstractDiscreteRatesAcrossSitesTreeLikelihood::getTransitionProbabilities(int nodeId, size_t siteIndex) const
{
  VVVdouble p3 = getTransitionProbabilitiesPerRateClass(nodeId, siteIndex);
  VVdouble p2;
  Vdouble probs = rateDistribution_->getProbabilities();
  p2.resize(getNumberOfStates());
  for (size_t i = 0; i < p2.size(); i++)
  {
    p2[i].resize(getNumberOfStates());
    for (size_t j = 0; j < p2.size(); j++)
    {
      for (size_t k = 0; k < getNumberOfClasses(); k++)
      {
        p2[i][j] += p3[k][i][j] * probs[k];
      }
    }
  }
  return p2;
}

/******************************************************************************/
