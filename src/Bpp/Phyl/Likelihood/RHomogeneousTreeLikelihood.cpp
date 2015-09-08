//
// File: RHomogeneousTreeLikelihood.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 17 18:14:51 2003
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

#include "RHomogeneousTreeLikelihood.h"
#include "../PatternTools.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

RHomogeneousTreeLikelihood::RHomogeneousTreeLikelihood(
  const Tree& tree,
  SubstitutionModel* model,
  DiscreteDistribution* rDist,
  bool checkRooted,
  bool verbose,
  bool usePatterns)
throw (Exception) :
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  likelihoodData_(0),
  minusLogLik_(-1.)
{
  init_(usePatterns);
}

/******************************************************************************/

RHomogeneousTreeLikelihood::RHomogeneousTreeLikelihood(
  const Tree& tree,
  const SiteContainer& data,
  SubstitutionModel* model,
  DiscreteDistribution* rDist,
  bool checkRooted,
  bool verbose,
  bool usePatterns)
throw (Exception) :
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  likelihoodData_(0),
  minusLogLik_(-1.)
{
  init_(usePatterns);
  setData(data);
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::init_(bool usePatterns) throw (Exception)
{
  likelihoodData_ = new DRASRTreeLikelihoodData(
    tree_,
    rateDistribution_->getNumberOfCategories(),
    usePatterns);
}

/******************************************************************************/

RHomogeneousTreeLikelihood::RHomogeneousTreeLikelihood(
  const RHomogeneousTreeLikelihood& lik) :
  AbstractHomogeneousTreeLikelihood(lik),
  likelihoodData_(0),
  minusLogLik_(lik.minusLogLik_)
{
  likelihoodData_ = dynamic_cast<DRASRTreeLikelihoodData*>(lik.likelihoodData_->clone());
  likelihoodData_->setTree(tree_);
}

/******************************************************************************/

RHomogeneousTreeLikelihood& RHomogeneousTreeLikelihood::operator=(
  const RHomogeneousTreeLikelihood& lik)
{
  AbstractHomogeneousTreeLikelihood::operator=(lik);
  if (likelihoodData_) delete likelihoodData_;
  likelihoodData_ = dynamic_cast<DRASRTreeLikelihoodData*>(lik.likelihoodData_->clone());
  likelihoodData_->setTree(tree_);
  minusLogLik_ = lik.minusLogLik_;
  return *this;
}

/******************************************************************************/

RHomogeneousTreeLikelihood::~RHomogeneousTreeLikelihood()
{
  delete likelihoodData_;
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::setData(const SiteContainer& sites) throw (Exception)
{
  if (data_) delete data_;
  data_ = PatternTools::getSequenceSubset(sites, *tree_->getRootNode());
  if (verbose_) ApplicationTools::displayTask("Initializing data structure");
  likelihoodData_->initLikelihoods(*data_, *model_);
  if (verbose_) ApplicationTools::displayTaskDone();

  nbSites_ = likelihoodData_->getNumberOfSites();
  nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();
  nbStates_ = likelihoodData_->getNumberOfStates();

  if (verbose_) ApplicationTools::displayResult("Number of distinct sites",
                                                TextTools::toString(nbDistinctSites_));
  initialized_ = false;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLikelihood() const
{
  double l = 1.;
  for (size_t i = 0; i < nbSites_; i++)
  {
    l *= getLikelihoodForASite(i);
  }
  return l;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLogLikelihood() const
{
  double ll = 0;
  vector<double> la(nbSites_);
  for (size_t i = 0; i < nbSites_; i++)
  {
    la[i] = getLogLikelihoodForASite(i);
  }
  sort(la.begin(), la.end());
  for (size_t i = nbSites_; i > 0; i--)
  {
    ll += la[i - 1];
  }
  return ll;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLikelihoodForASite(size_t site) const
{
  double l = 0;
  for (size_t i = 0; i < nbClasses_; i++)
  {
    l += getLikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
  }
  return l;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLogLikelihoodForASite(size_t site) const
{
  double l = 0;
  for (size_t i = 0; i < nbClasses_; i++)
  {
    double li = getLikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
    if (li > 0) l+= li; //Corrects for numerical instabilities leading to slightly negative likelihoods
  }
  return log(l);
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  double l = 0;
  Vdouble* la = &likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
  for (size_t i = 0; i < nbStates_; i++)
  {
    //cout << (*la)[i] << "\t" << rootFreqs_[i] << endl;
    double li = (*la)[i] * rootFreqs_[i];
    if (li > 0) l+= li; //Corrects for numerical instabilities leading to slightly negative likelihoods
  }
  return l;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  double l = 0;
  Vdouble* la = &likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
  for (size_t i = 0; i < nbStates_; i++)
  {
    l += (*la)[i] * rootFreqs_[i];
  }
  //if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
  return log(l);
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  return likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass][static_cast<size_t>(state)];
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  return log(likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass][static_cast<size_t>(state)]);
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::setParameters(const ParameterList& parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::fireParameterChanged(const ParameterList& params)
{
  applyParameters();

  if (rateDistribution_->getParameters().getCommonParametersWith(params).size() > 0
      || model_->getParameters().getCommonParametersWith(params).size() > 0)
  {
    //Rate parameter changed, need to recompute all probs:
    computeAllTransitionProbabilities();
  }
  else if (params.size() > 0)
  {
    //We may save some computations:
    for (size_t i = 0; i < params.size(); i++)
    {
      string s = params[i].getName();
      if (s.substr(0, 5) == "BrLen")
      {
        //Branch length parameter:
        computeTransitionProbabilitiesForNode(nodes_[TextTools::to<size_t>(s.substr(5))]);
      }
    }
    rootFreqs_ = model_->getFrequencies();
  }

  computeTreeLikelihood();

  minusLogLik_ = -getLogLikelihood();
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getValue() const
throw (Exception)
{
  if (!isInitialized()) throw Exception("RHomogeneousTreeLikelihood::getValue(). Instance is not initialized.");
  return minusLogLik_;
}

/******************************************************************************
*                           First Order Derivatives                          *
******************************************************************************/

double RHomogeneousTreeLikelihood::getDLikelihoodForASiteForARateClass(
  size_t site,
  size_t rateClass) const
{
  double dl = 0;
  Vdouble* dla = &likelihoodData_->getDLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
  for (size_t i = 0; i < nbStates_; i++)
  {
    dl += (*dla)[i] * rootFreqs_[i];
  }
  return dl;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getDLikelihoodForASite(size_t site) const
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t i = 0; i < nbClasses_; i++)
  {
    dl += getDLikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
  }
  return dl;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getDLogLikelihoodForASite(size_t site) const
{
  // d(f(g(x)))/dx = dg(x)/dx . df(g(x))/dg :
  return getDLikelihoodForASite(site) / getLikelihoodForASite(site);
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getDLogLikelihood() const
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t i = 0; i < nbSites_; i++)
  {
    dl += getDLogLikelihoodForASite(i);
  }

  return dl;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getFirstOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("RHomogeneousTreeLikelihood::getFirstOrderDerivative().", variable);
  if (getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameter are not implemented.");
  }
  if (getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }

  const_cast<RHomogeneousTreeLikelihood*>(this)->computeTreeDLikelihood(variable);
  return -getDLogLikelihood();
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::computeTreeDLikelihood(const string& variable)
{
  // Get the node with the branch whose length must be derivated:
  size_t brI = TextTools::to<size_t>(variable.substr(5));
  const Node* branch = nodes_[brI];
  const Node* father = branch->getFather();
  VVVdouble* _dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(father->getId());

  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  size_t nbSites  = _dLikelihoods_father->size();
  for (size_t i = 0; i < nbSites; i++)
  {
    VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_dLikelihoods_father_i_c)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* son = father->getSon(l);

    vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
    VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

    if (son == branch)
    {
      VVVdouble* dpxy__son = &dpxy_[son->getId()];

      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
          Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
          VVdouble* dpxy__son_c = &(*dpxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* dpxy__son_c_x = &(*dpxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*dpxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
            }
            (*_dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
    else
    {
      VVVdouble* pxy__son = &pxy_[son->getId()];
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
          Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
          VVdouble* pxy__son_c = &(*pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
            }
            (*_dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
  }

  // Now we go down the tree toward the root node:
  computeDownSubtreeDLikelihood(father);
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::computeDownSubtreeDLikelihood(const Node* node)
{
  const Node* father = node->getFather();
  // We assume that the _dLikelihoods array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if (father == NULL) return; // We reached the root!

  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  VVVdouble* _dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(father->getId());
  size_t nbSites  = _dLikelihoods_father->size();
  for (size_t i = 0; i < nbSites; i++)
  {
    VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_dLikelihoods_father_i_c)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* son = father->getSon(l);
    VVVdouble* pxy__son = &pxy_[son->getId()];
    
    vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());

    if (son == node)
    {
      VVVdouble* _dLikelihoods_son = &likelihoodData_->getDLikelihoodArray(son->getId());
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _dLikelihoods_son_i = &(*_dLikelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _dLikelihoods_son_i_c = &(*_dLikelihoods_son_i)[c];
          Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
          VVdouble* pxy__son_c = &(*pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*pxy__son_c_x)[y] * (*_dLikelihoods_son_i_c)[y];
            }
            (*_dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
    else
    {
      VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
          Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
          VVdouble* pxy__son_c = &(*pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
            }
            (*_dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
  }

  //Next step: move toward grand father...
  computeDownSubtreeDLikelihood(father);
}

/******************************************************************************
*                           Second Order Derivatives                         *
******************************************************************************/

double RHomogeneousTreeLikelihood::getD2LikelihoodForASiteForARateClass(
  size_t site,
  size_t rateClass) const
{
  double d2l = 0;
  Vdouble* d2la = &likelihoodData_->getD2LikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
  for (size_t i = 0; i < nbStates_; i++)
  {
    d2l += (*d2la)[i] * rootFreqs_[i];
  }
  return d2l;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getD2LikelihoodForASite(size_t site) const
{
  // Derivative of the sum is the sum of derivatives:
  double d2l = 0;
  for (size_t i = 0; i < nbClasses_; i++)
  {
    d2l += getD2LikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
  }
  return d2l;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getD2LogLikelihoodForASite(size_t site) const
{
  return getD2LikelihoodForASite(site) / getLikelihoodForASite(site)
         - pow( getDLikelihoodForASite(site) / getLikelihoodForASite(site), 2);
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getD2LogLikelihood() const
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t i = 0; i < nbSites_; i++)
  {
    dl += getD2LogLikelihoodForASite(i);
  }
  return dl;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getSecondOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("RHomogeneousTreeLikelihood::getSecondOrderDerivative().", variable);
  if (getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameter are not implemented.");
  }
  if (getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }

  const_cast<RHomogeneousTreeLikelihood*>(this)->computeTreeD2Likelihood(variable);
  return -getD2LogLikelihood();
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::computeTreeD2Likelihood(const string& variable)
{
  // Get the node with the branch whose length must be derivated:
  size_t brI = TextTools::to<size_t>(variable.substr(5));
  const Node* branch = nodes_[brI];
  const Node* father = branch->getFather();

  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  VVVdouble* _d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
  size_t nbSites  = _d2Likelihoods_father->size();
  for (size_t i = 0; i < nbSites; i++)
  {
    VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_d2Likelihoods_father_i_c)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* son = father->getSon(l);

    vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
    VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

    if (son == branch)
    {
      VVVdouble* d2pxy__son = &d2pxy_[son->getId()];
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
          Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
          VVdouble* d2pxy__son_c = &(*d2pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double d2l = 0;
            Vdouble* d2pxy__son_c_x = &(*d2pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              d2l += (*d2pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
            }
            (*_d2Likelihoods_father_i_c)[x] *= d2l;
          }
        }
      }
    }
    else
    {
      VVVdouble* pxy__son = &pxy_[son->getId()];
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
          Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
          VVdouble* pxy__son_c = &(*pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double d2l = 0;
            Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              d2l += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
            }
            (*_d2Likelihoods_father_i_c)[x] *= d2l;
          }
        }
      }
    }
  }

  // Now we go down the tree toward the root node:
  computeDownSubtreeD2Likelihood(father);
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::computeDownSubtreeD2Likelihood(const Node* node)
{
  const Node* father = node->getFather();
  // We assume that the _dLikelihoods array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if (father == NULL) return; // We reached the root!

  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  VVVdouble* _d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
  size_t nbSites  = _d2Likelihoods_father->size();
  for (size_t i = 0; i < nbSites; i++)
  {
    VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_d2Likelihoods_father_i_c)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* son = father->getSon(l);

    VVVdouble* pxy__son = &pxy_[son->getId()];
    vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());

    if (son == node)
    {
      VVVdouble* _d2Likelihoods_son = &likelihoodData_->getD2LikelihoodArray(son->getId());
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _d2Likelihoods_son_i = &(*_d2Likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _d2Likelihoods_son_i_c = &(*_d2Likelihoods_son_i)[c];
          Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
          VVdouble* pxy__son_c = &(*pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double d2l = 0;
            Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              d2l += (*pxy__son_c_x)[y] * (*_d2Likelihoods_son_i_c)[y];
            }
            (*_d2Likelihoods_father_i_c)[x] *= d2l;
          }
        }
      }
    }
    else
    {
      VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
        VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
          Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
          VVdouble* pxy__son_c = &(*pxy__son)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
            for (size_t y = 0; y < nbStates_; y++)
            {
              dl += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
            }
            (*_d2Likelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
  }

  //Next step: move toward grand father...
  computeDownSubtreeD2Likelihood(father);
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::computeTreeLikelihood()
{
  computeSubtreeLikelihood(tree_->getRootNode());
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::computeSubtreeLikelihood(const Node* node)
{
  if (node->isLeaf()) return;

  size_t nbSites = likelihoodData_->getLikelihoodArray(node->getId()).size();
  size_t nbNodes = node->getNumberOfSons();

  // Must reset the likelihood array first (i.e. set all of them to 1):
  VVVdouble* _likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());
  for (size_t i = 0; i < nbSites; i++)
  {
    //For each site in the sequence,
    VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
    for (size_t c = 0; c < nbClasses_; c++)
    {
      //For each rate classe,
      Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
      for (size_t x = 0; x < nbStates_; x++)
      {
        //For each initial state,
        (*_likelihoods_node_i_c)[x] = 1.;
      }
    }
  }

  for (size_t l = 0; l < nbNodes; l++)
  {
    //For each son node,

    const Node* son = node->getSon(l);

    computeSubtreeLikelihood(son); //Recursive method:

    VVVdouble* pxy__son = &pxy_[son->getId()];
    vector<size_t> * _patternLinks_node_son = &likelihoodData_->getArrayPositions(node->getId(), son->getId());
    VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

    for (size_t i = 0; i < nbSites; i++)
    {
      //For each site in the sequence,
      VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_node_son)[i]];
      VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
      for (size_t c = 0; c < nbClasses_; c++)
      {
        //For each rate classe,
        Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
        Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
        VVdouble* pxy__son_c = &(*pxy__son)[c];
        for (size_t x = 0; x < nbStates_; x++)
        {
          //For each initial state,
          Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
          double likelihood = 0;
          for (size_t y = 0; y < nbStates_; y++)
            likelihood += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
          
          (*_likelihoods_node_i_c)[x] *= likelihood;
        }
      }
    }
  }
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::displayLikelihood(const Node* node)
{
  cout << "Likelihoods at node " << node->getName() << ": " << endl;
  displayLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId()));
  cout << "                                         ***" << endl;
}

/*******************************************************************************/

