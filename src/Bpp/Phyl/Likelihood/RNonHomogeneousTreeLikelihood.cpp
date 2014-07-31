//
// File: RNonHomogeneousTreeLikelihood.cpp
// Created by: Julien Dutheil
// Created on: Tue Oct 09 16:07 2007
// From file: RHomogeneousTreeLikelihood.cpp
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

#include "RNonHomogeneousTreeLikelihood.h"
#include "../PatternTools.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

RNonHomogeneousTreeLikelihood::RNonHomogeneousTreeLikelihood(
  const Tree& tree,
  SubstitutionModelSet* modelSet,
  DiscreteDistribution* rDist,
  bool verbose,
  bool usePatterns,
  bool reparametrizeRoot)
throw (Exception) :
  AbstractNonHomogeneousTreeLikelihood(tree, modelSet, rDist, verbose, reparametrizeRoot),
  likelihoodData_(0),
  minusLogLik_(-1.)
{
  if (!modelSet->isFullySetUpFor(tree))
    throw Exception("RNonHomogeneousTreeLikelihood(constructor). Model set is not fully specified.");
  init_(usePatterns);
}

/******************************************************************************/

RNonHomogeneousTreeLikelihood::RNonHomogeneousTreeLikelihood(
  const Tree& tree,
  const SiteContainer& data,
  SubstitutionModelSet* modelSet,
  DiscreteDistribution* rDist,
  bool verbose,
  bool usePatterns,
  bool reparametrizeRoot)
throw (Exception) :
  AbstractNonHomogeneousTreeLikelihood(tree, modelSet, rDist, verbose, reparametrizeRoot),
  likelihoodData_(0),
  minusLogLik_(-1.)
{
  if (!modelSet->isFullySetUpFor(tree))
    throw Exception("RNonHomogeneousTreeLikelihood(constructor). Model set is not fully specified.");
  init_(usePatterns);
  setData(data);
}

/******************************************************************************/

void RNonHomogeneousTreeLikelihood::init_(bool usePatterns) throw (Exception)
{
  likelihoodData_ = new DRASRTreeLikelihoodData(
    tree_,
    rateDistribution_->getNumberOfCategories(),
    usePatterns);
}

/******************************************************************************/

RNonHomogeneousTreeLikelihood::RNonHomogeneousTreeLikelihood(
  const RNonHomogeneousTreeLikelihood& lik) :
  AbstractNonHomogeneousTreeLikelihood(lik),
  likelihoodData_(0),
  minusLogLik_(lik.minusLogLik_)
{
  likelihoodData_ = dynamic_cast<DRASRTreeLikelihoodData*>(lik.likelihoodData_->clone());
  likelihoodData_->setTree(tree_);
}

/******************************************************************************/

RNonHomogeneousTreeLikelihood& RNonHomogeneousTreeLikelihood::operator=(
  const RNonHomogeneousTreeLikelihood& lik)
{
  AbstractNonHomogeneousTreeLikelihood::operator=(lik);
  if (likelihoodData_) delete likelihoodData_;
  likelihoodData_ = dynamic_cast<DRASRTreeLikelihoodData*>(lik.likelihoodData_->clone());
  likelihoodData_->setTree(tree_);
  minusLogLik_ = lik.minusLogLik_;
  return *this;
}

/******************************************************************************/

RNonHomogeneousTreeLikelihood::~RNonHomogeneousTreeLikelihood()
{
  delete likelihoodData_;
}

/******************************************************************************/

void RNonHomogeneousTreeLikelihood::setData(const SiteContainer& sites) throw (Exception)
{
  if (data_) delete data_;
  data_ = PatternTools::getSequenceSubset(sites, *tree_->getRootNode());
  if (verbose_) ApplicationTools::displayTask("Initializing data structure");
  likelihoodData_->initLikelihoods(*data_, *modelSet_->getModel(0)); //We assume here that all models have the same number of states, and that they have the same 'init' method,
                                                                     //Which is a reasonable assumption as long as they share the same alphabet.
  if (verbose_) ApplicationTools::displayTaskDone();

  nbSites_         = likelihoodData_->getNumberOfSites();
  nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();
  nbStates_        = likelihoodData_->getNumberOfStates();

  if (verbose_) ApplicationTools::displayResult("Number of distinct sites",
                                                TextTools::toString(nbDistinctSites_));
  initialized_ = false;
}

/******************************************************************************/

double RNonHomogeneousTreeLikelihood::getLikelihood() const
{
  double l = 1.;
  for (size_t i = 0; i < nbSites_; i++)
  {
    l *= getLikelihoodForASite(i);
  }
  return l;
}

/******************************************************************************/

double RNonHomogeneousTreeLikelihood::getLogLikelihood() const
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

double RNonHomogeneousTreeLikelihood::getLikelihoodForASite(size_t site) const
{
  double l = 0;
  for (size_t i = 0; i < nbClasses_; i++)
  {
    l += getLikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
  }
  return l;
}

/******************************************************************************/

double RNonHomogeneousTreeLikelihood::getLogLikelihoodForASite(size_t site) const
{
  double l = 0;
  for (size_t i = 0; i < nbClasses_; i++)
  {
    l += getLikelihoodForASiteForARateClass(site, i) * rateDistribution_->getProbability(i);
  }
  //if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
  if (l < 0) l = 0; //May happen because of numerical errors.
  return log(l);
}

/******************************************************************************/

double RNonHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  double l = 0;
  Vdouble* la = &likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
  for (size_t i = 0; i < nbStates_; i++)
  {
    l += (*la)[i] * rootFreqs_[i];
  }
  return l;
}

/******************************************************************************/

double RNonHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  double l = 0;
  Vdouble* la = &likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass];
  for (size_t i = 0; i < nbStates_; i++)
  {
    l += (*la)[i] * rootFreqs_[i];
  }
  return log(l);
}

/******************************************************************************/

double RNonHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  return likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass][static_cast<size_t>(state)];
}

/******************************************************************************/

double RNonHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  return log(likelihoodData_->getLikelihoodArray(tree_->getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][rateClass][static_cast<size_t>(state)]);
}

/******************************************************************************/

void RNonHomogeneousTreeLikelihood::setParameters(const ParameterList& parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

void RNonHomogeneousTreeLikelihood::fireParameterChanged(const ParameterList& params)
{
  applyParameters();

  if (params.getCommonParametersWith(rateDistribution_->getIndependentParameters()).size() > 0)
  {
    computeAllTransitionProbabilities();
  }
  else
  {
    vector<int> ids;
    vector<string> tmp = params.getCommonParametersWith(modelSet_->getNodeParameters()).getParameterNames();
    for (size_t i = 0; i < tmp.size(); i++)
    {
      vector<int> tmpv = modelSet_->getNodesWithParameter(tmp[i]);
      ids = VectorTools::vectorUnion(ids, tmpv);
    }
    tmp = params.getCommonParametersWith(brLenParameters_).getParameterNames();
    vector<const Node*> nodes;
    for (size_t i = 0; i < ids.size(); i++)
    {
      nodes.push_back(idToNode_[ids[i]]);
    }
    vector<const Node*> tmpv;
    bool test = false;
    for (size_t i = 0; i < tmp.size(); i++)
    {
      if (tmp[i] == "BrLenRoot" || tmp[i] == "RootPosition")
      {
        if (!test)
        {
          tmpv.push_back(tree_->getRootNode()->getSon(0));
          tmpv.push_back(tree_->getRootNode()->getSon(1));
          test = true; //Add only once.
        }
      }
      else
        tmpv.push_back(nodes_[TextTools::to < size_t > (tmp[i].substr(5))]);
    }
    nodes = VectorTools::vectorUnion(nodes, tmpv);

    for (size_t i = 0; i < nodes.size(); i++)
    {
      computeTransitionProbabilitiesForNode(nodes[i]);
    }
    rootFreqs_ = modelSet_->getRootFrequencies();
  }
  computeTreeLikelihood();

  minusLogLik_ = -getLogLikelihood();
}

/******************************************************************************/

double RNonHomogeneousTreeLikelihood::getValue() const
throw (Exception)
{
  if (!isInitialized()) throw Exception("RNonHomogeneousTreeLikelihood::getValue(). Instance is not initialized.");
  return minusLogLik_;
}

/******************************************************************************
*                           First Order Derivatives                          *
******************************************************************************/

double RNonHomogeneousTreeLikelihood::getDLikelihoodForASiteForARateClass(
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

double RNonHomogeneousTreeLikelihood::getDLikelihoodForASite(size_t site) const
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

double RNonHomogeneousTreeLikelihood::getDLogLikelihoodForASite(size_t site) const
{
  // d(f(g(x)))/dx = dg(x)/dx . df(g(x))/dg :
  return getDLikelihoodForASite(site) / getLikelihoodForASite(site);
}

/******************************************************************************/

double RNonHomogeneousTreeLikelihood::getDLogLikelihood() const
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

double RNonHomogeneousTreeLikelihood::getFirstOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("RNonHomogeneousTreeLikelihood::getFirstOrderDerivative().", variable);
  if (getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameter are not implemented.");
  }
  if (getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }

  const_cast<RNonHomogeneousTreeLikelihood*>(this)->computeTreeDLikelihood(variable);
  return -getDLogLikelihood();
}

/******************************************************************************/

void RNonHomogeneousTreeLikelihood::computeTreeDLikelihood(const string& variable)
{
  if (variable == "BrLenRoot")
  {
    const Node* father = tree_->getRootNode();

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

      if (son->getId() == root1_)
      {
        const Node* root1 = father->getSon(0);
        const Node* root2 = father->getSon(1);
        vector<size_t> * _patternLinks_fatherroot1_ = &likelihoodData_->getArrayPositions(father->getId(), root1->getId());
        vector<size_t> * _patternLinks_fatherroot2_ = &likelihoodData_->getArrayPositions(father->getId(), root2->getId());
        VVVdouble* _likelihoodsroot1_ = &likelihoodData_->getLikelihoodArray(root1->getId());
        VVVdouble* _likelihoodsroot2_ = &likelihoodData_->getLikelihoodArray(root2->getId());
        double pos = getParameterValue("RootPosition");

        VVVdouble* dpxy_root1_  = &dpxy_[root1_];
        VVVdouble* dpxy_root2_  = &dpxy_[root2_];
        VVVdouble* pxy_root1_   = &pxy_[root1_];
        VVVdouble* pxy_root2_   = &pxy_[root2_];
        for (size_t i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoodsroot1__i = &(*_likelihoodsroot1_)[(*_patternLinks_fatherroot1_)[i]];
          VVdouble* _likelihoodsroot2__i = &(*_likelihoodsroot2_)[(*_patternLinks_fatherroot2_)[i]];
          VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
          for (size_t c = 0; c < nbClasses_; c++)
          {
            Vdouble* _likelihoodsroot1__i_c = &(*_likelihoodsroot1__i)[c];
            Vdouble* _likelihoodsroot2__i_c = &(*_likelihoodsroot2__i)[c];
            Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
            VVdouble* dpxy_root1__c  = &(*dpxy_root1_)[c];
            VVdouble* dpxy_root2__c  = &(*dpxy_root2_)[c];
            VVdouble* pxy_root1__c   = &(*pxy_root1_)[c];
            VVdouble* pxy_root2__c   = &(*pxy_root2_)[c];
            for (size_t x = 0; x < nbStates_; x++)
            {
              Vdouble* dpxy_root1__c_x  = &(*dpxy_root1__c)[x];
              Vdouble* dpxy_root2__c_x  = &(*dpxy_root2__c)[x];
              Vdouble* pxy_root1__c_x   = &(*pxy_root1__c)[x];
              Vdouble* pxy_root2__c_x   = &(*pxy_root2__c)[x];
              double dl1 = 0, dl2 = 0, l1 = 0, l2 = 0;
              for (size_t y = 0; y < nbStates_; y++)
              {
                dl1  += (*dpxy_root1__c_x)[y]  * (*_likelihoodsroot1__i_c)[y];
                dl2  += (*dpxy_root2__c_x)[y]  * (*_likelihoodsroot2__i_c)[y];
                l1   += (*pxy_root1__c_x)[y]   * (*_likelihoodsroot1__i_c)[y];
                l2   += (*pxy_root2__c_x)[y]   * (*_likelihoodsroot2__i_c)[y];
              }
              double dl = pos * dl1 * l2 + (1. - pos) * dl2 * l1;
              (*_dLikelihoods_father_i_c)[x] *= dl;
            }
          }
        }
      }
      else if (son->getId() == root2_)
      {
        //Do nothing, this was accounted when dealing with root1_
      }
      else
      {
        //Account for a putative multifurcation:
        vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
        VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

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
    return;
  }
  else if (variable == "RootPosition")
  {
    const Node* father = tree_->getRootNode();

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

      if (son->getId() == root1_)
      {
        const Node* root1 = father->getSon(0);
        const Node* root2 = father->getSon(1);
        vector<size_t> * _patternLinks_fatherroot1_ = &likelihoodData_->getArrayPositions(father->getId(), root1->getId());
        vector<size_t> * _patternLinks_fatherroot2_ = &likelihoodData_->getArrayPositions(father->getId(), root2->getId());
        VVVdouble* _likelihoodsroot1_ = &likelihoodData_->getLikelihoodArray(root1->getId());
        VVVdouble* _likelihoodsroot2_ = &likelihoodData_->getLikelihoodArray(root2->getId());
        double len = getParameterValue("BrLenRoot");

        VVVdouble* dpxy_root1_  = &dpxy_[root1_];
        VVVdouble* dpxy_root2_  = &dpxy_[root2_];
        VVVdouble* pxy_root1_   = &pxy_[root1_];
        VVVdouble* pxy_root2_   = &pxy_[root2_];
        for (size_t i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoodsroot1__i = &(*_likelihoodsroot1_)[(*_patternLinks_fatherroot1_)[i]];
          VVdouble* _likelihoodsroot2__i = &(*_likelihoodsroot2_)[(*_patternLinks_fatherroot2_)[i]];
          VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
          for (size_t c = 0; c < nbClasses_; c++)
          {
            Vdouble* _likelihoodsroot1__i_c = &(*_likelihoodsroot1__i)[c];
            Vdouble* _likelihoodsroot2__i_c = &(*_likelihoodsroot2__i)[c];
            Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
            VVdouble* dpxy_root1__c  = &(*dpxy_root1_)[c];
            VVdouble* dpxy_root2__c  = &(*dpxy_root2_)[c];
            VVdouble* pxy_root1__c   = &(*pxy_root1_)[c];
            VVdouble* pxy_root2__c   = &(*pxy_root2_)[c];
            for (size_t x = 0; x < nbStates_; x++)
            {
              Vdouble* dpxy_root1__c_x  = &(*dpxy_root1__c)[x];
              Vdouble* dpxy_root2__c_x  = &(*dpxy_root2__c)[x];
              Vdouble* pxy_root1__c_x   = &(*pxy_root1__c)[x];
              Vdouble* pxy_root2__c_x   = &(*pxy_root2__c)[x];
              double dl1 = 0, dl2 = 0, l1 = 0, l2 = 0;
              for (size_t y = 0; y < nbStates_; y++)
              {
                dl1  += (*dpxy_root1__c_x)[y]  * (*_likelihoodsroot1__i_c)[y];
                dl2  += (*dpxy_root2__c_x)[y]  * (*_likelihoodsroot2__i_c)[y];
                l1   += (*pxy_root1__c_x)[y]   * (*_likelihoodsroot1__i_c)[y];
                l2   += (*pxy_root2__c_x)[y]   * (*_likelihoodsroot2__i_c)[y];
              }
              double dl = len * (dl1 * l2 - dl2 * l1);
              (*_dLikelihoods_father_i_c)[x] *= dl;
            }
          }
        }
      }
      else if (son->getId() == root2_)
      {
        //Do nothing, this was accounted when dealing with root1_
      }
      else
      {
        //Account for a putative multifurcation:
        vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
        VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

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
    return;
  }

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

void RNonHomogeneousTreeLikelihood::computeDownSubtreeDLikelihood(const Node* node)
{
  const Node* father = node->getFather();
  // We assume that the _dLikelihoods array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if (father == 0) return; // We reached the root!

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

double RNonHomogeneousTreeLikelihood::getD2LikelihoodForASiteForARateClass(
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

double RNonHomogeneousTreeLikelihood::getD2LikelihoodForASite(size_t site) const
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

double RNonHomogeneousTreeLikelihood::getD2LogLikelihoodForASite(size_t site) const
{
  return getD2LikelihoodForASite(site) / getLikelihoodForASite(site)
         - pow( getDLikelihoodForASite(site) / getLikelihoodForASite(site), 2);
}

/******************************************************************************/

double RNonHomogeneousTreeLikelihood::getD2LogLikelihood() const
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

double RNonHomogeneousTreeLikelihood::getSecondOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("RNonHomogeneousTreeLikelihood::getSecondOrderDerivative().", variable);
  if (getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameter are not implemented.");
  }
  if (getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }

  const_cast<RNonHomogeneousTreeLikelihood*>(this)->computeTreeD2Likelihood(variable);
  return -getD2LogLikelihood();
}

/******************************************************************************/

void RNonHomogeneousTreeLikelihood::computeTreeD2Likelihood(const string& variable)
{
  if (variable == "BrLenRoot")
  {
    const Node* father = tree_->getRootNode();

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

      if (son->getId() == root1_)
      {
        const Node* root1 = father->getSon(0);
        const Node* root2 = father->getSon(1);
        vector<size_t> * _patternLinks_fatherroot1_ = &likelihoodData_->getArrayPositions(father->getId(), root1->getId());
        vector<size_t> * _patternLinks_fatherroot2_ = &likelihoodData_->getArrayPositions(father->getId(), root2->getId());
        VVVdouble* _likelihoodsroot1_ = &likelihoodData_->getLikelihoodArray(root1->getId());
        VVVdouble* _likelihoodsroot2_ = &likelihoodData_->getLikelihoodArray(root2->getId());
        double pos = getParameterValue("RootPosition");

        VVVdouble* d2pxy_root1_ = &d2pxy_[root1_];
        VVVdouble* d2pxy_root2_ = &d2pxy_[root2_];
        VVVdouble* dpxy_root1_  = &dpxy_[root1_];
        VVVdouble* dpxy_root2_  = &dpxy_[root2_];
        VVVdouble* pxy_root1_   = &pxy_[root1_];
        VVVdouble* pxy_root2_   = &pxy_[root2_];
        for (size_t i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoodsroot1__i = &(*_likelihoodsroot1_)[(*_patternLinks_fatherroot1_)[i]];
          VVdouble* _likelihoodsroot2__i = &(*_likelihoodsroot2_)[(*_patternLinks_fatherroot2_)[i]];
          VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
          for (size_t c = 0; c < nbClasses_; c++)
          {
            Vdouble* _likelihoodsroot1__i_c = &(*_likelihoodsroot1__i)[c];
            Vdouble* _likelihoodsroot2__i_c = &(*_likelihoodsroot2__i)[c];
            Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
            VVdouble* d2pxy_root1__c = &(*d2pxy_root1_)[c];
            VVdouble* d2pxy_root2__c = &(*d2pxy_root2_)[c];
            VVdouble* dpxy_root1__c  = &(*dpxy_root1_)[c];
            VVdouble* dpxy_root2__c  = &(*dpxy_root2_)[c];
            VVdouble* pxy_root1__c   = &(*pxy_root1_)[c];
            VVdouble* pxy_root2__c   = &(*pxy_root2_)[c];
            for (size_t x = 0; x < nbStates_; x++)
            {
              Vdouble* d2pxy_root1__c_x = &(*d2pxy_root1__c)[x];
              Vdouble* d2pxy_root2__c_x = &(*d2pxy_root2__c)[x];
              Vdouble* dpxy_root1__c_x  = &(*dpxy_root1__c)[x];
              Vdouble* dpxy_root2__c_x  = &(*dpxy_root2__c)[x];
              Vdouble* pxy_root1__c_x   = &(*pxy_root1__c)[x];
              Vdouble* pxy_root2__c_x   = &(*pxy_root2__c)[x];
              double d2l1 = 0, d2l2 = 0, dl1 = 0, dl2 = 0, l1 = 0, l2 = 0;
              for (size_t y = 0; y < nbStates_; y++)
              {
                d2l1 += (*d2pxy_root1__c_x)[y] * (*_likelihoodsroot1__i_c)[y];
                d2l2 += (*d2pxy_root2__c_x)[y] * (*_likelihoodsroot2__i_c)[y];
                dl1  += (*dpxy_root1__c_x)[y]  * (*_likelihoodsroot1__i_c)[y];
                dl2  += (*dpxy_root2__c_x)[y]  * (*_likelihoodsroot2__i_c)[y];
                l1   += (*pxy_root1__c_x)[y]   * (*_likelihoodsroot1__i_c)[y];
                l2   += (*pxy_root2__c_x)[y]   * (*_likelihoodsroot2__i_c)[y];
              }
              double d2l = pos * pos * d2l1 * l2 + (1. - pos) * (1. - pos) * d2l2 * l1 + 2 * pos * (1. - pos) * dl1 * dl2;
              (*_d2Likelihoods_father_i_c)[x] *= d2l;
            }
          }
        }
      }
      else if (son->getId() == root2_)
      {
        //Do nothing, this was accounted when dealing with root1_
      }
      else
      {
        //Account for a putative multifurcation:
        vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
        VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

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
    return;
  }
  else if (variable == "RootPosition")
  {
    const Node* father = tree_->getRootNode();

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

      if (son->getId() == root1_)
      {
        const Node* root1 = father->getSon(0);
        const Node* root2 = father->getSon(1);
        vector<size_t> * _patternLinks_fatherroot1_ = &likelihoodData_->getArrayPositions(father->getId(), root1->getId());
        vector<size_t> * _patternLinks_fatherroot2_ = &likelihoodData_->getArrayPositions(father->getId(), root2->getId());
        VVVdouble* _likelihoodsroot1_ = &likelihoodData_->getLikelihoodArray(root1->getId());
        VVVdouble* _likelihoodsroot2_ = &likelihoodData_->getLikelihoodArray(root2->getId());
        double len = getParameterValue("BrLenRoot");

        VVVdouble* d2pxy_root1_ = &d2pxy_[root1_];
        VVVdouble* d2pxy_root2_ = &d2pxy_[root2_];
        VVVdouble* dpxy_root1_  = &dpxy_[root1_];
        VVVdouble* dpxy_root2_  = &dpxy_[root2_];
        VVVdouble* pxy_root1_   = &pxy_[root1_];
        VVVdouble* pxy_root2_   = &pxy_[root2_];
        for (size_t i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoodsroot1__i = &(*_likelihoodsroot1_)[(*_patternLinks_fatherroot1_)[i]];
          VVdouble* _likelihoodsroot2__i = &(*_likelihoodsroot2_)[(*_patternLinks_fatherroot2_)[i]];
          VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
          for (size_t c = 0; c < nbClasses_; c++)
          {
            Vdouble* _likelihoodsroot1__i_c = &(*_likelihoodsroot1__i)[c];
            Vdouble* _likelihoodsroot2__i_c = &(*_likelihoodsroot2__i)[c];
            Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
            VVdouble* d2pxy_root1__c = &(*d2pxy_root1_)[c];
            VVdouble* d2pxy_root2__c = &(*d2pxy_root2_)[c];
            VVdouble* dpxy_root1__c  = &(*dpxy_root1_)[c];
            VVdouble* dpxy_root2__c  = &(*dpxy_root2_)[c];
            VVdouble* pxy_root1__c   = &(*pxy_root1_)[c];
            VVdouble* pxy_root2__c   = &(*pxy_root2_)[c];
            for (size_t x = 0; x < nbStates_; x++)
            {
              Vdouble* d2pxy_root1__c_x = &(*d2pxy_root1__c)[x];
              Vdouble* d2pxy_root2__c_x = &(*d2pxy_root2__c)[x];
              Vdouble* dpxy_root1__c_x  = &(*dpxy_root1__c)[x];
              Vdouble* dpxy_root2__c_x  = &(*dpxy_root2__c)[x];
              Vdouble* pxy_root1__c_x   = &(*pxy_root1__c)[x];
              Vdouble* pxy_root2__c_x   = &(*pxy_root2__c)[x];
              double d2l1 = 0, d2l2 = 0, dl1 = 0, dl2 = 0, l1 = 0, l2 = 0;
              for (size_t y = 0; y < nbStates_; y++)
              {
                d2l1 += (*d2pxy_root1__c_x)[y] * (*_likelihoodsroot1__i_c)[y];
                d2l2 += (*d2pxy_root2__c_x)[y] * (*_likelihoodsroot2__i_c)[y];
                dl1  += (*dpxy_root1__c_x)[y]  * (*_likelihoodsroot1__i_c)[y];
                dl2  += (*dpxy_root2__c_x)[y]  * (*_likelihoodsroot2__i_c)[y];
                l1   += (*pxy_root1__c_x)[y]   * (*_likelihoodsroot1__i_c)[y];
                l2   += (*pxy_root2__c_x)[y]   * (*_likelihoodsroot2__i_c)[y];
              }
              double d2l = len * len * (d2l1 * l2 + d2l2 * l1 - 2 * dl1 * dl2);
              (*_d2Likelihoods_father_i_c)[x] *= d2l;
            }
          }
        }
      }
      else if (son->getId() == root2_)
      {
        //Do nothing, this was accounted when dealing with root1_
      }
      else
      {
        //Account for a putative multifurcation:
        vector<size_t> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
        VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

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
    return;
  }

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

void RNonHomogeneousTreeLikelihood::computeDownSubtreeD2Likelihood(const Node* node)
{
  const Node* father = node->getFather();
  // We assume that the _dLikelihoods array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if (father == 0) return; // We reached the root!

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

void RNonHomogeneousTreeLikelihood::computeTreeLikelihood()
{
  computeSubtreeLikelihood(tree_->getRootNode());
}

/******************************************************************************/

void RNonHomogeneousTreeLikelihood::computeSubtreeLikelihood(const Node* node)
{
  if (node->isLeaf()) return;

  size_t nbSites  = likelihoodData_->getLikelihoodArray(node->getId()).size();
  size_t nbNodes  = node->getNumberOfSons();

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
          {
            likelihood += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
          }
          (*_likelihoods_node_i_c)[x] *= likelihood;
        }
      }
    }
  }
}


/******************************************************************************/

void RNonHomogeneousTreeLikelihood::displayLikelihood(const Node* node)
{
  cout << "Likelihoods at node " << node->getName() << ": " << endl;
  displayLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId()));
  cout << "                                         ***" << endl;
}

