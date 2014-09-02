//
// File: DRNonHomogeneousTreeLikelihood.cpp
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

#include "DRNonHomogeneousTreeLikelihood.h"
#include "../PatternTools.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

DRNonHomogeneousTreeLikelihood::DRNonHomogeneousTreeLikelihood(
  const Tree& tree,
  SubstitutionModelSet* modelSet,
  DiscreteDistribution* rDist,
  bool verbose,
  bool reparametrizeRoot)
throw (Exception) :
  AbstractNonHomogeneousTreeLikelihood(tree, modelSet, rDist, verbose, reparametrizeRoot),
  likelihoodData_(0),
  minusLogLik_(-1.)
{
  if (!modelSet->isFullySetUpFor(tree))
    throw Exception("DRNonHomogeneousTreeLikelihood(constructor). Model set is not fully specified.");
  init_();
}

/******************************************************************************/

DRNonHomogeneousTreeLikelihood::DRNonHomogeneousTreeLikelihood(
  const Tree& tree,
  const SiteContainer& data,
  SubstitutionModelSet* modelSet,
  DiscreteDistribution* rDist,
  bool verbose,
  bool reparametrizeRoot)
throw (Exception) :
  AbstractNonHomogeneousTreeLikelihood(tree, modelSet, rDist, verbose, reparametrizeRoot),
  likelihoodData_(0),
  minusLogLik_(-1.)
{
  if (!modelSet->isFullySetUpFor(tree))
    throw Exception("DRNonHomogeneousTreeLikelihood(constructor). Model set is not fully specified.");
  init_();
  setData(data);
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::init_() throw (Exception)
{
  likelihoodData_ = new DRASDRTreeLikelihoodData(
    tree_,
    rateDistribution_->getNumberOfCategories());
}

/******************************************************************************/

DRNonHomogeneousTreeLikelihood::DRNonHomogeneousTreeLikelihood(const DRNonHomogeneousTreeLikelihood& lik) :
  AbstractNonHomogeneousTreeLikelihood(lik),
  likelihoodData_(0),
  minusLogLik_(lik.minusLogLik_)
{
  likelihoodData_ = dynamic_cast<DRASDRTreeLikelihoodData*>(lik.likelihoodData_->clone());
  likelihoodData_->setTree(tree_);
}

/******************************************************************************/

DRNonHomogeneousTreeLikelihood& DRNonHomogeneousTreeLikelihood::operator=(const DRNonHomogeneousTreeLikelihood& lik)
{
  AbstractNonHomogeneousTreeLikelihood::operator=(lik);
  if (likelihoodData_)
    delete likelihoodData_;
  likelihoodData_ = dynamic_cast<DRASDRTreeLikelihoodData*>(lik.likelihoodData_->clone());
  likelihoodData_->setTree(tree_);
  minusLogLik_ = lik.minusLogLik_;
  return *this;
}

/******************************************************************************/

DRNonHomogeneousTreeLikelihood::~DRNonHomogeneousTreeLikelihood()
{
  delete likelihoodData_;
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::setData(const SiteContainer& sites) throw (Exception)
{
  if (data_)
    delete data_;
  data_ = PatternTools::getSequenceSubset(sites, *tree_->getRootNode());
  if (verbose_)
    ApplicationTools::displayTask("Initializing data structure");
  likelihoodData_->initLikelihoods(*data_, *modelSet_->getModel(0)); // We assume here that all models have the same number of states, and that they have the same 'init' method,
                                                                     // Which is a reasonable assumption as long as they share the same alphabet.
  if (verbose_)
    ApplicationTools::displayTaskDone();

  nbSites_         = likelihoodData_->getNumberOfSites();
  nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();
  nbStates_        = likelihoodData_->getNumberOfStates();

  if (verbose_)
    ApplicationTools::displayResult("Number of distinct sites",
                                    TextTools::toString(nbDistinctSites_));
  initialized_ = false;
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getLikelihood() const
{
  double l = 1.;
  Vdouble* lik = &likelihoodData_->getRootRateSiteLikelihoodArray();
  const vector<unsigned int>* w = &likelihoodData_->getWeights();
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    l *= std::pow((*lik)[i], (int)(*w)[i]);
  }
  return l;
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getLogLikelihood() const
{
  double ll = 0;
  Vdouble* lik = &likelihoodData_->getRootRateSiteLikelihoodArray();
  const vector<unsigned int>* w = &likelihoodData_->getWeights();
  vector<double> la(nbDistinctSites_);
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    la[i] = (*w)[i] * log((*lik)[i]);
  }
  sort(la.begin(), la.end());
  for (size_t i = nbDistinctSites_; i > 0; i--)
  {
    ll += la[i - 1];
  }
  return ll;
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getLikelihoodForASite(size_t site) const
{
  return likelihoodData_->getRootRateSiteLikelihoodArray()[likelihoodData_->getRootArrayPosition(site)];
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getLogLikelihoodForASite(size_t site) const
{
  return log(likelihoodData_->getRootRateSiteLikelihoodArray()[likelihoodData_->getRootArrayPosition(site)]);
}

/******************************************************************************/
double DRNonHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  return likelihoodData_->getRootSiteLikelihoodArray()[likelihoodData_->getRootArrayPosition(site)][rateClass];
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  return log(likelihoodData_->getRootSiteLikelihoodArray()[likelihoodData_->getRootArrayPosition(site)][rateClass]);
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  return likelihoodData_->getRootLikelihoodArray()[likelihoodData_->getRootArrayPosition(site)][rateClass][static_cast<size_t>(state)];
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  return log(likelihoodData_->getRootLikelihoodArray()[likelihoodData_->getRootArrayPosition(site)][rateClass][static_cast<size_t>(state)]);
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::setParameters(const ParameterList& parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::fireParameterChanged(const ParameterList& params)
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
          test = true; // Add only once.
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
  if (computeFirstOrderDerivatives_)
  {
    computeTreeDLikelihoods();
  }
  if (computeSecondOrderDerivatives_)
  {
    computeTreeD2Likelihoods();
  }
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getValue() const
throw (Exception)
{
  if (!isInitialized())
    throw Exception("DRNonHomogeneousTreeLikelihood::getValue(). Instance is not initialized.");
  return -getLogLikelihood();
}

/******************************************************************************
*                           First Order Derivatives                          *
******************************************************************************/
void DRNonHomogeneousTreeLikelihood::computeTreeDLikelihoodAtNode(const Node* node)
{
  const Node* father = node->getFather();
  VVVdouble* _likelihoods_father_node = &likelihoodData_->getLikelihoodArray(father->getId(), node->getId());
  Vdouble* _dLikelihoods_node = &likelihoodData_->getDLikelihoodArray(node->getId());
  VVVdouble*  pxy__node = &pxy_[node->getId()];
  VVVdouble* dpxy__node = &dpxy_[node->getId()];
  VVVdouble larray;
  computeLikelihoodAtNode_(father, larray);
  Vdouble* rootLikelihoodsSR = &likelihoodData_->getRootRateSiteLikelihoodArray();

  double dLi, dLic, dLicx, numerator, denominator;
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    VVdouble* _likelihoods_father_node_i = &(*_likelihoods_father_node)[i];
    VVdouble* larray_i = &larray[i];
    dLi = 0;
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* _likelihoods_father_node_i_c = &(*_likelihoods_father_node_i)[c];
      Vdouble* larray_i_c = &(*larray_i)[c];
      VVdouble*  pxy__node_c = &(*pxy__node)[c];
      VVdouble* dpxy__node_c = &(*dpxy__node)[c];
      dLic = 0;
      for (size_t x = 0; x < nbStates_; x++)
      {
        numerator = 0;
        denominator = 0;
        Vdouble*  pxy__node_c_x = &(*pxy__node_c)[x];
        Vdouble* dpxy__node_c_x = &(*dpxy__node_c)[x];
        dLicx = 0;
        for (size_t y = 0; y < nbStates_; y++)
        {
          numerator   += (*dpxy__node_c_x)[y] * (*_likelihoods_father_node_i_c)[y];
          denominator += (*pxy__node_c_x)[y] * (*_likelihoods_father_node_i_c)[y];
        }
        dLicx = denominator == 0. ? 0. : (*larray_i_c)[x] * numerator / denominator;
        dLic += dLicx;
      }
      dLi += rateDistribution_->getProbability(c) * dLic;
    }
    (*_dLikelihoods_node)[i] = dLi / (*rootLikelihoodsSR)[i];
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::computeTreeDLikelihoods()
{
  for (size_t k = 0; k < nbNodes_; k++)
  {
    computeTreeDLikelihoodAtNode(nodes_[k]);
  }
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getFirstOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("DRNonHomogeneousTreeLikelihood::getFirstOrderDerivative().", variable);
  if (getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameters are not implemented.");
  }
  if (getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }

  //
  // Computation for branch lengths:
  //
  const vector<unsigned int>* w = &likelihoodData_->getWeights();
  Vdouble* _dLikelihoods_branch;
  if (variable == "BrLenRoot")
  {
    _dLikelihoods_branch = &likelihoodData_->getDLikelihoodArray(root1_);
    double d1 = 0;
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      d1 -= (*w)[i] * (*_dLikelihoods_branch)[i];
    }
    _dLikelihoods_branch = &likelihoodData_->getDLikelihoodArray(root2_);
    double d2 = 0;
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      d2 -= (*w)[i] * (*_dLikelihoods_branch)[i];
    }
    double pos = getParameterValue("RootPosition");
    return pos * d1 + (1. - pos) * d2;
  }
  else if (variable == "RootPosition")
  {
    _dLikelihoods_branch = &likelihoodData_->getDLikelihoodArray(root1_);
    double d1 = 0;
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      d1 -= (*w)[i] * (*_dLikelihoods_branch)[i];
    }
    _dLikelihoods_branch = &likelihoodData_->getDLikelihoodArray(root2_);
    double d2 = 0;
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      d2 -= (*w)[i] * (*_dLikelihoods_branch)[i];
    }
    double len = getParameterValue("BrLenRoot");
    return len * (d1 - d2);
  }
  else
  {
    // Get the node with the branch whose length must be derivated:
    size_t brI = TextTools::to<size_t>(variable.substr(5));
    const Node* branch = nodes_[brI];
    _dLikelihoods_branch = &likelihoodData_->getDLikelihoodArray(branch->getId());
    double d = 0;
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      d += (*w)[i] * (*_dLikelihoods_branch)[i];
    }
    return -d;
  }
}

/******************************************************************************
*                           Second Order Derivatives                         *
******************************************************************************/
void DRNonHomogeneousTreeLikelihood::computeTreeD2LikelihoodAtNode(const Node* node)
{
  const Node* father = node->getFather();
  VVVdouble* _likelihoods_father_node = &likelihoodData_->getLikelihoodArray(father->getId(), node->getId());
  Vdouble* _d2Likelihoods_node = &likelihoodData_->getD2LikelihoodArray(node->getId());
  VVVdouble*   pxy__node = &pxy_[node->getId()];
  VVVdouble* d2pxy__node = &d2pxy_[node->getId()];
  VVVdouble larray;
  computeLikelihoodAtNode_(father, larray);
  Vdouble* rootLikelihoodsSR = &likelihoodData_->getRootRateSiteLikelihoodArray();

  double d2Li, d2Lic, d2Licx, numerator, denominator;

  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    VVdouble* _likelihoods_father_node_i = &(*_likelihoods_father_node)[i];
    VVdouble* larray_i = &larray[i];
    d2Li = 0;
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* _likelihoods_father_node_i_c = &(*_likelihoods_father_node_i)[c];
      Vdouble* larray_i_c = &(*larray_i)[c];
      VVdouble*   pxy__node_c = &(*pxy__node)[c];
      VVdouble* d2pxy__node_c = &(*d2pxy__node)[c];
      d2Lic = 0;
      for (size_t x = 0; x < nbStates_; x++)
      {
        numerator = 0;
        denominator = 0;
        Vdouble*   pxy__node_c_x = &(*pxy__node_c)[x];
        Vdouble* d2pxy__node_c_x = &(*d2pxy__node_c)[x];
        d2Licx = 0;
        for (size_t y = 0; y < nbStates_; y++)
        {
          numerator   += (*d2pxy__node_c_x)[y] * (*_likelihoods_father_node_i_c)[y];
          denominator += (*pxy__node_c_x)[y] * (*_likelihoods_father_node_i_c)[y];
        }
        d2Licx = denominator == 0. ? 0. : (*larray_i_c)[x] * numerator / denominator;
        d2Lic += d2Licx;
      }
      d2Li += rateDistribution_->getProbability(c) * d2Lic;
    }
    (*_d2Likelihoods_node)[i] = d2Li / (*rootLikelihoodsSR)[i];
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::computeTreeD2Likelihoods()
{
  for (size_t k = 0; k < nbNodes_; k++)
  {
    computeTreeD2LikelihoodAtNode(nodes_[k]);
  }
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getSecondOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("DRNonHomogeneousTreeLikelihood::getSecondOrderDerivative().", variable);
  if (getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameters are not implemented.");
  }
  if (getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }

  //
  // Computation for branch lengths:
  //

  const vector<unsigned int>* w = &likelihoodData_->getWeights();
  // We can't deduce second order derivatives regarding BrLenRoot and RootPosition from the
  // branch length derivatives. We need a bit more calculations...
  // NB: we could save a few calculations here...
  if (variable == "BrLenRoot")
  {
    const Node* father = tree_->getRootNode();

    // Compute dLikelihoods array for the father node.
    // Fist initialize to 1:
    VVVdouble dLikelihoods_father(nbDistinctSites_);
    VVVdouble d2Likelihoods_father(nbDistinctSites_);
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble* dLikelihoods_father_i = &dLikelihoods_father[i];
      VVdouble* d2Likelihoods_father_i = &d2Likelihoods_father[i];
      dLikelihoods_father_i->resize(nbClasses_);
      d2Likelihoods_father_i->resize(nbClasses_);
      for (size_t c = 0; c < nbClasses_; c++)
      {
        Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
        Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
        dLikelihoods_father_i_c->resize(nbStates_);
        d2Likelihoods_father_i_c->resize(nbStates_);
        for (size_t s = 0; s < nbStates_; s++)
        {
          (*dLikelihoods_father_i_c)[s] = 1.;
          (*d2Likelihoods_father_i_c)[s] = 1.;
        }
      }
    }

    size_t nbNodes = father->getNumberOfSons();
    for (size_t l = 0; l < nbNodes; l++)
    {
      const Node* son = father->getSon(l);

      if (son->getId() == root1_)
      {
        VVVdouble* _likelihoodsroot1_ = &likelihoodData_->getLikelihoodArray(father->getId(), root1_);
        VVVdouble* _likelihoodsroot2_ = &likelihoodData_->getLikelihoodArray(father->getId(), root2_);
        double pos = getParameterValue("RootPosition");

        VVVdouble* d2pxy_root1_ = &d2pxy_[root1_];
        VVVdouble* d2pxy_root2_ = &d2pxy_[root2_];
        VVVdouble* dpxy_root1_  = &dpxy_[root1_];
        VVVdouble* dpxy_root2_  = &dpxy_[root2_];
        VVVdouble* pxy_root1_   = &pxy_[root1_];
        VVVdouble* pxy_root2_   = &pxy_[root2_];
        for (size_t i = 0; i < nbDistinctSites_; i++)
        {
          VVdouble* _likelihoodsroot1__i = &(*_likelihoodsroot1_)[i];
          VVdouble* _likelihoodsroot2__i = &(*_likelihoodsroot2_)[i];
          VVdouble* dLikelihoods_father_i = &dLikelihoods_father[i];
          VVdouble* d2Likelihoods_father_i = &d2Likelihoods_father[i];
          for (size_t c = 0; c < nbClasses_; c++)
          {
            Vdouble* _likelihoodsroot1__i_c = &(*_likelihoodsroot1__i)[c];
            Vdouble* _likelihoodsroot2__i_c = &(*_likelihoodsroot2__i)[c];
            Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
            Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
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
              double dl = pos * dl1 * l2 + (1. - pos) * dl2 * l1;
              double d2l = pos * pos * d2l1 * l2 + (1. - pos) * (1. - pos) * d2l2 * l1 + 2 * pos * (1. - pos) * dl1 * dl2;
              (*dLikelihoods_father_i_c)[x] *= dl;
              (*d2Likelihoods_father_i_c)[x] *= d2l;
            }
          }
        }
      }
      else if (son->getId() == root2_)
      {
        // Do nothing, this was accounted when dealing with root1_
      }
      else
      {
        // Account for a putative multifurcation:
        VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(father->getId(), son->getId());

        VVVdouble* pxy__son = &pxy_[son->getId()];
        for (size_t i = 0; i < nbDistinctSites_; i++)
        {
          VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[i];
          VVdouble* dLikelihoods_father_i = &dLikelihoods_father[i];
          VVdouble* d2Likelihoods_father_i = &d2Likelihoods_father[i];
          for (size_t c = 0; c < nbClasses_; c++)
          {
            Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
            Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
            Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
            VVdouble* pxy__son_c = &(*pxy__son)[c];
            for (size_t x = 0; x < nbStates_; x++)
            {
              double dl = 0;
              Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
              for (size_t y = 0; y < nbStates_; y++)
              {
                dl += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
              }
              (*dLikelihoods_father_i_c)[x] *= dl;
              (*d2Likelihoods_father_i_c)[x] *= dl;
            }
          }
        }
      }
    }
    Vdouble* rootLikelihoodsSR = &likelihoodData_->getRootRateSiteLikelihoodArray();
    double d2l = 0, dlx, d2lx;
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble* dLikelihoods_father_i = &dLikelihoods_father[i];
      VVdouble* d2Likelihoods_father_i = &d2Likelihoods_father[i];
      dlx = 0, d2lx = 0;
      for (size_t c = 0; c < nbClasses_; c++)
      {
        Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
        Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
        for (size_t x = 0; x < nbStates_; x++)
        {
          dlx += rateDistribution_->getProbability(c) * rootFreqs_[x] * (*dLikelihoods_father_i_c)[x];
          d2lx += rateDistribution_->getProbability(c) * rootFreqs_[x] * (*d2Likelihoods_father_i_c)[x];
        }
      }
      d2l += (*w)[i] * (d2lx / (*rootLikelihoodsSR)[i] - pow(dlx / (*rootLikelihoodsSR)[i], 2));
    }
    return -d2l;
  }
  else if (variable == "RootPosition")
  {
    const Node* father = tree_->getRootNode();

    // Compute dLikelihoods array for the father node.
    // Fist initialize to 1:
    VVVdouble dLikelihoods_father(nbDistinctSites_);
    VVVdouble d2Likelihoods_father(nbDistinctSites_);
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble* dLikelihoods_father_i = &dLikelihoods_father[i];
      VVdouble* d2Likelihoods_father_i = &d2Likelihoods_father[i];
      dLikelihoods_father_i->resize(nbClasses_);
      d2Likelihoods_father_i->resize(nbClasses_);
      for (size_t c = 0; c < nbClasses_; c++)
      {
        Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
        Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
        dLikelihoods_father_i_c->resize(nbStates_);
        d2Likelihoods_father_i_c->resize(nbStates_);
        for (size_t s = 0; s < nbStates_; s++)
        {
          (*dLikelihoods_father_i_c)[s] = 1.;
          (*d2Likelihoods_father_i_c)[s] = 1.;
        }
      }
    }

    size_t nbNodes = father->getNumberOfSons();
    for (size_t l = 0; l < nbNodes; l++)
    {
      const Node* son = father->getSon(l);

      if (son->getId() == root1_)
      {
        VVVdouble* _likelihoodsroot1_ = &likelihoodData_->getLikelihoodArray(father->getId(), root1_);
        VVVdouble* _likelihoodsroot2_ = &likelihoodData_->getLikelihoodArray(father->getId(), root2_);
        double len = getParameterValue("BrLenRoot");

        VVVdouble* d2pxy_root1_ = &d2pxy_[root1_];
        VVVdouble* d2pxy_root2_ = &d2pxy_[root2_];
        VVVdouble* dpxy_root1_  = &dpxy_[root1_];
        VVVdouble* dpxy_root2_  = &dpxy_[root2_];
        VVVdouble* pxy_root1_   = &pxy_[root1_];
        VVVdouble* pxy_root2_   = &pxy_[root2_];
        for (size_t i = 0; i < nbDistinctSites_; i++)
        {
          VVdouble* _likelihoodsroot1__i = &(*_likelihoodsroot1_)[i];
          VVdouble* _likelihoodsroot2__i = &(*_likelihoodsroot2_)[i];
          VVdouble* dLikelihoods_father_i = &dLikelihoods_father[i];
          VVdouble* d2Likelihoods_father_i = &d2Likelihoods_father[i];
          for (size_t c = 0; c < nbClasses_; c++)
          {
            Vdouble* _likelihoodsroot1__i_c = &(*_likelihoodsroot1__i)[c];
            Vdouble* _likelihoodsroot2__i_c = &(*_likelihoodsroot2__i)[c];
            Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
            Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
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
              double dl = len * (dl1 * l2 - dl2 * l1);
              double d2l = len * len * (d2l1 * l2 + d2l2 * l1 - 2 * dl1 * dl2);
              (*dLikelihoods_father_i_c)[x] *= dl;
              (*d2Likelihoods_father_i_c)[x] *= d2l;
            }
          }
        }
      }
      else if (son->getId() == root2_)
      {
        // Do nothing, this was accounted when dealing with root1_
      }
      else
      {
        // Account for a putative multifurcation:
        VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(father->getId(), son->getId());

        VVVdouble* pxy__son = &pxy_[son->getId()];
        for (size_t i = 0; i < nbDistinctSites_; i++)
        {
          VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[i];
          VVdouble* dLikelihoods_father_i = &dLikelihoods_father[i];
          VVdouble* d2Likelihoods_father_i = &d2Likelihoods_father[i];
          for (size_t c = 0; c < nbClasses_; c++)
          {
            Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
            Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
            Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
            VVdouble* pxy__son_c = &(*pxy__son)[c];
            for (size_t x = 0; x < nbStates_; x++)
            {
              double dl = 0;
              Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
              for (size_t y = 0; y < nbStates_; y++)
              {
                dl += (*pxy__son_c_x)[y] * (*_likelihoods_son_i_c)[y];
              }
              (*dLikelihoods_father_i_c)[x] *= dl;
              (*d2Likelihoods_father_i_c)[x] *= dl;
            }
          }
        }
      }
    }
    Vdouble* rootLikelihoodsSR = &likelihoodData_->getRootRateSiteLikelihoodArray();
    double d2l = 0, dlx, d2lx;
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble* dLikelihoods_father_i = &dLikelihoods_father[i];
      VVdouble* d2Likelihoods_father_i = &d2Likelihoods_father[i];
      dlx = 0, d2lx = 0;
      for (size_t c = 0; c < nbClasses_; c++)
      {
        Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
        Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
        for (size_t x = 0; x < nbStates_; x++)
        {
          dlx += rateDistribution_->getProbability(c) * rootFreqs_[x] * (*dLikelihoods_father_i_c)[x];
          d2lx += rateDistribution_->getProbability(c) * rootFreqs_[x] * (*d2Likelihoods_father_i_c)[x];
        }
      }
      d2l += (*w)[i] * (d2lx / (*rootLikelihoodsSR)[i] - pow(dlx / (*rootLikelihoodsSR)[i], 2));
    }
    return -d2l;
  }
  else
  {
    Vdouble* _dLikelihoods_branch;
    Vdouble* _d2Likelihoods_branch;
    // Get the node with the branch whose length must be derivated:
    size_t brI = TextTools::to<size_t>(variable.substr(5));
    const Node* branch = nodes_[brI];
    _dLikelihoods_branch = &likelihoodData_->getDLikelihoodArray(branch->getId());
    _d2Likelihoods_branch = &likelihoodData_->getD2LikelihoodArray(branch->getId());
    double d2l = 0;
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      d2l += (*w)[i] * ((*_d2Likelihoods_branch)[i] - pow((*_dLikelihoods_branch)[i], 2));
    }
    return -d2l;
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::resetLikelihoodArrays(const Node* node)
{
  for (size_t n = 0; n < node->getNumberOfSons(); n++)
  {
    const Node* subNode = node->getSon(n);
    resetLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId(), subNode->getId()));
  }
  if (node->hasFather())
  {
    const Node* father = node->getFather();
    resetLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId(), father->getId()));
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::computeTreeLikelihood()
{
  computeSubtreeLikelihoodPostfix(tree_->getRootNode());
  computeSubtreeLikelihoodPrefix(tree_->getRootNode());
  computeRootLikelihood();
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::computeSubtreeLikelihoodPostfix(const Node* node)
{
//  if(node->isLeaf()) return;
// cout << node->getId() << "\t" << (node->hasName()?node->getName():"") << endl;
  if (node->getNumberOfSons() == 0)
    return;

  // Set all likelihood arrays to 1 for a start:
  resetLikelihoodArrays(node);

  map<int, VVVdouble>* _likelihoods_node = &likelihoodData_->getLikelihoodArrays(node->getId());
  size_t nbNodes = node->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    // For each son node...

    const Node* son = node->getSon(l);
    VVVdouble* _likelihoods_node_son = &(*_likelihoods_node)[son->getId()];

    if (son->isLeaf())
    {
      VVdouble* _likelihoods_leaf = &likelihoodData_->getLeafLikelihoods(son->getId());
      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        // For each site in the sequence,
        Vdouble* _likelihoods_leaf_i = &(*_likelihoods_leaf)[i];
        VVdouble* _likelihoods_node_son_i = &(*_likelihoods_node_son)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          // For each rate classe,
          Vdouble* _likelihoods_node_son_i_c = &(*_likelihoods_node_son_i)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            // For each initial state,
            (*_likelihoods_node_son_i_c)[x] = (*_likelihoods_leaf_i)[x];
          }
        }
      }
    }
    else
    {
      computeSubtreeLikelihoodPostfix(son); // Recursive method:
      size_t nbSons = son->getNumberOfSons();
      map<int, VVVdouble>* _likelihoods_son = &likelihoodData_->getLikelihoodArrays(son->getId());

      vector<const VVVdouble*> iLik(nbSons);
      vector<const VVVdouble*> tProb(nbSons);
      for (size_t n = 0; n < nbSons; n++)
      {
        const Node* sonSon = son->getSon(n);
        tProb[n] = &pxy_[sonSon->getId()];
        iLik[n] = &(*_likelihoods_son)[sonSon->getId()];
      }
      computeLikelihoodFromArrays(iLik, tProb, *_likelihoods_node_son, nbSons, nbDistinctSites_, nbClasses_, nbStates_, false);
    }
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::computeSubtreeLikelihoodPrefix(const Node* node)
{
  if (!node->hasFather())
  {
    // 'node' is the root of the tree.
    // Just call the method on each son node:
    size_t nbSons = node->getNumberOfSons();
    for (size_t n = 0; n < nbSons; n++)
    {
      computeSubtreeLikelihoodPrefix(node->getSon(n));
    }
    return;
  }
  else
  {
    const Node* father = node->getFather();
    map<int, VVVdouble>* _likelihoods_node = &likelihoodData_->getLikelihoodArrays(node->getId());
    map<int, VVVdouble>* _likelihoods_father = &likelihoodData_->getLikelihoodArrays(father->getId());
    VVVdouble* _likelihoods_node_father = &(*_likelihoods_node)[father->getId()];
    if (node->isLeaf())
    {
      resetLikelihoodArray(*_likelihoods_node_father);
    }

    if (father->isLeaf())
    {
      // If the tree is rooted by a leaf
      VVdouble* _likelihoods_leaf = &likelihoodData_->getLeafLikelihoods(father->getId());
      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        // For each site in the sequence,
        Vdouble* _likelihoods_leaf_i = &(*_likelihoods_leaf)[i];
        VVdouble* _likelihoods_node_father_i = &(*_likelihoods_node_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          // For each rate classe,
          Vdouble* _likelihoods_node_father_i_c = &(*_likelihoods_node_father_i)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            // For each initial state,
            (*_likelihoods_node_father_i_c)[x] = (*_likelihoods_leaf_i)[x];
          }
        }
      }
    }
    else
    {
      vector<const Node*> nodes;
      // Add brothers:
      size_t nbFatherSons = father->getNumberOfSons();
      for (size_t n = 0; n < nbFatherSons; n++)
      {
        const Node* son = father->getSon(n);
        if (son->getId() != node->getId())
          nodes.push_back(son);  // This is a real brother, not current node!
      }
      // Now the real stuff... We've got to compute the likelihoods for the
      // subtree defined by node 'father'.
      // This is the same as postfix method, but with different subnodes.

      size_t nbSons = nodes.size(); // In case of a bifurcating tree this is equal to 1.

      vector<const VVVdouble*> iLik(nbSons);
      vector<const VVVdouble*> tProb(nbSons);
      for (size_t n = 0; n < nbSons; n++)
      {
        const Node* fatherSon = nodes[n];
        tProb[n] = &pxy_[fatherSon->getId()];
        iLik[n] = &(*_likelihoods_father)[fatherSon->getId()];
      }

      if (father->hasFather())
      {
        const Node* fatherFather = father->getFather();
        computeLikelihoodFromArrays(iLik, tProb, &(*_likelihoods_father)[fatherFather->getId()], &pxy_[father->getId()], *_likelihoods_node_father, nbSons, nbDistinctSites_, nbClasses_, nbStates_, false);
      }
      else
      {
        computeLikelihoodFromArrays(iLik, tProb, *_likelihoods_node_father, nbSons, nbDistinctSites_, nbClasses_, nbStates_, false);
      }
    }

    if (!father->hasFather())
    {
      // We have to account for the root frequencies:
      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        VVdouble* _likelihoods_node_father_i = &(*_likelihoods_node_father)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* _likelihoods_node_father_i_c = &(*_likelihoods_node_father_i)[c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            (*_likelihoods_node_father_i_c)[x] *= rootFreqs_[x];
          }
        }
      }
    }

    // Call the method on each son node:
    size_t nbNodeSons = node->getNumberOfSons();
    for (size_t i = 0; i < nbNodeSons; i++)
    {
      computeSubtreeLikelihoodPrefix(node->getSon(i)); // Recursive method.
    }
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::computeRootLikelihood()
{
  const Node* root = tree_->getRootNode();
  VVVdouble* rootLikelihoods = &likelihoodData_->getRootLikelihoodArray();
  // Set all likelihoods to 1 for a start:
  if (root->isLeaf())
  {
    VVdouble* leavesLikelihoods_root = &likelihoodData_->getLeafLikelihoods(root->getId());
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble* rootLikelihoods_i = &(*rootLikelihoods)[i];
      Vdouble* leavesLikelihoods_root_i = &(*leavesLikelihoods_root)[i];
      for (size_t c = 0; c < nbClasses_; c++)
      {
        Vdouble* rootLikelihoods_i_c = &(*rootLikelihoods_i)[c];
        for (size_t x = 0; x < nbStates_; x++)
        {
          (*rootLikelihoods_i_c)[x] = (*leavesLikelihoods_root_i)[x];
        }
      }
    }
  }
  else
  {
    resetLikelihoodArray(*rootLikelihoods);
  }

  map<int, VVVdouble>* likelihoods_root = &likelihoodData_->getLikelihoodArrays(root->getId());
  size_t nbNodes = root->getNumberOfSons();
  vector<const VVVdouble*> iLik(nbNodes);
  vector<const VVVdouble*> tProb(nbNodes);
  for (size_t n = 0; n < nbNodes; n++)
  {
    const Node* son = root->getSon(n);
    tProb[n] = &pxy_[son->getId()];
    iLik[n] = &(*likelihoods_root)[son->getId()];
  }
  computeLikelihoodFromArrays(iLik, tProb, *rootLikelihoods, nbNodes, nbDistinctSites_, nbClasses_, nbStates_, false);

  Vdouble p = rateDistribution_->getProbabilities();
  VVdouble* rootLikelihoodsS  = &likelihoodData_->getRootSiteLikelihoodArray();
  Vdouble* rootLikelihoodsSR = &likelihoodData_->getRootRateSiteLikelihoodArray();
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    // For each site in the sequence,
    VVdouble* rootLikelihoods_i = &(*rootLikelihoods)[i];
    Vdouble* rootLikelihoodsS_i = &(*rootLikelihoodsS)[i];
    (*rootLikelihoodsSR)[i] = 0;
    for (size_t c = 0; c < nbClasses_; c++)
    {
      // For each rate classe,
      Vdouble* rootLikelihoods_i_c = &(*rootLikelihoods_i)[c];
      double* rootLikelihoodsS_i_c = &(*rootLikelihoodsS_i)[c];
      (*rootLikelihoodsS_i_c) = 0;
      for (size_t x = 0; x < nbStates_; x++)
      {
        // For each initial state,
        (*rootLikelihoodsS_i_c) += rootFreqs_[x] * (*rootLikelihoods_i_c)[x];
      }
      (*rootLikelihoodsSR)[i] += p[c] * (*rootLikelihoodsS_i_c);
    }

    // Final checking (for numerical errors):
    if ((*rootLikelihoodsSR)[i] < 0)
      (*rootLikelihoodsSR)[i] = 0.;
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::computeLikelihoodAtNode_(const Node* node, VVVdouble& likelihoodArray) const
{
//  const Node * node = tree_->getNode(nodeId);
  int nodeId = node->getId();
  likelihoodArray.resize(nbDistinctSites_);
  map<int, VVVdouble>* likelihoods_node = &likelihoodData_->getLikelihoodArrays(node->getId());

  // Initialize likelihood array:
  if (node->isLeaf())
  {
    VVdouble* leavesLikelihoods_node = &likelihoodData_->getLeafLikelihoods(nodeId);
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble* likelihoodArray_i = &likelihoodArray[i];
      Vdouble* leavesLikelihoods_node_i = &(*leavesLikelihoods_node)[i];
      likelihoodArray_i->resize(nbClasses_);
      for (size_t c = 0; c < nbClasses_; c++)
      {
        Vdouble* likelihoodArray_i_c = &(*likelihoodArray_i)[c];
        likelihoodArray_i_c->resize(nbStates_);
        for (size_t x = 0; x < nbStates_; x++)
        {
          (*likelihoodArray_i_c)[x] = (*leavesLikelihoods_node_i)[x];
        }
      }
    }
  }
  else
  {
    // Otherwise:
    // Set all likelihoods to 1 for a start:
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble* likelihoodArray_i = &likelihoodArray[i];
      likelihoodArray_i->resize(nbClasses_);
      for (size_t c = 0; c < nbClasses_; c++)
      {
        Vdouble* likelihoodArray_i_c = &(*likelihoodArray_i)[c];
        likelihoodArray_i_c->resize(nbStates_);
        for (size_t x = 0; x < nbStates_; x++)
        {
          (*likelihoodArray_i_c)[x] = 1.;
        }
      }
    }
  }

  size_t nbNodes = node->getNumberOfSons();

  vector<const VVVdouble*> iLik(nbNodes);
  vector<const VVVdouble*> tProb(nbNodes);
  for (size_t n = 0; n < nbNodes; n++)
  {
    const Node* son = node->getSon(n);
    tProb[n] = &pxy_[son->getId()];
    iLik[n] = &(*likelihoods_node)[son->getId()];
  }

  if (node->hasFather())
  {
    const Node* father = node->getFather();
    computeLikelihoodFromArrays(iLik, tProb, &(*likelihoods_node)[father->getId()], &pxy_[nodeId], likelihoodArray, nbNodes, nbDistinctSites_, nbClasses_, nbStates_, false);
  }
  else
  {
    computeLikelihoodFromArrays(iLik, tProb, likelihoodArray, nbNodes, nbDistinctSites_, nbClasses_, nbStates_, false);

    // We have to account for the root frequencies:
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble* likelihoodArray_i = &likelihoodArray[i];
      for (size_t c = 0; c < nbClasses_; c++)
      {
        Vdouble* likelihoodArray_i_c = &(*likelihoodArray_i)[c];
        for (size_t x = 0; x < nbStates_; x++)
        {
          (*likelihoodArray_i_c)[x] *= rootFreqs_[x];
        }
      }
    }
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::computeLikelihoodFromArrays(
  const vector<const VVVdouble*>& iLik,
  const vector<const VVVdouble*>& tProb,
  VVVdouble& oLik,
  size_t nbNodes,
  size_t nbDistinctSites,
  size_t nbClasses,
  size_t nbStates,
  bool reset)
{
  if (reset)
    resetLikelihoodArray(oLik);

  for (size_t n = 0; n < nbNodes; n++)
  {
    const VVVdouble* pxy_n = tProb[n];
    const VVVdouble* iLik_n = iLik[n];

    for (size_t i = 0; i < nbDistinctSites; i++)
    {
      // For each site in the sequence,
      const VVdouble* iLik_n_i = &(*iLik_n)[i];
      VVdouble* oLik_i = &(oLik)[i];

      for (size_t c = 0; c < nbClasses; c++)
      {
        // For each rate classe,
        const Vdouble* iLik_n_i_c = &(*iLik_n_i)[c];
        Vdouble* oLik_i_c = &(*oLik_i)[c];
        const VVdouble* pxy_n_c = &(*pxy_n)[c];
        for (size_t x = 0; x < nbStates; x++)
        {
          // For each initial state,
          const Vdouble* pxy_n_c_x = &(*pxy_n_c)[x];
          double likelihood = 0;
          for (size_t y = 0; y < nbStates; y++)
          {
            likelihood += (*pxy_n_c_x)[y] * (*iLik_n_i_c)[y];
          }
          // We store this conditionnal likelihood into the corresponding array:
          (*oLik_i_c)[x] *= likelihood;
        }
      }
    }
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::computeLikelihoodFromArrays(
  const vector<const VVVdouble*>& iLik,
  const vector<const VVVdouble*>& tProb,
  const VVVdouble* iLikR,
  const VVVdouble* tProbR,
  VVVdouble& oLik,
  size_t nbNodes,
  size_t nbDistinctSites,
  size_t nbClasses,
  size_t nbStates,
  bool reset)
{
  if (reset)
    resetLikelihoodArray(oLik);

  for (size_t n = 0; n < nbNodes; n++)
  {
    const VVVdouble* pxy_n = tProb[n];
    const VVVdouble* iLik_n = iLik[n];

    for (size_t i = 0; i < nbDistinctSites; i++)
    {
      // For each site in the sequence,
      const VVdouble* iLik_n_i = &(*iLik_n)[i];
      VVdouble* oLik_i = &(oLik)[i];

      for (size_t c = 0; c < nbClasses; c++)
      {
        // For each rate classe,
        const Vdouble* iLik_n_i_c = &(*iLik_n_i)[c];
        Vdouble* oLik_i_c = &(*oLik_i)[c];
        const VVdouble* pxy_n_c = &(*pxy_n)[c];
        for (size_t x = 0; x < nbStates; x++)
        {
          // For each initial state,
          const Vdouble* pxy_n_c_x = &(*pxy_n_c)[x];
          double likelihood = 0;
          for (size_t y = 0; y < nbStates; y++)
          {
            // cout << "1:" << (* pxy_n_c_x)[y]  << endl;
            // cout << "2:" << (* iLik_n_i_c)[y] << endl;
            likelihood += (*pxy_n_c_x)[y] * (*iLik_n_i_c)[y];
            // cout << i << "\t" << c << "\t" << x << "\t" << y << "\t" <<  (* pxy__son_c_x)[y] << "\t" << (* likelihoods_root_son_i_c)[y] << endl;
          }
          // We store this conditionnal likelihood into the corresponding array:
          (*oLik_i_c)[x] *= likelihood;
        }
      }
    }
  }

  // Now deal with the subtree containing the root:
  for (size_t i = 0; i < nbDistinctSites; i++)
  {
    // For each site in the sequence,
    const VVdouble* iLikR_i = &(*iLikR)[i];
    VVdouble* oLik_i = &(oLik)[i];

    for (size_t c = 0; c < nbClasses; c++)
    {
      // For each rate classe,
      const Vdouble* iLikR_i_c = &(*iLikR_i)[c];
      Vdouble* oLik_i_c = &(*oLik_i)[c];
      const VVdouble* pxyR_c = &(*tProbR)[c];
      for (size_t x = 0; x < nbStates; x++)
      {
        double likelihood = 0;
        for (size_t y = 0; y < nbStates; y++)
        {
          // For each final state,
          const Vdouble* pxyR_c_y = &(*pxyR_c)[y];
          likelihood += (*pxyR_c_y)[x] * (*iLikR_i_c)[y];
        }
        // We store this conditionnal likelihood into the corresponding array:
        (*oLik_i_c)[x] *= likelihood;
      }
    }
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::displayLikelihood(const Node* node)
{
  cout << "Likelihoods at node " << node->getId() << ": " << endl;
  for (size_t n = 0; n < node->getNumberOfSons(); n++)
  {
    const Node* subNode = node->getSon(n);
    cout << "Array for sub-node " << subNode->getId() << endl;
    displayLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId(), subNode->getId()));
  }
  if (node->hasFather())
  {
    const Node* father = node->getFather();
    cout << "Array for father node " << father->getId() << endl;
    displayLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId(), father->getId()));
  }
  cout << "                                         ***" << endl;
}

/*******************************************************************************/

