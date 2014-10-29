//
// File: DRHomogeneousTreeLikelihood.cpp
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

#include "DRHomogeneousTreeLikelihood.h"
#include "../PatternTools.h"

// From SeqLib:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

DRHomogeneousTreeLikelihood::DRHomogeneousTreeLikelihood(
  const Tree& tree,
  SubstitutionModel* model,
  DiscreteDistribution* rDist,
  bool checkRooted,
  bool verbose)
throw (Exception) :
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  likelihoodData_(0),
  minusLogLik_(-1.)
{
  init_();
}

/******************************************************************************/

DRHomogeneousTreeLikelihood::DRHomogeneousTreeLikelihood(
  const Tree& tree,
  const SiteContainer& data,
  SubstitutionModel* model,
  DiscreteDistribution* rDist,
  bool checkRooted,
  bool verbose)
throw (Exception) :
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  likelihoodData_(0),
  minusLogLik_(-1.)
{
  init_();
  setData(data);
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::init_() throw (Exception)
{
  likelihoodData_ = new DRASDRTreeLikelihoodData(
    tree_,
    rateDistribution_->getNumberOfCategories());
}

/******************************************************************************/

DRHomogeneousTreeLikelihood::DRHomogeneousTreeLikelihood(const DRHomogeneousTreeLikelihood& lik) :
  AbstractHomogeneousTreeLikelihood(lik),
  likelihoodData_(0),
  minusLogLik_(-1.)
{
  likelihoodData_ = dynamic_cast<DRASDRTreeLikelihoodData*>(lik.likelihoodData_->clone());
  likelihoodData_->setTree(tree_);
  minusLogLik_ = lik.minusLogLik_;
}

/******************************************************************************/

DRHomogeneousTreeLikelihood& DRHomogeneousTreeLikelihood::operator=(const DRHomogeneousTreeLikelihood& lik)
{
  AbstractHomogeneousTreeLikelihood::operator=(lik);
  if (likelihoodData_)
    delete likelihoodData_;
  likelihoodData_ = dynamic_cast<DRASDRTreeLikelihoodData*>(lik.likelihoodData_->clone());
  likelihoodData_->setTree(tree_);
  minusLogLik_ = lik.minusLogLik_;
  return *this;
}

/******************************************************************************/

DRHomogeneousTreeLikelihood::~DRHomogeneousTreeLikelihood()
{
  delete likelihoodData_;
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::setData(const SiteContainer& sites) throw (Exception)
{
  if (data_)
    delete data_;
  data_ = PatternTools::getSequenceSubset(sites, *tree_->getRootNode());
  if (verbose_)
    ApplicationTools::displayTask("Initializing data structure");
  likelihoodData_->initLikelihoods(*data_, *model_);
  if (verbose_)
    ApplicationTools::displayTaskDone();

  nbSites_ = likelihoodData_->getNumberOfSites();
  nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();
  nbStates_ = likelihoodData_->getNumberOfStates();

  if (verbose_)
    ApplicationTools::displayResult("Number of distinct sites",
                                    TextTools::toString(nbDistinctSites_));
  initialized_ = false;
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLikelihood() const
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

double DRHomogeneousTreeLikelihood::getLogLikelihood() const
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

double DRHomogeneousTreeLikelihood::getLikelihoodForASite(size_t site) const
{
  return likelihoodData_->getRootRateSiteLikelihoodArray()[likelihoodData_->getRootArrayPosition(site)];
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihoodForASite(size_t site) const
{
  return log(likelihoodData_->getRootRateSiteLikelihoodArray()[likelihoodData_->getRootArrayPosition(site)]);
}

/******************************************************************************/
double DRHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  return likelihoodData_->getRootSiteLikelihoodArray()[likelihoodData_->getRootArrayPosition(site)][rateClass];
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  return log(likelihoodData_->getRootSiteLikelihoodArray()[likelihoodData_->getRootArrayPosition(site)][rateClass]);
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  return likelihoodData_->getRootLikelihoodArray()[likelihoodData_->getRootArrayPosition(site)][rateClass][static_cast<size_t>(state)];
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  return log(likelihoodData_->getRootLikelihoodArray()[likelihoodData_->getRootArrayPosition(site)][rateClass][static_cast<size_t>(state)]);
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::setParameters(const ParameterList& parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::fireParameterChanged(const ParameterList& params)
{
  applyParameters();

  if (rateDistribution_->getParameters().getCommonParametersWith(params).size() > 0
      || model_->getParameters().getCommonParametersWith(params).size() > 0)
  {
    // Rate parameter changed, need to recompute all probs:
    computeAllTransitionProbabilities();
  }
  else if (params.size() > 0)
  {
    // We may save some computations:
    for (size_t i = 0; i < params.size(); i++)
    {
      string s = params[i].getName();
      if (s.substr(0, 5) == "BrLen")
      {
        // Branch length parameter:
        computeTransitionProbabilitiesForNode(nodes_[TextTools::to < size_t > (s.substr(5))]);
      }
    }
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

  minusLogLik_ = -getLogLikelihood();
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getValue() const
throw (Exception)
{
  if (!isInitialized())
    throw Exception("DRHomogeneousTreeLikelihood::getValue(). Instance is not initialized.");
  return minusLogLik_;
}

/******************************************************************************
*                           First Order Derivatives                          *
******************************************************************************/
void DRHomogeneousTreeLikelihood::computeTreeDLikelihoodAtNode(const Node* node)
{
  const Node* father = node->getFather();
  VVVdouble* likelihoods_father_node = &likelihoodData_->getLikelihoodArray(father->getId(), node->getId());
  Vdouble* dLikelihoods_node = &likelihoodData_->getDLikelihoodArray(node->getId());
  VVVdouble* dpxy_node = &dpxy_[node->getId()];
  VVVdouble larray;
  computeLikelihoodAtNode_(father, larray, node);

  Vdouble* rootLikelihoodsSR = &likelihoodData_->getRootRateSiteLikelihoodArray();

  double dLi, dLic, dLicx;

  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    VVdouble* likelihoods_father_node_i = &(*likelihoods_father_node)[i];
    VVdouble* larray_i = &larray[i];
    dLi = 0;
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* likelihoods_father_node_i_c = &(*likelihoods_father_node_i)[c];
      Vdouble* larray_i_c = &(*larray_i)[c];
      VVdouble* dpxy_node_c = &(*dpxy_node)[c];
      dLic = 0;
      for (size_t x = 0; x < nbStates_; x++)
      {
        Vdouble* dpxy_node_c_x = &(*dpxy_node_c)[x];
        dLicx = 0;
        for (size_t y = 0; y < nbStates_; y++)
        {
          dLicx += (*dpxy_node_c_x)[y] * (*likelihoods_father_node_i_c)[y];
        }
        dLicx *= (*larray_i_c)[x];
        dLic += dLicx;
      }
      dLi += rateDistribution_->getProbability(c) * dLic;
    }

    (*dLikelihoods_node)[i] = dLi / (*rootLikelihoodsSR)[i];
    // cout << dLi << "\t" << (*rootLikelihoodsSR)[i] << endl;
  }
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeTreeDLikelihoods()
{
  for (size_t k = 0; k < nbNodes_; k++)
  {
    computeTreeDLikelihoodAtNode(nodes_[k]);
  }
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getFirstOrderDerivative(const std::string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("DRHomogeneousTreeLikelihood::getFirstOrderDerivative().", variable);
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

  // Get the node with the branch whose length must be derivated:
  size_t brI = TextTools::to<size_t>(variable.substr(5));
  const Node* branch = nodes_[brI];
  Vdouble* dLikelihoods_branch = &likelihoodData_->getDLikelihoodArray(branch->getId());
  double d = 0;
  const vector<unsigned int>* w = &likelihoodData_->getWeights();
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    d += (*w)[i] * (*dLikelihoods_branch)[i];
  }

  return -d;
}

/******************************************************************************
*                           Second Order Derivatives                         *
******************************************************************************/
void DRHomogeneousTreeLikelihood::computeTreeD2LikelihoodAtNode(const Node* node)
{
  const Node* father = node->getFather();
  VVVdouble* likelihoods_father_node = &likelihoodData_->getLikelihoodArray(father->getId(), node->getId());
  Vdouble* d2Likelihoods_node = &likelihoodData_->getD2LikelihoodArray(node->getId());
  VVVdouble* d2pxy_node = &d2pxy_[node->getId()];
  VVVdouble larray;
  computeLikelihoodAtNode_(father, larray, node);
  Vdouble* rootLikelihoodsSR = &likelihoodData_->getRootRateSiteLikelihoodArray();

  double d2Li, d2Lic, d2Licx;

  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    VVdouble* likelihoods_father_node_i = &(*likelihoods_father_node)[i];
    VVdouble* larray_i = &larray[i];
    d2Li = 0;
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* likelihoods_father_node_i_c = &(*likelihoods_father_node_i)[c];
      Vdouble* larray_i_c = &(*larray_i)[c];
      VVdouble* d2pxy_node_c = &(*d2pxy_node)[c];
      d2Lic = 0;
      for (size_t x = 0; x < nbStates_; x++)
      {
        Vdouble* d2pxy_node_c_x = &(*d2pxy_node_c)[x];
        d2Licx = 0;
        for (size_t y = 0; y < nbStates_; y++)
        {
          d2Licx += (*d2pxy_node_c_x)[y] * (*likelihoods_father_node_i_c)[y];
        }
        d2Licx *= (*larray_i_c)[x];
        d2Lic += d2Licx;
      }
      d2Li += rateDistribution_->getProbability(c) * d2Lic;
    }
    (*d2Likelihoods_node)[i] = d2Li / (*rootLikelihoodsSR)[i];
  }
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeTreeD2Likelihoods()
{
  for (size_t k = 0; k < nbNodes_; k++)
  {
    computeTreeD2LikelihoodAtNode(nodes_[k]);
  }
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getSecondOrderDerivative(const std::string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("DRHomogeneousTreeLikelihood::getSecondOrderDerivative().", variable);
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

  // Get the node with the branch whose length must be derivated:
  size_t brI = TextTools::to<size_t>(variable.substr(5));
  const Node* branch = nodes_[brI];
  Vdouble* _dLikelihoods_branch = &likelihoodData_->getDLikelihoodArray(branch->getId());
  Vdouble* _d2Likelihoods_branch = &likelihoodData_->getD2LikelihoodArray(branch->getId());
  double d2 = 0;
  const vector<unsigned int>* w = &likelihoodData_->getWeights();
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    d2 += (*w)[i] * ((*_d2Likelihoods_branch)[i] - pow((*_dLikelihoods_branch)[i], 2));
  }
  return -d2;
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::resetLikelihoodArrays(const Node* node)
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

void DRHomogeneousTreeLikelihood::computeTreeLikelihood()
{
  computeSubtreeLikelihoodPostfix(tree_->getRootNode());
  computeSubtreeLikelihoodPrefix(tree_->getRootNode());
  computeRootLikelihood();
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeSubtreeLikelihoodPostfix(const Node* node)
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

void DRHomogeneousTreeLikelihood::computeSubtreeLikelihoodPrefix(const Node* node)
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

      size_t nbSons = nodes.size(); // In case of a bifurcating tree, this is equal to 1, excepted for the root.

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

void DRHomogeneousTreeLikelihood::computeRootLikelihood()
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

void DRHomogeneousTreeLikelihood::computeLikelihoodAtNode_(const Node* node, VVVdouble& likelihoodArray, const Node* sonNode) const
{
  // const Node * node = tree_->getNode(nodeId);
  int nodeId = node->getId();
  likelihoodArray.resize(nbDistinctSites_);
  map<int, VVVdouble>* likelihoods_node = &likelihoodData_->getLikelihoodArrays(nodeId);

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

  vector<const VVVdouble*> iLik;
  vector<const VVVdouble*> tProb;
  bool test = false;
  for (size_t n = 0; n < nbNodes; n++)
  {
    const Node* son = node->getSon(n);
    if (son != sonNode) {
      tProb.push_back(&pxy_[son->getId()]);
      iLik.push_back(&(*likelihoods_node)[son->getId()]);
    } else {
      test = true;
    }
  }
  if (sonNode) {
    if (test)
      nbNodes--;
    else
      throw Exception("DRHomogeneousTreeLikelihood::computeLikelihoodAtNode_(...). 'sonNode' not found as a son of 'node'.");
  }

  if (node->hasFather())
  {
    const Node* father = node->getFather();
    computeLikelihoodFromArrays(iLik, tProb, &(*likelihoods_node)[father->getId()], &pxy_[nodeId], likelihoodArray, nbNodes, nbDistinctSites_, nbClasses_, nbStates_, false);
  }
  else
  {
    computeLikelihoodFromArrays(iLik, tProb, likelihoodArray, nbNodes, nbDistinctSites_, nbClasses_, nbStates_, false);

    // We have to account for the equilibrium frequencies:
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

void DRHomogeneousTreeLikelihood::computeLikelihoodFromArrays(
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

void DRHomogeneousTreeLikelihood::computeLikelihoodFromArrays(
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

void DRHomogeneousTreeLikelihood::displayLikelihood(const Node* node)
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

