//
// File: SingleRecursiveTreeLikelihoodCalculation.cpp
// Created by: Julien Dutheil
// Created on: Tue May 15 14:30 2012
// From file: RNonHomogeneousTreeLikelihood.cpp
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

#include "SingleRecursiveTreeLikelihoodCalculation.h"
#include "../PatternTools.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/App/ApplicationTools.h>

using namespace bpp;
using namespace newlik;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

SingleRecursiveTreeLikelihoodCalculation::SingleRecursiveTreeLikelihoodCalculation(
  SubstitutionProcess* process,
  bool verbose,
  bool usePatterns)
throw (Exception) :
  AbstractTreeLikelihoodCalculation(process),
  likelihoodData_(0),
  root1_(-1),
  root2_(-1)
{
  init_(usePatterns);
}

/******************************************************************************/

SingleRecursiveTreeLikelihoodCalculation::SingleRecursiveTreeLikelihoodCalculation(
  const SiteContainer& data,
  SubstitutionProcess* process,
  bool verbose,
  bool usePatterns)
throw (Exception) :
  AbstractTreeLikelihoodCalculation(process),
  likelihoodData_(0),
  root1_(-1),
  root2_(-1)
{
  init_(usePatterns);
  setData(data);
}

/******************************************************************************/

void SingleRecursiveTreeLikelihoodCalculation::init_(bool usePatterns) throw (Exception)
{
  likelihoodData_.reset(new SingleRecursiveTreeLikelihoodData(
                          process_->getNumberOfClasses(),
                          usePatterns));
}

/******************************************************************************/

SingleRecursiveTreeLikelihoodCalculation::SingleRecursiveTreeLikelihoodCalculation(const SingleRecursiveTreeLikelihoodCalculation& tlc) :
  AbstractTreeLikelihoodCalculation(tlc),
  likelihoodData_(0),
  root1_(tlc.root1_),
  root2_(tlc.root2_)
{
  likelihoodData_.reset(tlc.likelihoodData_->clone());
}

/******************************************************************************/

SingleRecursiveTreeLikelihoodCalculation& SingleRecursiveTreeLikelihoodCalculation::operator=(
  const SingleRecursiveTreeLikelihoodCalculation& tlc)
{
  AbstractTreeLikelihoodCalculation::operator=(tlc);
  likelihoodData_.reset(tlc.likelihoodData_->clone());
  root1_ = tlc.root1_;
  root2_ = tlc.root2_;
  return *this;
}

/******************************************************************************/

void SingleRecursiveTreeLikelihoodCalculation::setData(const SiteContainer& sites) throw (Exception)
{
  // try
  // {
    const TreeTemplate<Node>& tt = dynamic_cast<const TreeTemplate<Node>&>(process_->getTree());
    data_.reset(PatternTools::getSequenceSubset(sites, *tt.getRootNode()));
  // }
  //     catch (std::bad_cast)
  // {
  //   throw Exception("DEBUG. SingleRecursiveTreeLikelihoodCalculation::setData. The SubstitutionProcess does not use a TreeTemplate object.");
  // }
  if (verbose_)
    ApplicationTools::displayTask("Initializing data structure");
  likelihoodData_->initLikelihoods(*data_, *process_); // We assume here that all models have the same number of states, and that they have the same 'init' method,
                                                       // Which is a reasonable assumption as long as they share the same alphabet.
  if (verbose_)
    ApplicationTools::displayTaskDone();

  nbSites_         = likelihoodData_->getNumberOfSites();
  nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();
  nbStates_        = likelihoodData_->getNumberOfStates();

  if (verbose_)
    ApplicationTools::displayResult("Number of distinct sites", TextTools::toString(nbDistinctSites_));

  initialized_ = true;

  // Recompute likelihood:
  computeTreeLikelihood();
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getLikelihoodForASite(size_t site) const
{
  double l = 0;
  VVdouble* la = &likelihoodData_->getLikelihoodArray(
    process_->getTree().getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)];
  for (size_t i = 0; i < nbClasses_; ++i)
  {
    for (size_t j = 0; j < nbStates_; ++j)
    {
      l += (*la)[i][j] * process_->getProbabilityForModel(i) * process_->getRootFrequencies()[j];
    }
  }
  if (l < 0) l = 0; //May happen because of numerical errors.
  return l;
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getLikelihoodForASiteForAClass(size_t site, size_t classIndex) const
{
  double l = 0;
  Vdouble* la = &likelihoodData_->getLikelihoodArray(
    process_->getTree().getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][classIndex];
  for (size_t i = 0; i < nbStates_; ++i)
  {
    l += (*la)[i] * process_->getRootFrequencies()[i];
  }
  return l;
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getLikelihoodForASiteForAState(size_t site, int state) const
{
  double l = 0;
  VVdouble* la = &likelihoodData_->getLikelihoodArray(
    process_->getTree().getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)];
  for (size_t c = 0; c < nbClasses_; ++c)
  {
    l += (*la)[c][state] * process_->getProbabilityForModel(c);
  }
  return l;
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state) const
{
  return likelihoodData_->getLikelihoodArray(
           process_->getTree().getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][classIndex][state];
}


/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getDLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state) const
{
  return likelihoodData_->getDLikelihoodArray(
           process_->getTree().getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][classIndex][state];
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getD2LikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state) const
{
  return likelihoodData_->getD2LikelihoodArray(
           process_->getTree().getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][classIndex][state];
}

/******************************************************************************
 *                           Likelihood computation                           *
 ******************************************************************************/

void SingleRecursiveTreeLikelihoodCalculation::computeTreeLikelihood()
{
  computeSubtreeLikelihood_(process_->getTree().getRootNode());
}

/******************************************************************************/

void SingleRecursiveTreeLikelihoodCalculation::computeSubtreeLikelihood_(const Node* node)
{
  if (node->isLeaf())
    return;
  size_t nbSites  = likelihoodData_->getLikelihoodArray(node->getId()).size();
  size_t nbNodes  = node->getNumberOfSons();

  // Must reset the likelihood array first (i.e. set all of them to 1):
  VVVdouble* likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());
  for (size_t i = 0; i < nbSites; i++)
  {
    // For each site in the sequence,
    VVdouble* likelihoods_node_i = &(*likelihoods_node)[i];
    for (size_t c = 0; c < nbClasses_; c++)
    {
      // For each rate classe,
      Vdouble* likelihoods_node_i_c = &(*likelihoods_node_i)[c];
      for (size_t x = 0; x < nbStates_; x++)
      {
        // For each initial state,
        (*likelihoods_node_i_c)[x] = 1.;
      }
    }
  }

  for (size_t l = 0; l < nbNodes; l++)
  {
    // For each son node,
    const Node* son = node->getSon(l);

    computeSubtreeLikelihood_(son); // Recursive method.

    vector<size_t>* patternLinks_node_son = &likelihoodData_->getArrayPositions(node->getId(), son->getId());
    VVVdouble* likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

    // Get all transition probabilities:

    vector<const Matrix<double>*> pxy_son(nbClasses_);
    for (size_t c = 0; c < nbClasses_; ++c)
    {
      pxy_son[c] = &process_->getTransitionProbabilities(son->getId(), c);
    }

    // Loop over sites:
    for (size_t i = 0; i < nbSites; i++)
    {
      // For each site in the sequence,
      VVdouble* likelihoods_son_i = &(*likelihoods_son)[(*patternLinks_node_son)[i]];
      VVdouble* likelihoods_node_i = &(*likelihoods_node)[i];
      for (size_t c = 0; c < nbClasses_; c++)
      {
        // For each rate classe,
        Vdouble* likelihoods_son_i_c = &(*likelihoods_son_i)[c];
        Vdouble* likelihoods_node_i_c = &(*likelihoods_node_i)[c];
        for (size_t x = 0; x < nbStates_; x++)
        {
          // For each initial state,
          double likelihood = 0;
          for (size_t y = 0; y < nbStates_; y++)
          {
            likelihood += (*pxy_son[c])(x, y) * (*likelihoods_son_i_c)[y];
          }
          (*likelihoods_node_i_c)[x] *= likelihood;
        }
      }
    }
  }
}


/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/

void SingleRecursiveTreeLikelihoodCalculation::computeTreeDLikelihood(const string& variable)
{
/*  if (variable == "BrLenRoot")
   {
    const Node* father = process_->getTree().getRootNode();

    // Compute dLikelihoods array for the father node.
    // First initialize to 1:
    VVVdouble* dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(father->getId());
    size_t nbSites = dLikelihoods_father->size();
    for (size_t i = 0; i < nbSites; i++)
    {
      VVdouble* dLikelihoods_father_i = &(*dLikelihoods_father)[i];
      for (size_t c = 0; c < nbClasses_; c++)
      {
        Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
        for (size_t s = 0; s < nbStates_; s++)
        {
          (*dLikelihoods_father_i_c)[s] = 1.;
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
        vector<size_t>* _patternLinks_fatherroot1_ = &likelihoodData_->getArrayPositions(father->getId(), root1->getId());
        vector<size_t>* _patternLinks_fatherroot2_ = &likelihoodData_->getArrayPositions(father->getId(), root2->getId());
        VVVdouble* likelihoodsroot1_ = &likelihoodData_->getLikelihoodArray(root1->getId());
        VVVdouble* likelihoodsroot2_ = &likelihoodData_->getLikelihoodArray(root2->getId());
        double pos = getParameterValue("RootPosition");

        VVVdouble* dpxyRoot1_  = &dpxy_[root1_];
        VVVdouble* dpxyRoot2_  = &dpxy_[root2_];
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
        vector<size_t>* _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
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
    // First initialize to 1:
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
        vector<size_t>* _patternLinks_fatherroot1_ = &likelihoodData_->getArrayPositions(father->getId(), root1->getId());
        vector<size_t>* _patternLinks_fatherroot2_ = &likelihoodData_->getArrayPositions(father->getId(), root2->getId());
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
        vector<size_t>* _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
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
   }*/

  // Get the node with the branch whose length must be derivated:
  int brId = atoi(variable.substr(5).c_str());
  const Node* branch = process_->getTree().getNode(brId);
  const Node* father = branch->getFather();

  // Compute dLikelihoods array for the father node.
  // First initialize to 1:
  VVVdouble* dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(father->getId());
  size_t nbSites = dLikelihoods_father->size();
  for (size_t i = 0; i < nbSites; ++i)
  {
    VVdouble* dLikelihoods_father_i = &(*dLikelihoods_father)[i];
    for (size_t c = 0; c < nbClasses_; ++c)
    {
      Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
      for (size_t s = 0; s < nbStates_; ++s)
      {
        (*dLikelihoods_father_i_c)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; ++l)
  {
    const Node* son = father->getSon(l);

    vector<size_t>* patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
    VVVdouble* likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

    if (son == branch)
    {
      // Get all derivatives of transition probabilities:
      vector<const Matrix<double>*> dpxy_son(nbClasses_);
      for (size_t c = 0; c < nbClasses_; ++c)
      {
        dpxy_son[c] = &process_->getTransitionProbabilitiesD1(son->getId(), c);
      }

      // Loop over sites:
      for (size_t i = 0; i < nbSites; ++i)
      {
        VVdouble* likelihoods_son_i = &(*likelihoods_son)[(*patternLinks_father_son)[i]];
        VVdouble* dLikelihoods_father_i = &(*dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; ++c)
        {
          Vdouble* likelihoods_son_i_c = &(*likelihoods_son_i)[c];
          Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
          for (size_t x = 0; x < nbStates_; ++x)
          {
            double dl = 0;
            for (size_t y = 0; y < nbStates_; ++y)
            {
              dl += (*dpxy_son[c])(x, y) * (*likelihoods_son_i_c)[y];
            }
            (*dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
    else
    {
      // Get all transition probabilities:
      vector<const Matrix<double>*> pxy_son(nbClasses_);
      for (size_t c = 0; c < nbClasses_; ++c)
      {
        pxy_son[c] = &process_->getTransitionProbabilities(son->getId(), c);
      }

      // Loop over sites:
      for (size_t i = 0; i < nbSites; ++i)
      {
        VVdouble* likelihoods_son_i = &(*likelihoods_son)[(*patternLinks_father_son)[i]];
        VVdouble* dLikelihoods_father_i = &(*dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; ++c)
        {
          Vdouble* likelihoods_son_i_c = &(*likelihoods_son_i)[c];
          Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
          for (size_t x = 0; x < nbStates_; ++x)
          {
            double dl = 0;
            for (size_t y = 0; y < nbStates_; ++y)
            {
              dl += (*pxy_son[c])(x, y) * (*likelihoods_son_i_c)[y];
            }
            (*dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
  }

  // Now we go down the tree toward the root node:
  computeDownSubtreeDLikelihood_(father);
}

/******************************************************************************/

void SingleRecursiveTreeLikelihoodCalculation::computeDownSubtreeDLikelihood_(const Node* node)
{
  const Node* father = node->getFather();
  // We assume that the _dLikelihoods array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if (father == 0)
    return;  // We reached the root!

  // Compute dLikelihoods array for the father node.
  // First initialize to 1:
  VVVdouble* dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(father->getId());
  size_t nbSites = dLikelihoods_father->size();
  for (size_t i = 0; i < nbSites; ++i)
  {
    VVdouble* dLikelihoods_father_i = &(*dLikelihoods_father)[i];
    for (size_t c = 0; c < nbClasses_; ++c)
    {
      Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
      for (size_t s = 0; s < nbStates_; ++s)
      {
        (*dLikelihoods_father_i_c)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; ++l)
  {
    const Node* son = father->getSon(l);

    vector<size_t>* patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());

    // Get all transition probabilities:
    vector<const Matrix<double>*> pxy_son(nbClasses_);
    for (size_t c = 0; c < nbClasses_; ++c)
    {
      pxy_son[c] = &process_->getTransitionProbabilities(son->getId(), c);
    }

    if (son == node)
    {
      VVVdouble* dLikelihoods_son = &likelihoodData_->getDLikelihoodArray(son->getId());

      // Loop over sites:
      for (size_t i = 0; i < nbSites; ++i)
      {
        VVdouble* dLikelihoods_son_i = &(*dLikelihoods_son)[(*patternLinks_father_son)[i]];
        VVdouble* dLikelihoods_father_i = &(*dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; ++c)
        {
          Vdouble* dLikelihoods_son_i_c = &(*dLikelihoods_son_i)[c];
          Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
          for (size_t x = 0; x < nbStates_; ++x)
          {
            double dl = 0;
            for (size_t y = 0; y < nbStates_; ++y)
            {
              dl += (*pxy_son[c])(x, y) * (*dLikelihoods_son_i_c)[y];
            }
            (*dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
    else
    {
      VVVdouble* likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

      // Loop over sites:
      for (size_t i = 0; i < nbSites; ++i)
      {
        VVdouble* likelihoods_son_i = &(*likelihoods_son)[(*patternLinks_father_son)[i]];
        VVdouble* dLikelihoods_father_i = &(*dLikelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; ++c)
        {
          Vdouble* likelihoods_son_i_c = &(*likelihoods_son_i)[c];
          Vdouble* dLikelihoods_father_i_c = &(*dLikelihoods_father_i)[c];
          for (size_t x = 0; x < nbStates_; ++x)
          {
            double dl = 0;
            for (size_t y = 0; y < nbStates_; ++y)
            {
              dl += (*pxy_son[c])(x, y) * (*likelihoods_son_i_c)[y];
            }
            (*dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
  }

  // Next step: move toward grand father...
  computeDownSubtreeDLikelihood_(father);
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getDLikelihoodForASiteForAClass(size_t site, size_t classIndex) const
{
  vector<double> rootFreqs = process_->getRootFrequencies();
  Vdouble* dla = &likelihoodData_->getDLikelihoodArray(
           process_->getTree().getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][classIndex];
  double dl = 0;
  for (size_t i = 0; i < nbStates_; ++i)
  {
    dl += (*dla)[i] * rootFreqs[i];
  }
  return dl;
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getDLikelihoodForASite(size_t site) const
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t i = 0; i < nbClasses_; i++)
  {
    dl += getDLikelihoodForASiteForAClass(site, i) * process_->getProbabilityForModel(i);
  }
  return dl;
//  vector<double> rootFreqs = process_->getRootFrequencies();
//  VVdouble* dla = &likelihoodData_->getDLikelihoodArray(
//           process_->getTree().getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)];
//  // Derivative of the sum is the sum of derivatives:
//  double dl = 0;
//  for (size_t i = 0; i < nbClasses_; ++i)
//  {
//    for (size_t j = 0; j < nbStates_; ++j)
//    {
//      dl += (*dla)[i][j] * process_->getProbabilityForModel(i) * rootFreqs[j];
//    }
//  }
//  cout << site << "\t" << dl << endl;
//  return dl;
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getDLogLikelihoodForASite(size_t site) const
{
  // d(f(g(x)))/dx = dg(x)/dx . df(g(x))/dg :
  return getDLikelihoodForASite(site) / getLikelihoodForASite(site);
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getDLogLikelihood() const
{
  // Derivative of the sum is the sum of derivatives:
  vector<double> dla(nbSites_);
  for (size_t i = 0; i < nbSites_; ++i)
  {
    dla[i] = getDLogLikelihoodForASite(i);
  }
  sort(dla.begin(), dla.end());
  double dl = 0;
  for (size_t i = nbSites_; i > 0; --i)
  {
    dl += dla[i - 1];
  }
  return dl;
}

/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/

void SingleRecursiveTreeLikelihoodCalculation::computeTreeD2Likelihood(const string& variable)
{
  /*if (variable == "BrLenRoot")
     {
     const Node* father = tree_->getRootNode();

     // Compute dLikelihoods array for the father node.
     // First initialize to 1:
     VVVdouble* _d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
     unsigned int nbSites  = _d2Likelihoods_father->size();
     for (unsigned int i = 0; i < nbSites; i++)
     {
      VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
      for (unsigned int c = 0; c < nbClasses_; c++)
      {
        Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
        for (unsigned int s = 0; s < nbStates_; s++)
        {
          (*_d2Likelihoods_father_i_c)[s] = 1.;
        }
      }
     }

     unsigned int nbNodes = father->getNumberOfSons();
     for (unsigned int l = 0; l < nbNodes; l++)
     {
      const Node* son = father->getSon(l);

      if (son->getId() == root1_)
      {
        const Node* root1 = father->getSon(0);
        const Node* root2 = father->getSon(1);
        vector<unsigned int> * _patternLinks_fatherroot1_ = &likelihoodData_->getArrayPositions(father->getId(), root1->getId());
        vector<unsigned int> * _patternLinks_fatherroot2_ = &likelihoodData_->getArrayPositions(father->getId(), root2->getId());
        VVVdouble* _likelihoodsroot1_ = &likelihoodData_->getLikelihoodArray(root1->getId());
        VVVdouble* _likelihoodsroot2_ = &likelihoodData_->getLikelihoodArray(root2->getId());
        double pos = getParameterValue("RootPosition");

        VVVdouble* d2pxy_root1_ = &d2pxy_[root1_];
        VVVdouble* d2pxy_root2_ = &d2pxy_[root2_];
        VVVdouble* dpxy_root1_  = &dpxy_[root1_];
        VVVdouble* dpxy_root2_  = &dpxy_[root2_];
        VVVdouble* pxy_root1_   = &pxy_[root1_];
        VVVdouble* pxy_root2_   = &pxy_[root2_];
        for (unsigned int i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoodsroot1__i = &(*_likelihoodsroot1_)[(*_patternLinks_fatherroot1_)[i]];
          VVdouble* _likelihoodsroot2__i = &(*_likelihoodsroot2_)[(*_patternLinks_fatherroot2_)[i]];
          VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
          for (unsigned int c = 0; c < nbClasses_; c++)
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
            for (unsigned int x = 0; x < nbStates_; x++)
            {
              Vdouble* d2pxy_root1__c_x = &(*d2pxy_root1__c)[x];
              Vdouble* d2pxy_root2__c_x = &(*d2pxy_root2__c)[x];
              Vdouble* dpxy_root1__c_x  = &(*dpxy_root1__c)[x];
              Vdouble* dpxy_root2__c_x  = &(*dpxy_root2__c)[x];
              Vdouble* pxy_root1__c_x   = &(*pxy_root1__c)[x];
              Vdouble* pxy_root2__c_x   = &(*pxy_root2__c)[x];
              double d2l1 = 0, d2l2 = 0, dl1 = 0, dl2 = 0, l1 = 0, l2 = 0;
              for (unsigned int y = 0; y < nbStates_; y++)
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
        vector<unsigned int> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
        VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

        VVVdouble* pxy__son = &pxy_[son->getId()];
        for (unsigned int i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
          VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
          for (unsigned int c = 0; c < nbClasses_; c++)
          {
            Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
            Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
            VVdouble* pxy__son_c = &(*pxy__son)[c];
            for (unsigned int x = 0; x < nbStates_; x++)
            {
              double d2l = 0;
              Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
              for (unsigned int y = 0; y < nbStates_; y++)
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
     // First initialize to 1:
     VVVdouble* _d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
     unsigned int nbSites  = _d2Likelihoods_father->size();
     for (unsigned int i = 0; i < nbSites; i++)
     {
      VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
      for (unsigned int c = 0; c < nbClasses_; c++)
      {
        Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
        for (unsigned int s = 0; s < nbStates_; s++)
        {
          (*_d2Likelihoods_father_i_c)[s] = 1.;
        }
      }
     }

     unsigned int nbNodes = father->getNumberOfSons();
     for (unsigned int l = 0; l < nbNodes; l++)
     {
      const Node* son = father->getSon(l);

      if (son->getId() == root1_)
      {
        const Node* root1 = father->getSon(0);
        const Node* root2 = father->getSon(1);
        vector<unsigned int> * _patternLinks_fatherroot1_ = &likelihoodData_->getArrayPositions(father->getId(), root1->getId());
        vector<unsigned int> * _patternLinks_fatherroot2_ = &likelihoodData_->getArrayPositions(father->getId(), root2->getId());
        VVVdouble* _likelihoodsroot1_ = &likelihoodData_->getLikelihoodArray(root1->getId());
        VVVdouble* _likelihoodsroot2_ = &likelihoodData_->getLikelihoodArray(root2->getId());
        double len = getParameterValue("BrLenRoot");

        VVVdouble* d2pxy_root1_ = &d2pxy_[root1_];
        VVVdouble* d2pxy_root2_ = &d2pxy_[root2_];
        VVVdouble* dpxy_root1_  = &dpxy_[root1_];
        VVVdouble* dpxy_root2_  = &dpxy_[root2_];
        VVVdouble* pxy_root1_   = &pxy_[root1_];
        VVVdouble* pxy_root2_   = &pxy_[root2_];
        for (unsigned int i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoodsroot1__i = &(*_likelihoodsroot1_)[(*_patternLinks_fatherroot1_)[i]];
          VVdouble* _likelihoodsroot2__i = &(*_likelihoodsroot2_)[(*_patternLinks_fatherroot2_)[i]];
          VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
          for (unsigned int c = 0; c < nbClasses_; c++)
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
            for (unsigned int x = 0; x < nbStates_; x++)
            {
              Vdouble* d2pxy_root1__c_x = &(*d2pxy_root1__c)[x];
              Vdouble* d2pxy_root2__c_x = &(*d2pxy_root2__c)[x];
              Vdouble* dpxy_root1__c_x  = &(*dpxy_root1__c)[x];
              Vdouble* dpxy_root2__c_x  = &(*dpxy_root2__c)[x];
              Vdouble* pxy_root1__c_x   = &(*pxy_root1__c)[x];
              Vdouble* pxy_root2__c_x   = &(*pxy_root2__c)[x];
              double d2l1 = 0, d2l2 = 0, dl1 = 0, dl2 = 0, l1 = 0, l2 = 0;
              for (unsigned int y = 0; y < nbStates_; y++)
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
        vector<unsigned int> * _patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
        VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

        VVVdouble* pxy__son = &pxy_[son->getId()];
        for (unsigned int i = 0; i < nbSites; i++)
        {
          VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_father_son)[i]];
          VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
          for (unsigned int c = 0; c < nbClasses_; c++)
          {
            Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
            Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
            VVdouble* pxy__son_c = &(*pxy__son)[c];
            for (unsigned int x = 0; x < nbStates_; x++)
            {
              double d2l = 0;
              Vdouble* pxy__son_c_x = &(*pxy__son_c)[x];
              for (unsigned int y = 0; y < nbStates_; y++)
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
     }*/

  // Get the node with the branch whose length must be derivated:
  int brId = atoi(variable.substr(5).c_str());
  const Node* branch = process_->getTree().getNode(brId);
  const Node* father = branch->getFather();

  // Compute dLikelihoods array for the father node.
  // First initialize to 1:
  VVVdouble* d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
  size_t nbSites  = d2Likelihoods_father->size();
  for (size_t i = 0; i < nbSites; ++i)
  {
    VVdouble* d2Likelihoods_father_i = &(*d2Likelihoods_father)[i];
    for (size_t c = 0; c < nbClasses_; ++c)
    {
      Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
      for (size_t s = 0; s < nbStates_; ++s)
      {
        (*d2Likelihoods_father_i_c)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; ++l)
  {
    const Node* son = father->getSon(l);

    vector<size_t>* patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());
    VVVdouble* likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

    if (son == branch)
    {
      // Get all derivatives of transition probabilities:
      vector<const Matrix<double>*> d2pxy_son(nbClasses_);
      for (size_t c = 0; c < nbClasses_; ++c)
      {
        d2pxy_son[c] = &process_->getTransitionProbabilitiesD2(son->getId(), c);
      }

      // Loop over sites:
      for (size_t i = 0; i < nbSites; ++i)
      {
        VVdouble* likelihoods_son_i = &(*likelihoods_son)[(*patternLinks_father_son)[i]];
        VVdouble* d2Likelihoods_father_i = &(*d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; ++c)
        {
          Vdouble* likelihoods_son_i_c = &(*likelihoods_son_i)[c];
          Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
          for (size_t x = 0; x < nbStates_; ++x)
          {
            double d2l = 0;
            for (size_t y = 0; y < nbStates_; ++y)
            {
              d2l += (*d2pxy_son[c])(x, y) * (*likelihoods_son_i_c)[y];
            }
            (*d2Likelihoods_father_i_c)[x] *= d2l;
          }
        }
      }
    }
    else
    {
      // Get all transition probabilities:
      vector<const Matrix<double>*> pxy_son(nbClasses_);
      for (size_t c = 0; c < nbClasses_; ++c)
      {
        pxy_son[c] = &process_->getTransitionProbabilities(son->getId(), c);
      }

      // Loop over sites:
      for (size_t i = 0; i < nbSites; ++i)
      {
        VVdouble* likelihoods_son_i = &(*likelihoods_son)[(*patternLinks_father_son)[i]];
        VVdouble* d2Likelihoods_father_i = &(*d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; ++c)
        {
          Vdouble* likelihoods_son_i_c = &(*likelihoods_son_i)[c];
          Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
          for (size_t x = 0; x < nbStates_; ++x)
          {
            double d2l = 0;
            for (size_t y = 0; y < nbStates_; ++y)
            {
              d2l += (*pxy_son[c])(x, y) * (*likelihoods_son_i_c)[y];
            }
            (*d2Likelihoods_father_i_c)[x] *= d2l;
          }
        }
      }
    }
  }

  // Now we go down the tree toward the root node:
  computeDownSubtreeD2Likelihood_(father);
}

/******************************************************************************/

void SingleRecursiveTreeLikelihoodCalculation::computeDownSubtreeD2Likelihood_(const Node* node)
{
  const Node* father = node->getFather();
  // We assume that the dLikelihoods_ array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if (father == 0)
    return;  // We reached the root!

  // Compute dLikelihoods array for the father node.
  // First initialize to 1:
  VVVdouble* d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());
  size_t nbSites  = d2Likelihoods_father->size();
  for (size_t i = 0; i < nbSites; ++i)
  {
    VVdouble* d2Likelihoods_father_i = &(*d2Likelihoods_father)[i];
    for (size_t c = 0; c < nbClasses_; ++c)
    {
      Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
      for (size_t s = 0; s < nbStates_; ++s)
      {
        (*d2Likelihoods_father_i_c)[s] = 1.;
      }
    }
  }

  size_t nbNodes = father->getNumberOfSons();
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* son = father->getSon(l);

    vector<size_t>* patternLinks_father_son = &likelihoodData_->getArrayPositions(father->getId(), son->getId());

    // Get all transition probabilities:
    vector<const Matrix<double>*> pxy_son(nbClasses_);
    for (size_t c = 0; c < nbClasses_; ++c)
    {
      pxy_son[c] = &process_->getTransitionProbabilities(son->getId(), c);
    }

    if (son == node)
    {
      VVVdouble* d2Likelihoods_son = &likelihoodData_->getD2LikelihoodArray(son->getId());

      // Loop over sites:
      for (size_t i = 0; i < nbSites; i++)
      {
        VVdouble* d2Likelihoods_son_i = &(*d2Likelihoods_son)[(*patternLinks_father_son)[i]];
        VVdouble* d2Likelihoods_father_i = &(*d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; ++c)
        {
          Vdouble* d2Likelihoods_son_i_c = &(*d2Likelihoods_son_i)[c];
          Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
          for (size_t x = 0; x < nbStates_; ++x)
          {
            double d2l = 0;
            for (size_t y = 0; y < nbStates_; ++y)
            {
              d2l += (*pxy_son[c])(x, y) * (*d2Likelihoods_son_i_c)[y];
            }
            (*d2Likelihoods_father_i_c)[x] *= d2l;
          }
        }
      }
    }
    else
    {
      VVVdouble* likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());

      // Loop over sites:
      for (size_t i = 0; i < nbSites; ++i)
      {
        VVdouble* likelihoods_son_i = &(*likelihoods_son)[(*patternLinks_father_son)[i]];
        VVdouble* d2Likelihoods_father_i = &(*d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; ++c)
        {
          Vdouble* likelihoods_son_i_c = &(*likelihoods_son_i)[c];
          Vdouble* d2Likelihoods_father_i_c = &(*d2Likelihoods_father_i)[c];
          for (size_t x = 0; x < nbStates_; ++x)
          {
            double dl = 0;
            for (unsigned int y = 0; y < nbStates_; ++y)
            {
              dl += (*pxy_son[c])(x, y) * (*likelihoods_son_i_c)[y];
            }
            (*d2Likelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
  }

  // Next step: move toward grand father...
  computeDownSubtreeD2Likelihood_(father);
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getD2LikelihoodForASiteForAClass(size_t site, size_t classIndex) const
{
  vector<double> rootFreqs = process_->getRootFrequencies();
  Vdouble* d2la = &likelihoodData_->getD2LikelihoodArray(
           process_->getTree().getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][classIndex];
  double d2l = 0;
  for (size_t i = 0; i < nbStates_; ++i)
  {
    d2l += (*d2la)[i] * rootFreqs[i];
  }
  return d2l;
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getD2LikelihoodForASite(size_t site) const
{
  vector<double> rootFreqs = process_->getRootFrequencies();
  VVdouble* d2la = &likelihoodData_->getD2LikelihoodArray(
           process_->getTree().getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)];
  // Derivative of the sum is the sum of derivatives:
  double d2l = 0;
  for (size_t i = 0; i < nbClasses_; ++i)
  {
    for (size_t j = 0; j < nbStates_; ++j)
    {
      d2l += (*d2la)[i][j] * process_->getProbabilityForModel(i) * rootFreqs[j];
    }
  }
  return d2l;
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getD2LogLikelihoodForASite(size_t site) const
{
  return getD2LikelihoodForASite(site) / getLikelihoodForASite(site)
         - pow( getDLikelihoodForASite(site) / getLikelihoodForASite(site), 2);
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getD2LogLikelihood() const
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t i = 0; i < nbSites_; ++i)
  {
    dl += getD2LogLikelihoodForASite(i);
  }
  return dl;
}

/******************************************************************************/


void SingleRecursiveTreeLikelihoodCalculation::displayLikelihood(const Node* node)
{
  cout << "Likelihoods at node " << node->getName() << ": " << endl;
  displayLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId()));
  cout << "                                         ***" << endl;
}

/******************************************************************************/

