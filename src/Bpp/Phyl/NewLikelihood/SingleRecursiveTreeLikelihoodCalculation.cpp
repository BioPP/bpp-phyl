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
#include "ComputingNode.h"

using namespace bpp;
using namespace newlik;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

SingleRecursiveTreeLikelihoodCalculation::SingleRecursiveTreeLikelihoodCalculation(
  const SubstitutionProcess* process,
  bool verbose,
  bool usePatterns)
throw (Exception) :
  AbstractTreeLikelihoodCalculation(process, verbose),
  likelihoodData_(0),
  root1_(-1),
  root2_(-1),
  nullDLikelihood_(true),
  nullD2Likelihood_(true)
{
  init_(usePatterns);
}

 /******************************************************************************/

SingleRecursiveTreeLikelihoodCalculation::SingleRecursiveTreeLikelihoodCalculation(
  const SiteContainer& data,
  const SubstitutionProcess* process,
  bool verbose,
  bool usePatterns)
throw (Exception) :
  AbstractTreeLikelihoodCalculation(process, verbose),
  likelihoodData_(0),
  root1_(-1),
  root2_(-1),
  nullDLikelihood_(true),
  nullD2Likelihood_(true)
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
  root2_(tlc.root2_),
  nullDLikelihood_(tlc.nullDLikelihood_),
  nullD2Likelihood_(tlc.nullD2Likelihood_)
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
  nullDLikelihood_ = tlc.nullDLikelihood_;
  nullD2Likelihood_ = tlc.nullD2Likelihood_;
  return *this;
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getLikelihoodForASite(size_t site)
{
  if (computeLikelihoods_)
    computeTreeLikelihood();
  
  double l = 0;
  VVVdouble* lla = &likelihoodData_->getLikelihoodArray(process_->getTree().getRootNode()->getId());
  size_t posR=likelihoodData_->getRootArrayPosition(site);

  
  for (size_t c = 0; c < nbClasses_; ++c)
  {
    for (size_t j = 0; j < nbStates_; ++j)
    {
      l += (*lla)[c][posR][j] * process_->getProbabilityForModel(c) * process_->getRootFrequencies()[j];
    }
  }

  if (l < 0) l = 0; //May happen because of numerical errors.
  return l;
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getLikelihoodForASiteForAClass(size_t site, size_t classIndex)
{
  if (computeLikelihoods_)
  {
    computeTreeLikelihood();
    computeLikelihoods_=false;
  }

  Vdouble* la = &likelihoodData_->getLikelihoodArray(
    process_->getTree().getRootNode()->getId())[classIndex][likelihoodData_->getRootArrayPosition(site)];

  double l = 0;
  for (size_t i = 0; i < nbStates_; ++i)
  {
    l += (*la)[i] * process_->getRootFrequencies()[i];
  }
  return l;
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getLikelihoodForASiteForAState(size_t site, int state)
{
  if (computeLikelihoods_)
    computeTreeLikelihood();

  double l = 0;
  VVVdouble* lla = &likelihoodData_->getLikelihoodArray(process_->getTree().getRootNode()->getId());
  size_t posR=likelihoodData_->getRootArrayPosition(site);
  for (size_t c = 0; c < nbClasses_; ++c)
  {
    l += (*lla)[c][posR][state] * process_->getProbabilityForModel(c);
  }
  return l;
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state)
{
  if (computeLikelihoods_)
    computeTreeLikelihood();

  return likelihoodData_->getLikelihoodArray(
    process_->getTree().getRootNode()->getId())[classIndex][likelihoodData_->getRootArrayPosition(site)][state];
}


/******************************************************************************
 *                           Likelihood computation                           *
 ******************************************************************************/

void SingleRecursiveTreeLikelihoodCalculation::computeTreeLikelihood()
{
  computeSubtreeLikelihood_(process_->getTree().getRootNode());
  computeLikelihoods_=false;  
}

/******************************************************************************/

void SingleRecursiveTreeLikelihoodCalculation::computeSubtreeLikelihood_(const Node* node)
{
  if (node->isLeaf())
    return;

  size_t nbNodes  = node->getNumberOfSons();

  // Must reset the likelihood array first (i.e. set all of them to 1):
  TreeLikelihoodData::resetLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId()));

  VVVdouble* likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());

  vector<const vector<size_t>* > vPatt;
  vector<const VVVdouble*> vLikArr;
  
  for (size_t l = 0; l < nbNodes; l++){
    const Node* son = node->getSon(l);
    
    computeSubtreeLikelihood_(son);
    vPatt.push_back(&likelihoodData_->getArrayPositions(node->getId(), son->getId()));
    vLikArr.push_back(&likelihoodData_->getLikelihoodArray(son->getId()));
  }

  process_->multiplyUpwardPartialLikelihoods(likelihoods_node, vLikArr, node->getId(), vPatt, ComputingNode::D0);
}


/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/

void SingleRecursiveTreeLikelihoodCalculation::computeTreeDLogLikelihood(const string& variable)
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

  int brId;
  try {
    brId = atoi(variable.substr(5).c_str());
  }
  catch (std::exception const& e)
  {
    nullDLikelihood_=true;
    return;
  }

  nullDLikelihood_=false;
  const Node* branch = process_->getTree().getNode(brId);
  const Node* father = branch->getFather();

  // Compute dLikelihoods array for the father node.
  // First initialize to 1:
  TreeLikelihoodData::resetLikelihoodArray(likelihoodData_->getDLikelihoodArray(father->getId()));

  VVVdouble* dLikelihoods_father = &(likelihoodData_->getDLikelihoodArray(father->getId()));
  size_t nbNodes = father->getNumberOfSons();

  vector<const vector<size_t>* > vPatt(nbNodes);
  vector<const VVVdouble*> vLikArr(nbNodes);

  size_t nBrId=0;
  
  for (size_t l = 0; l < nbNodes; ++l)
  {
    int Sid=father->getSon(l)->getId();
    vPatt[l]=&likelihoodData_->getArrayPositions(father->getId(), Sid);
    
    if (Sid!=brId)
      vLikArr[l]=&likelihoodData_->getLikelihoodArray(Sid);
    else{
      vLikArr[l]=0;
      nBrId=l;
    }
    
  }

  process_->multiplyUpwardPartialLikelihoods(dLikelihoods_father, vLikArr, father->getId(), vPatt, ComputingNode::D0);

  process_->multiplyUpwardPartialLikelihoods(dLikelihoods_father, &(likelihoodData_->getLikelihoodArray(brId)), brId, *vPatt[nBrId], ComputingNode::D1);
  
  // Now we go down the tree toward the root node:
  computeDownSubtreeDLikelihood_(father);
  computeLikelihoodsD1_=false;
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
  TreeLikelihoodData::resetLikelihoodArray(likelihoodData_->getDLikelihoodArray(father->getId()));

  VVVdouble* dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(father->getId());

  int brId=node->getId();

  size_t nbNodes = father->getNumberOfSons();

  vector<const vector<size_t>* > vPatt;
  vector<const VVVdouble*> vLikArr;

  for (size_t l = 0; l < nbNodes; ++l)
  {
    int Sid=father->getSon(l)->getId();
    
    vPatt.push_back(&likelihoodData_->getArrayPositions(father->getId(), Sid));
    if (Sid==brId)
      vLikArr.push_back(&likelihoodData_->getDLikelihoodArray(Sid));
    else
      vLikArr.push_back(&likelihoodData_->getLikelihoodArray(Sid));
  }

  process_->multiplyUpwardPartialLikelihoods(dLikelihoods_father, vLikArr, father->getId(), vPatt, ComputingNode::D0);
     
  // Next step: move toward grand father...
  computeDownSubtreeDLikelihood_(father);
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getDLikelihoodForASite(size_t site)
{
  if (nullDLikelihood_)
    return 0;

  if (computeLikelihoodsD1_){
    throw Exception("getDLogLikelihoodForASite() : DLogLikelihood not computed");
    computeLikelihoodsD1_=false;
  }

  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  VVVdouble* ldla = &likelihoodData_->getDLikelihoodArray(process_->getTree().getRootNode()->getId());
  size_t posR=likelihoodData_->getRootArrayPosition(site);

  for (size_t c = 0; c < nbClasses_; c++)
  {
    for (size_t j = 0; j < nbStates_; ++j)
    {
      dl += (*ldla)[c][posR][j] * process_->getProbabilityForModel(c) * process_->getRootFrequencies()[j];
    }
  }
  return dl;
}

/******************************************************************************/

// double SingleRecursiveTreeLikelihoodCalculation::getDLikelihoodForASiteForAClass(size_t site, size_t classIndex) const
// {
//   if (nullDLikelihood_)
//     return 0;
  
//   Vdouble* dla = &likelihoodData_->getDLikelihoodArray(
//     process_->getTree().getRootNode()->getId())[classIndex][likelihoodData_->getRootArrayPosition(site)];
 
//   double dl = 0;
//   for (size_t i = 0; i < nbStates_; ++i)
//   {
//     dl += (*dla)[i] * process_->getRootFrequencies()[i];
//   }
//   return dl;
// }


/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/

void SingleRecursiveTreeLikelihoodCalculation::computeTreeD2LogLikelihood(const string& variable)
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
  int brId;
  
  try {
    brId = atoi(variable.substr(5).c_str());
  }
  catch (std::exception const& e)
  {
    nullD2Likelihood_=true;
    return;
  }

  nullD2Likelihood_=false;

  const Node* branch = process_->getTree().getNode(brId);
  const Node* father = branch->getFather();

  // Compute dLikelihoods array for the father node.
  // First initialize to 1:
  TreeLikelihoodData::resetLikelihoodArray(likelihoodData_->getD2LikelihoodArray(father->getId()));

  VVVdouble* d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());

  size_t nbNodes = father->getNumberOfSons();
 
  vector<const vector<size_t>* > vPatt(nbNodes);
  vector<const VVVdouble*> vLikArr(nbNodes);

  size_t nBrId=0;
  
  for (size_t l = 0; l < nbNodes; ++l)
  {
    int Sid=father->getSon(l)->getId();
    vPatt[l]=&likelihoodData_->getArrayPositions(father->getId(), Sid);
    
    if (Sid!=brId)
      vLikArr[l]=&likelihoodData_->getLikelihoodArray(Sid);
    else{
      vLikArr[l]=0;
      nBrId=l;
    }
    
  }

  process_->multiplyUpwardPartialLikelihoods(d2Likelihoods_father, vLikArr, father->getId(), vPatt, ComputingNode::D0);

  process_->multiplyUpwardPartialLikelihoods(d2Likelihoods_father, &(likelihoodData_->getLikelihoodArray(brId)), brId, *vPatt[nBrId], ComputingNode::D2);
  
  // Now we go down the tree toward the root node:
  computeDownSubtreeD2Likelihood_(father);

  computeLikelihoodsD2_=false;
}

/******************************************************************************/

void SingleRecursiveTreeLikelihoodCalculation::computeDownSubtreeD2Likelihood_(const Node* node)
{
  const Node* father = node->getFather();
  // We assume that the dLikelihoods_ array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if (father == 0)
    return;  // We reached the root!

  int brId=node->getId();

  // Compute dLikelihoods array for the father node.
  // First initialize to 1:
  TreeLikelihoodData::resetLikelihoodArray(likelihoodData_->getD2LikelihoodArray(father->getId()));

  VVVdouble* d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(father->getId());

  size_t nbNodes = father->getNumberOfSons();
  vector<const vector<size_t>* > vPatt;
  vector<const VVVdouble*> vLikArr;
  
  for (size_t l = 0; l < nbNodes; ++l)
  {
    int Sid=father->getSon(l)->getId();

    vPatt.push_back(&likelihoodData_->getArrayPositions(father->getId(), Sid));
    if (Sid==brId)
      vLikArr.push_back(&likelihoodData_->getD2LikelihoodArray(Sid));
    else
      vLikArr.push_back(&likelihoodData_->getLikelihoodArray(Sid));
  }
  
  process_->multiplyUpwardPartialLikelihoods(d2Likelihoods_father, vLikArr, father->getId(), vPatt, ComputingNode::D0);
 
  // Next step: move toward grand father...
  computeDownSubtreeD2Likelihood_(father);
}

/******************************************************************************/

double SingleRecursiveTreeLikelihoodCalculation::getD2LikelihoodForASite(size_t site)
{
  if (nullD2Likelihood_)
    return 0;

  if (computeLikelihoodsD2_){
    throw Exception("getD2LogLikelihoodForASite() : D2LogLikelihood not computed");
    computeLikelihoodsD2_=false;
  }

  VVVdouble* d2la = &likelihoodData_->getD2LikelihoodArray(process_->getTree().getRootNode()->getId());
  size_t posR=likelihoodData_->getRootArrayPosition(site);

  // Derivative of the sum is the sum of derivatives:
  double d2l = 0;
  for (size_t c = 0; c < nbClasses_; c++)
  {
    for (size_t j = 0; j < nbStates_; ++j)
    {
      d2l += (*d2la)[c][posR][j] * process_->getProbabilityForModel(c) * process_->getRootFrequencies()[j];
    }
  }
  return d2l;
}

/******************************************************************************/

// double SingleRecursiveTreeLikelihoodCalculation::getD2LikelihoodForASiteForAClass(size_t site, size_t classIndex) const
// {
//   if (nullD2Likelihood_)
//     return 0;
  
//   Vdouble* d2la = &likelihoodData_->getD2LikelihoodArray(
//     process_->getTree().getRootNode()->getId())[classIndex][likelihoodData_->getRootArrayPosition(site)];
//   double d2l = 0;
//   for (size_t i = 0; i < nbStates_; ++i)
//   {
//     d2l += (*d2la)[i] * process_->getRootFrequencies()[i];
//   }
//   return d2l;
// }

/******************************************************************************/


void SingleRecursiveTreeLikelihoodCalculation::displayLikelihood(const Node* node)
{
  cout << "Likelihoods at node " << node->getName() << ": " << endl;
  displayLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId()));
  cout << "                                         ***" << endl;
}

/******************************************************************************/

