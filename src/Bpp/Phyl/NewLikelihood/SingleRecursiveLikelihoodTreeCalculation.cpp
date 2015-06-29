//
// File: SingleRecursiveLikelihoodTreeCalculation.cpp
// Created by: Julien Dutheil, Laurent Guéguen
// Created on: mardi 23 juin 2015, à 18h 53
// From file: RNonHomogeneousLikelihoodTree.cpp
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

#include "SingleRecursiveLikelihoodTreeCalculation.h"
#include "SingleRecursiveLikelihoodTree.h"
//#include "ComputingNode.h"

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

SingleRecursiveLikelihoodTreeCalculation::SingleRecursiveLikelihoodTreeCalculation(
  const SubstitutionProcess* process,
  bool verbose,
  bool usePatterns)
throw (Exception) :
  AbstractLikelihoodTreeCalculation(process, verbose),
  likelihoodData_(0),
  root1_(-1),
  root2_(-1),
  nullDLikelihood_(true),
  nullD2Likelihood_(true)
{
  init_(usePatterns);
}

 /******************************************************************************/

SingleRecursiveLikelihoodTreeCalculation::SingleRecursiveLikelihoodTreeCalculation(
  const SiteContainer& data,
  const SubstitutionProcess* process,
  bool verbose,
  bool usePatterns)
throw (Exception) :
  AbstractLikelihoodTreeCalculation(process, verbose),
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

void SingleRecursiveLikelihoodTreeCalculation::init_(bool usePatterns) throw (Exception)
{
  likelihoodData_.reset(new SingleRecursiveLikelihoodTree(
                          *process_,
                          usePatterns));
}

/******************************************************************************/

SingleRecursiveLikelihoodTreeCalculation::SingleRecursiveLikelihoodTreeCalculation(const SingleRecursiveLikelihoodTreeCalculation& tlc) :
  AbstractLikelihoodTreeCalculation(tlc),
  likelihoodData_(0),
  root1_(tlc.root1_),
  root2_(tlc.root2_),
  nullDLikelihood_(tlc.nullDLikelihood_),
  nullD2Likelihood_(tlc.nullD2Likelihood_)
{
  likelihoodData_.reset(tlc.likelihoodData_->clone());
}

/******************************************************************************/

SingleRecursiveLikelihoodTreeCalculation& SingleRecursiveLikelihoodTreeCalculation::operator=(
  const SingleRecursiveLikelihoodTreeCalculation& tlc)
{
  AbstractLikelihoodTreeCalculation::operator=(tlc);
  likelihoodData_.reset(tlc.likelihoodData_->clone());
  root1_ = tlc.root1_;
  root2_ = tlc.root2_;
  nullDLikelihood_ = tlc.nullDLikelihood_;
  nullD2Likelihood_ = tlc.nullD2Likelihood_;
  return *this;
}

/******************************************************************************/

double SingleRecursiveLikelihoodTreeCalculation::getLikelihoodForASite(size_t site)
{
  if (computeLikelihoods_)
    computeTreeLikelihood();
  
  double l = 0;
  size_t posR=likelihoodData_->getRootArrayPosition(site);

  int Rid=process_->getTree().getRootNode()->getId();
  
  for (size_t c = 0; c < nbClasses_; ++c)
  {
    const VVdouble& lla = likelihoodData_->getLikelihoodArray(Rid, c);
    for (size_t j = 0; j < nbStates_; ++j)
    {
      l += lla[posR][j] * process_->getProbabilityForModel(c) * process_->getRootFrequencies()[j];
    }
  }

  if (l < 0) l = 0; //May happen because of numerical errors.
  return l;
}

/******************************************************************************/

double SingleRecursiveLikelihoodTreeCalculation::getLikelihoodForASiteForAClass(size_t site, size_t classIndex)
{
  if (computeLikelihoods_)
  {
    computeTreeLikelihood();
    computeLikelihoods_=false;
  }

  const Vdouble& la = likelihoodData_->getLikelihoodArray(
    process_->getTree().getRootNode()->getId(),classIndex)[likelihoodData_->getRootArrayPosition(site)];

  double l = 0;
  for (size_t i = 0; i < nbStates_; ++i)
  {
    l += la[i] * process_->getRootFrequencies()[i];
  }
  return l;
}

/******************************************************************************/

double SingleRecursiveLikelihoodTreeCalculation::getLikelihoodForASiteForAState(size_t site, int state)
{
  if (computeLikelihoods_)
    computeTreeLikelihood();

  double l = 0;
  int Rid=process_->getTree().getRootNode()->getId();

  size_t posR=likelihoodData_->getRootArrayPosition(site);
  for (size_t c = 0; c < nbClasses_; ++c)
  {
    const VVdouble& lla = likelihoodData_->getLikelihoodArray(Rid,c);
    l += lla[posR][state] * process_->getProbabilityForModel(c);
  }
  return l;
}

/******************************************************************************/

double SingleRecursiveLikelihoodTreeCalculation::getLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state)
{
  if (computeLikelihoods_)
    computeTreeLikelihood();

  return likelihoodData_->getLikelihoodArray(
    process_->getTree().getRootNode()->getId(),classIndex)[likelihoodData_->getRootArrayPosition(site)][state];
}


/******************************************************************************
 *                           Likelihood computation                           *
 ******************************************************************************/

void SingleRecursiveLikelihoodTreeCalculation::computeTreeLikelihood()
{
  computeSubtreeLikelihood_(process_->getTree().getRootNode());
  computeLikelihoods_=false;  
}

/******************************************************************************/

void SingleRecursiveLikelihoodTreeCalculation::computeSubtreeLikelihood_(const Node* node)
{
  if (node->isLeaf())
    return;

  size_t nbNodes  = node->getNumberOfSons();
  vector<const vector<size_t>* > vPatt(nbNodes);

  for (size_t l = 0; l < nbNodes; l++){
    const Node* son = node->getSon(l);    
    computeSubtreeLikelihood_(son);
    vPatt[l]=&likelihoodData_->getArrayPositions(node->getId(), son->getId());
  }

  likelihoodData_->computeUpwardPartialLikelihoods(process_->getComputingTree(), node->getId(), vPatt, ComputingNode::D0);
}


/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/

void SingleRecursiveLikelihoodTreeCalculation::computeTreeDLogLikelihood(const string& variable)
{

  // Get the node with the branch whose length must be derivated:

  Vint VbrId;
  try {
    VbrId.push_back(atoi(variable.substr(5).c_str()));
  }
  catch (std::exception const& e)
  {
    nullDLikelihood_=true;
    return;
  }

  nullDLikelihood_=false;
  const Node* branch = process_->getTree().getNode(VbrId[0]);
  const Node* father = branch->getFather();


  size_t nbNodes = father->getNumberOfSons();

  vector<const vector<size_t>* > vPatt(nbNodes);

  for (size_t l = 0; l < nbNodes; ++l)
  {
    int Sid=father->getSon(l)->getId();
    vPatt[l]=&likelihoodData_->getArrayPositions(father->getId(), Sid);
    likelihoodData_->resetLikelihoods(Sid,ComputingNode::D1);
  }
  
  likelihoodData_->computeUpwardPartialLikelihoods(process_->getComputingTree(), father->getId(), vPatt, ComputingNode::D1, &VbrId);

  // Now we go down the tree toward the root node:
  computeDownSubtreeDLikelihood_(father);
  computeLikelihoodsD1_=false;
}

/******************************************************************************/
      
void SingleRecursiveLikelihoodTreeCalculation::computeDownSubtreeDLikelihood_(const Node* node)
{
  const Node* father = node->getFather();

  if (father == 0)
    return;  // We reached the root!

  size_t nbNodes = father->getNumberOfSons();
  vector<const vector<size_t>* > vPatt;

  for (size_t l = 0; l < nbNodes; ++l)
  {
    int Sid=father->getSon(l)->getId();
    vPatt[l]=&likelihoodData_->getArrayPositions(father->getId(), Sid);
  }

  likelihoodData_->computeUpwardPartialLikelihoods(process_->getComputingTree(), father->getId(), vPatt, ComputingNode::D1);
     
  // Next step: move toward grand father...
  computeDownSubtreeDLikelihood_(father);
}

/******************************************************************************/

double SingleRecursiveLikelihoodTreeCalculation::getDLikelihoodForASite(size_t site)
{
  if (nullDLikelihood_)
    return 0;

  if (computeLikelihoodsD1_){
    throw Exception("getDLogLikelihoodForASite() : DLogLikelihood not computed");
    computeLikelihoodsD1_=false;
  }


  size_t posR=likelihoodData_->getRootArrayPosition(site);
  int Rid=process_->getTree().getRootNode()->getId();
  
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t c = 0; c < nbClasses_; c++)
  {
    VVdouble* ldla = &likelihoodData_->getDLikelihoodArray(Rid, c);
    for (size_t j = 0; j < nbStates_; ++j)
      dl += (*ldla)[posR][j] * process_->getProbabilityForModel(c) * process_->getRootFrequencies()[j];
  }
  
  return dl;
}

/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/

void SingleRecursiveLikelihoodTreeCalculation::computeTreeD2LogLikelihood(const string& variable)
{

  // Get the node with the branch whose length must be derivated:
  Vint VbrId;
  
  try {
    VbrId.push_back(atoi(variable.substr(5).c_str()));
  }
  catch (std::exception const& e)
  {
    nullD2Likelihood_=true;
    return;
  }

  nullD2Likelihood_=false;

  const Node* branch = process_->getTree().getNode(VbrId[0]);
  const Node* father = branch->getFather();

  size_t nbNodes = father->getNumberOfSons();
 
  vector<const vector<size_t>* > vPatt(nbNodes);
  
  for (size_t l = 0; l < nbNodes; ++l)
  {
    int Sid=father->getSon(l)->getId();
    vPatt[l]=&likelihoodData_->getArrayPositions(father->getId(), Sid);
    likelihoodData_->resetLikelihoods(Sid,ComputingNode::D2);
  }

  likelihoodData_->computeUpwardPartialLikelihoods(process_->getComputingTree(), father->getId(), vPatt, ComputingNode::D2, &VbrId);

  // Now we go down the tree toward the root node:
  computeDownSubtreeD2Likelihood_(father);
  computeLikelihoodsD2_=false;
}

/******************************************************************************/

void SingleRecursiveLikelihoodTreeCalculation::computeDownSubtreeD2Likelihood_(const Node* node)
{
  const Node* father = node->getFather();
  // We assume that the dLikelihoods_ array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if (father == 0)
    return;  // We reached the root!
  
  size_t nbNodes = father->getNumberOfSons();
  vector<const vector<size_t>* > vPatt;
  
  for (size_t l = 0; l < nbNodes; ++l)
  {
    int Sid=father->getSon(l)->getId();
    vPatt[l]=&likelihoodData_->getArrayPositions(father->getId(), Sid);
  }
  
  likelihoodData_->computeUpwardPartialLikelihoods(process_->getComputingTree(), father->getId(), vPatt, ComputingNode::D2);
 
  // Next step: move toward grand father...
  computeDownSubtreeD2Likelihood_(father);
}

/******************************************************************************/

double SingleRecursiveLikelihoodTreeCalculation::getD2LikelihoodForASite(size_t site)
{
  if (nullD2Likelihood_)
    return 0;

  if (computeLikelihoodsD2_){
    throw Exception("getD2LogLikelihoodForASite() : D2LogLikelihood not computed");
    computeLikelihoodsD2_=false;
  }

  size_t posR=likelihoodData_->getRootArrayPosition(site);
  int Rid=process_->getTree().getRootNode()->getId();

  // Derivative of the sum is the sum of derivatives:
  double d2l = 0;
  for (size_t c = 0; c < nbClasses_; c++)
  {
    VVdouble* d2la = &likelihoodData_->getD2LikelihoodArray(Rid, c);
    for (size_t j = 0; j < nbStates_; ++j)
      d2l += (*d2la)[posR][j] * process_->getProbabilityForModel(c) * process_->getRootFrequencies()[j];
  }

  return d2l;
}

/******************************************************************************/

