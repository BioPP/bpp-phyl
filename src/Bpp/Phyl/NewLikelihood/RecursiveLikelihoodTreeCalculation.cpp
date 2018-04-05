//
// File: RecursiveLikelihoodTreeCalculation.cpp
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

#include "RecursiveLikelihoodTreeCalculation.h"
#include "RecursiveLikelihoodTree.h"

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

RecursiveLikelihoodTreeCalculation::RecursiveLikelihoodTreeCalculation(
  const SubstitutionProcess* process,
  bool verbose,
  bool usePatterns) :
  AbstractLikelihoodTreeCalculation(process, verbose),
  likelihoodData_(),
  root1_(-1),
  root2_(-1)
{
  init_(usePatterns);
}

 /******************************************************************************/

RecursiveLikelihoodTreeCalculation::RecursiveLikelihoodTreeCalculation(
  const AlignedValuesContainer& data,
  const SubstitutionProcess* process,
  bool verbose,
  bool usePatterns) :
  AbstractLikelihoodTreeCalculation(process, verbose),
  likelihoodData_(),
  root1_(-1),
  root2_(-1)
{
  init_(usePatterns);
  setData(data);
}

/******************************************************************************/

void RecursiveLikelihoodTreeCalculation::init_(bool usePatterns) 
{
  likelihoodData_.reset(new RecursiveLikelihoodTree(
                          *process_,
                          usePatterns));
}

/******************************************************************************/

RecursiveLikelihoodTreeCalculation::RecursiveLikelihoodTreeCalculation(const RecursiveLikelihoodTreeCalculation& tlc) :
  AbstractLikelihoodTreeCalculation(tlc),
  likelihoodData_(),
  root1_(tlc.root1_),
  root2_(tlc.root2_)
{
  likelihoodData_.reset(tlc.likelihoodData_->clone());
}

/******************************************************************************/

RecursiveLikelihoodTreeCalculation& RecursiveLikelihoodTreeCalculation::operator=(
  const RecursiveLikelihoodTreeCalculation& tlc)
{
  AbstractLikelihoodTreeCalculation::operator=(tlc);
  likelihoodData_.reset(tlc.likelihoodData_->clone());
  root1_ = tlc.root1_;
  root2_ = tlc.root2_;
  return *this;
}

/******************************************************************************
 *                           Likelihood computation                           *
 ******************************************************************************/

void RecursiveLikelihoodTreeCalculation::updateLikelihoodFlags_()
{
  Vuint upId=process_->getComputingTree().toBeUpdatedNodes();

  unsigned int rootId=process_->getComputingTree()[0]->getNodeIndex(process_->getComputingTree()[0]->getRoot());

  if (upId.size()!=0)
  {
    for (size_t c=0; c<nbClasses_; c++){
      for (size_t i=0;i<upId.size();i++){
        shared_ptr<RecursiveLikelihoodNode> node=(*likelihoodData_)[c].getNode(upId[i]);
        if (upId[i]!=rootId)
        {
          node->updateFatherBelow_(false, ComputingNode::D0);
          node->updateAbove(false);
        }
        else{
          const Vdouble& rf=process_->getRootFrequencies();
          
          if ((!node->usesLog() && node->getAboveLikelihoodArray_()[0]!=rf)
              || (node->usesLog() && node->getAboveLikelihoodArray_()[0]!=VectorTools::log(rf)))
          {
            node->updateAbove(false);
            node->setAboveLikelihoods(process_->getRootFrequencies());
            node->updateAbove(true);
          }
        }
      }
    }
    up2date_=false;
  }
}

/******************************************************************************/

void RecursiveLikelihoodTreeCalculation::computeTreeLikelihood()
{
  if (!up2date_){
    likelihoodData_->computeLikelihoods(process_->getComputingTree(), ComputingNode::D0);
    up2date_=true;
  }
}


/******************************************************************************/

void RecursiveLikelihoodTreeCalculation::computeLikelihoodsAtNode(int nodeId)
{
  if (likelihoodData_->usePatterns())
    throw Exception("RecursiveLikelihoodTreeCalculation::computeLikelihoodsAtNode not available wth patterns.");
  
  if (!likelihoodData_->isAboveLikelihoodsInitialized())
  {
    likelihoodData_->resetInnerAboveLikelihoods();
    likelihoodData_->resetDownwardLikelihoods(nbDistinctSites_, nbStates_, ComputingNode::D0);
  }

  // recursion

  likelihoodData_->computeLikelihoodsAtNode(process_->getComputingTree(), nodeId);
}


/******************************************************************************/

void RecursiveLikelihoodTreeCalculation::computeLikelihoodsAtAllNodes()
{
  if (likelihoodData_->usePatterns())
    throw Exception("RecursiveLikelihoodTreeCalculation::computeLikelihoodsAtAllNodes not available wth patterns.");
  
  if (!likelihoodData_->isAboveLikelihoodsInitialized())
  {
    likelihoodData_->resetInnerAboveLikelihoods();
    likelihoodData_->resetDownwardLikelihoods(nbDistinctSites_, nbStates_, ComputingNode::D0);
  }

  // recursion

  vector<uint> nodesId=process_->getParametrizablePhyloTree().getAllNodesIndexes();

  for (auto nodeId: nodesId)
    likelihoodData_->computeLikelihoodsAtNode(process_->getComputingTree(), nodeId);
}


/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/

void RecursiveLikelihoodTreeCalculation::computeTreeDLogLikelihood(const Vuint& vbrId)
{
  if (vbrId.size()==0)
  {
    nullDLogLikelihood_=true;
    return;
  }

  nullDLogLikelihood_=false;
  
  vector<unsigned int> lId=(*process_->getComputingTree()[0]).getNodeIndexes((*process_->getComputingTree()[0]).getAllLeaves());
  
  for (size_t i=0;i<lId.size();i++)
    for (size_t c=0; c<nbClasses_; c++){
      shared_ptr<RecursiveLikelihoodNode> branch=(*likelihoodData_)[c].getNode(lId[i]);
      branch->updateFatherBelow_(false, ComputingNode::D1);
    }
  
  likelihoodData_->computeLikelihoods(process_->getComputingTree(), ComputingNode::D1, &vbrId);
}

/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/

void RecursiveLikelihoodTreeCalculation::computeTreeD2LogLikelihood(const Vuint& vbrId)
{
  if (vbrId.size()==0)
  {
    nullD2LogLikelihood_=true;
    return;
  }

  nullD2LogLikelihood_=false;

  vector<unsigned int> lId=(*process_->getComputingTree()[0]).getNodeIndexes((*process_->getComputingTree()[0]).getAllLeaves());
  
  for (size_t i=0;i<lId.size();i++)
    for (size_t c=0; c<nbClasses_; c++){
      shared_ptr<RecursiveLikelihoodNode> branch=(*likelihoodData_)[c].getNode(lId[i]);
      branch->updateFatherBelow_(false, ComputingNode::D2);
    }
  
  likelihoodData_->computeLikelihoods(process_->getComputingTree(), ComputingNode::D2, &vbrId);
}



