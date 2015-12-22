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
//#include "ComputingNode.h"

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

RecursiveLikelihoodTreeCalculation::RecursiveLikelihoodTreeCalculation(
  const SubstitutionProcess* process,
  bool verbose,
  bool usePatterns)
throw (Exception) :
  AbstractLikelihoodTreeCalculation(process, verbose),
  likelihoodData_(0),
  root1_(-1),
  root2_(-1)
{
  init_(usePatterns);
}

 /******************************************************************************/

RecursiveLikelihoodTreeCalculation::RecursiveLikelihoodTreeCalculation(
  const SiteContainer& data,
  const SubstitutionProcess* process,
  bool verbose,
  bool usePatterns)
throw (Exception) :
  AbstractLikelihoodTreeCalculation(process, verbose),
  likelihoodData_(0),
  root1_(-1),
  root2_(-1)
{
  init_(usePatterns);
  setData(data);
}

/******************************************************************************/

void RecursiveLikelihoodTreeCalculation::init_(bool usePatterns) throw (Exception)
{
  likelihoodData_.reset(new RecursiveLikelihoodTree(
                          *process_,
                          usePatterns));
}

/******************************************************************************/

RecursiveLikelihoodTreeCalculation::RecursiveLikelihoodTreeCalculation(const RecursiveLikelihoodTreeCalculation& tlc) :
  AbstractLikelihoodTreeCalculation(tlc),
  likelihoodData_(0),
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
  Vint upId=process_->getComputingTree().updatedNodes();

  int rootId=process_->getComputingTree()[0]->getRootId();

  if (upId.size()!=0)
  {
    for (size_t c=0; c<nbClasses_; c++){
      for (size_t i=0;i<upId.size();i++){
        RecursiveLikelihoodNode* node=(*likelihoodData_)[c].getNode(upId[i]);
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


/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/

void RecursiveLikelihoodTreeCalculation::computeTreeDLogLikelihood(const string& variable)
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

  /// This is correct only if one branch is derivated
  
  for (size_t c=0; c<nbClasses_; c++){
    RecursiveLikelihoodNode* branch=(*likelihoodData_)[c].getNode(VbrId[0]);
    branch->updateFatherBelow_(false, ComputingNode::D1);
  }
  
  likelihoodData_->computeLikelihoods(process_->getComputingTree(), ComputingNode::D1, &VbrId);

}

/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/

void RecursiveLikelihoodTreeCalculation::computeTreeD2LogLikelihood(const string& variable)
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


  /// This is correct only if one branch is derivated
  
  for (size_t c=0; c<nbClasses_; c++){
    RecursiveLikelihoodNode* branch=(*likelihoodData_)[c].getNode(VbrId[0]);
    branch->updateFatherBelow_(false, ComputingNode::D2);
  }
  
  likelihoodData_->computeLikelihoods(process_->getComputingTree(), ComputingNode::D2, &VbrId);

}



