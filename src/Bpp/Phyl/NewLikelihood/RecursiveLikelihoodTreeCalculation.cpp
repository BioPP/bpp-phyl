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
  root2_(-1),
  nullDLikelihood_(true),
  nullD2Likelihood_(true)
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
  root2_(-1),
  nullDLikelihood_(true),
  nullD2Likelihood_(true)
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
  root2_(tlc.root2_),
  nullDLikelihood_(tlc.nullDLikelihood_),
  nullD2Likelihood_(tlc.nullD2Likelihood_)
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
  nullDLikelihood_ = tlc.nullDLikelihood_;
  nullD2Likelihood_ = tlc.nullD2Likelihood_;
  return *this;
}

/******************************************************************************/

double RecursiveLikelihoodTreeCalculation::getLikelihoodForASite(size_t site)
{
  double l = 0;
  size_t posR=getLikelihoodData().getRootArrayPosition(site);
  int Rid=process_->getTree().getRootNode()->getId();
  
  for (size_t c = 0; c < nbClasses_; ++c)
  {
    const VVdouble& lla = likelihoodData_->getLikelihoodArray(Rid, c, ComputingNode::D0);

    for (size_t j = 0; j < nbStates_; ++j)
    {
      l += lla[posR][j] * process_->getProbabilityForModel(c);
    }
  }

  if (l < 0) l = 0; //May happen because of numerical errors.
  return l;
}

/******************************************************************************/

double RecursiveLikelihoodTreeCalculation::getLikelihoodForASiteForAClass(size_t site, size_t classIndex)
{
  const Vdouble& la = likelihoodData_->getLikelihoodArray(
    process_->getTree().getRootNode()->getId(),classIndex, ComputingNode::D0)[likelihoodData_->getRootArrayPosition(site)];

  return VectorTools::sum(la);
}

/******************************************************************************/

double RecursiveLikelihoodTreeCalculation::getLikelihoodForASiteForAState(size_t site, int state)
{
  double l = 0;
  int Rid=process_->getTree().getRootNode()->getId();

  size_t posR=likelihoodData_->getRootArrayPosition(site);
  for (size_t c = 0; c < nbClasses_; ++c)
  {
    const VVdouble& lla = likelihoodData_->getLikelihoodArray(Rid,c, ComputingNode::D0);
    l += lla[posR][state] * process_->getProbabilityForModel(c);
  }
  return l;
}

/******************************************************************************/

double RecursiveLikelihoodTreeCalculation::getLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state)
{
  return likelihoodData_->getLikelihoodArray(
    process_->getTree().getRootNode()->getId(),classIndex, ComputingNode::D0)[likelihoodData_->getRootArrayPosition(site)][state];
}


/******************************************************************************
 *                           Likelihood computation                           *
 ******************************************************************************/

bool RecursiveLikelihoodTreeCalculation::updateLikelihoodFlags_()
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
          node->updateFatherBelow(false, ComputingNode::D0);
          node->updateAbove(false);
        }
        else{
          node->updateAbove(false);
          node->setAboveLikelihoods(process_->getRootFrequencies());
          node->updateAbove(true);
        }
      }
    }
    return true;
  }
  return false;
}

/******************************************************************************/

void RecursiveLikelihoodTreeCalculation::computeTreeLikelihood()
{
  if (updateLikelihoodFlags_())
    likelihoodData_->computeLikelihoods(process_->getComputingTree(), ComputingNode::D0);
}


/******************************************************************************/

void RecursiveLikelihoodTreeCalculation::computeLikelihoodsAtNode(int nodeId)
{
  if (likelihoodData_->usePatterns())
    throw Exception("RecursiveLikelihoodTreeCalculation::computeLikelihoodsAtNode not available wth patterns.");
  
  if (!likelihoodData_->isAboveLikelihoodsInitialized())
  {
    likelihoodData_->resetInnerAboveLikelihoods();
    
    likelihoodData_->resetInnerLikelihoods(nbDistinctSites_, nbStates_, ComputingNode::D0);
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
    branch->updateFatherBelow(false, ComputingNode::D1);
  }
  
  likelihoodData_->computeLikelihoods(process_->getComputingTree(), ComputingNode::D1, &VbrId);

}

/******************************************************************************/

double RecursiveLikelihoodTreeCalculation::getDLikelihoodForASite(size_t site)
{
  if (nullDLikelihood_)
    return 0;

  size_t posR=likelihoodData_->getRootArrayPosition(site);
  int Rid=process_->getTree().getRootNode()->getId();
  
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t c = 0; c < nbClasses_; c++)
  {
    VVdouble& ldla = likelihoodData_->getLikelihoodArray(Rid, c, ComputingNode::D1);
    for (size_t j = 0; j < nbStates_; ++j)
      dl += ldla[posR][j] * process_->getProbabilityForModel(c);
  }
  
  return dl;
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
    branch->updateFatherBelow(false, ComputingNode::D2);
  }
  
  likelihoodData_->computeLikelihoods(process_->getComputingTree(), ComputingNode::D2, &VbrId);

}

/******************************************************************************/

double RecursiveLikelihoodTreeCalculation::getD2LikelihoodForASite(size_t site)
{
  if (nullD2Likelihood_)
    return 0;

  size_t posR=likelihoodData_->getRootArrayPosition(site);
  int Rid=process_->getTree().getRootNode()->getId();

  // Derivative of the sum is the sum of derivatives:
  double d2l = 0;
  for (size_t c = 0; c < nbClasses_; c++)
  {
    VVdouble& d2la = likelihoodData_->getLikelihoodArray(Rid, c, ComputingNode::D2);
    for (size_t j = 0; j < nbStates_; ++j)
      d2l += d2la[posR][j] * process_->getProbabilityForModel(c);
  }

  return d2l;
}


