//
// File: RHomogeneousClockTreeLikelihood.cpp
// Created by: Beno�t Nabholz
//             Julien Dutheil
// Created on: Fri Apr 06 14:11 2007
//

/*
Copyright or � or Copr. CNRS, (November 16, 2004)

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

#include "RHomogeneousClockTreeLikelihood.h"
#include "TreeTemplateTools.h"

#include <iostream>

using namespace std;

// From Utils:
#include <Utils/ApplicationTools.h>

using namespace bpp;

/******************************************************************************/

RHomogeneousClockTreeLikelihood::RHomogeneousClockTreeLikelihood(
  const Tree & tree,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose)
throw (Exception):
  RHomogeneousTreeLikelihood(tree, model, rDist, false, verbose, true)
{
  _init();
}

/******************************************************************************/

RHomogeneousClockTreeLikelihood::RHomogeneousClockTreeLikelihood(
  const Tree & tree,
  const SiteContainer & data,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose)
throw (Exception):
  RHomogeneousTreeLikelihood(tree, data, model, rDist, false, verbose, true)
{
  _init();
}

/******************************************************************************/

void RHomogeneousClockTreeLikelihood::_init()
{
  //Check if the tree is rooted:
  if(!_tree->isRooted()) throw Exception("RHomogeneousClockTreeLikelihood::init(). Tree is unrooted!");
  if(TreeTemplateTools::isMultifurcating(*_tree->getRootNode())) throw Exception("HomogeneousClockTreeLikelihood::init(). Tree is multifurcating.");
  setMinimumBranchLength(0.);
}

/******************************************************************************/

void RHomogeneousClockTreeLikelihood::applyParameters() throw (Exception)
{
  if(!_initialized) throw Exception("RHomogeneousClockTreeLikelihood::applyParameters(). Object not initialized.");
   //Apply branch lengths:
  brLenParameters_.matchParametersValues(getParameters());
  computeBranchLengthsFromHeights(_tree->getRootNode(), brLenParameters_.getParameter("TotalHeight").getValue());

  //Apply substitution model parameters:
  model_->matchParametersValues(getParameters());
  //Apply rate distribution parameters:
  _rateDistribution->matchParametersValues(getParameters());
}

/******************************************************************************/

void RHomogeneousClockTreeLikelihood::fireParameterChanged(const ParameterList & params)
{
  applyParameters();

  computeAllTransitionProbabilities();

  computeTreeLikelihood();
}

/******************************************************************************/

void RHomogeneousClockTreeLikelihood::initBranchLengthsParameters()
{
  //Check branch lengths first:
  for(unsigned int i = 0; i < nbNodes_; i++)
  {
    double d = minimumBrLen_;
    if(!nodes_[i]->hasDistanceToFather())
    {
      ApplicationTools::displayWarning("Missing branch length " + TextTools::toString(i) + ". Value is set to " + TextTools::toString(minimumBrLen_));
      nodes_[i]->setDistanceToFather(minimumBrLen_);
    }
    else
    {
      d = nodes_[i]->getDistanceToFather();
      if (d < minimumBrLen_)
      {
        ApplicationTools::displayWarning("Branch length " + TextTools::toString(i) + " is too small: " + TextTools::toString(d) + ". Value is set to " + TextTools::toString(minimumBrLen_));
        nodes_[i]->setDistanceToFather(minimumBrLen_);
      }
    }
  }

  brLenParameters_.reset();

  map<const Node *, double> heights;
  TreeTemplateTools::getHeights(*_tree->getRootNode(), heights);
  double totalHeight = heights[_tree->getRootNode()];
  brLenParameters_.addParameter(Parameter("TotalHeight", totalHeight, brLenConstraint_->clone(), true)); 
  for(map<const Node *, double>::iterator it = heights.begin(); it != heights.end(); it++)
  {
    if(!it->first->isLeaf() && it->first->hasFather())
    {
      double fatherHeight = heights[it->first->getFather()];
      brLenParameters_.addParameter(Parameter("HeightP" + TextTools::toString(it->first->getId()), it->second / fatherHeight, &Parameter::PROP_CONSTRAINT_IN));
    }
  }
}

/******************************************************************************/

void RHomogeneousClockTreeLikelihood::computeBranchLengthsFromHeights(Node * node, double height) throw (Exception)
{
  for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    Node * son = node->getSon(i);
    if(son->isLeaf())
    {
      son->setDistanceToFather(std::max(minimumBrLen_, height));
    }
    else
    {
      Parameter * p = &brLenParameters_.getParameter(string("HeightP") + TextTools::toString(son->getId()));
      double sonHeightP = p->getValue();
      double sonHeight = sonHeightP * height;
      son->setDistanceToFather(std::max(minimumBrLen_, height - sonHeight));
      computeBranchLengthsFromHeights(son, sonHeight);
    }
  }
}

/******************************************************************************/

ParameterList RHomogeneousClockTreeLikelihood::getDerivableParameters() const throw (Exception)
{
  if(!_initialized) throw Exception("RHomogeneousClockTreeLikelihood::getDerivableParameters(). Object is not initialized.");
  return ParameterList();
}

/******************************************************************************/

ParameterList RHomogeneousClockTreeLikelihood::getNonDerivableParameters() const throw (Exception)
{
  if(!_initialized) throw Exception("RHomogeneousClockTreeLikelihood::getNonDerivableParameters(). Object is not initialized.");
  return getParameters();
}

/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/  

double RHomogeneousClockTreeLikelihood::getFirstOrderDerivative(const string & variable) const
throw (Exception)
{ 
  throw Exception("No first order derivative is implemented for this function.");
}

/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/  

double RHomogeneousClockTreeLikelihood::getSecondOrderDerivative(const string & variable) const
throw (Exception)
{
  throw Exception("No second order derivative is implemented for this function.");
}

/******************************************************************************/

