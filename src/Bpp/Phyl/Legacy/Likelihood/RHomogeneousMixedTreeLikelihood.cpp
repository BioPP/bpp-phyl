//
// File: RHomogeneousMixedTreeLikelihood.cpp
// Authors:
//   David Fournier, Laurent Gueguen
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/


#include "RHomogeneousMixedTreeLikelihood.h"

// From the STL:
#include <iostream>

#include <cmath>
#include "../../PatternTools.h"

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/App/ApplicationTools.h>

using namespace bpp;
using namespace std;

RHomogeneousMixedTreeLikelihood::RHomogeneousMixedTreeLikelihood(
  const Tree& tree,
  shared_ptr<TransitionModelInterface> model,
  shared_ptr<DiscreteDistribution> rDist,
  bool checkRooted,
  bool verbose,
  bool usePatterns) :
  RHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose, usePatterns),
  treeLikelihoodsContainer_(),
  probas_()
{
  shared_ptr<MixedTransitionModelInterface> mixedmodel;
  if ((mixedmodel = dynamic_pointer_cast<MixedTransitionModelInterface>(model_)) == nullptr)
    throw Exception("Bad model: RHomogeneousMixedTreeLikelihood needs a MixedTransitionModel.");
  size_t s = mixedmodel->getNumberOfModels();
  for (size_t i = 0; i < s; ++i)
  {
    treeLikelihoodsContainer_.push_back(
      make_shared<RHomogeneousTreeLikelihood>(
	  tree,
	  unique_ptr<TransitionModelInterface>(mixedmodel->nModel(i).clone()),
	  rDist,
	  checkRooted,
	  false,
	  usePatterns));
    probas_.push_back(mixedmodel->getNProbability(i));
  }
}

RHomogeneousMixedTreeLikelihood::RHomogeneousMixedTreeLikelihood(
  const Tree& tree,
  const AlignmentDataInterface& data,
  std::shared_ptr<TransitionModelInterface> model,
  std::shared_ptr<DiscreteDistribution> rDist,
  bool checkRooted,
  bool verbose,
  bool usePatterns) :
  RHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose, usePatterns),
  treeLikelihoodsContainer_(),
  probas_()
{
  shared_ptr<MixedTransitionModelInterface> mixedmodel;

  if ((mixedmodel = dynamic_pointer_cast<MixedTransitionModelInterface>(model_)) == nullptr)
    throw Exception("Bad model: RHomogeneousMixedTreeLikelihood needs a MixedTransitionModel.");

  size_t s = mixedmodel->getNumberOfModels();
  for (size_t i = 0; i < s; i++)
  {
    treeLikelihoodsContainer_.push_back(
      make_shared<RHomogeneousTreeLikelihood>(
	  tree, 
	  unique_ptr<TransitionModelInterface>(mixedmodel->nModel(i).clone()),
	  rDist,
	  checkRooted,
	  false,
	  usePatterns));
    probas_.push_back(mixedmodel->getNProbability(i));
  }
  setData(data);
}

RHomogeneousMixedTreeLikelihood& RHomogeneousMixedTreeLikelihood::operator=(const RHomogeneousMixedTreeLikelihood& lik)
{
  RHomogeneousTreeLikelihood::operator=(lik);

  treeLikelihoodsContainer_.clear();
  probas_.clear();

  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_.push_back(shared_ptr<RHomogeneousTreeLikelihood>(lik.treeLikelihoodsContainer_[i]->clone()));
    probas_.push_back(lik.probas_[i]);
  }

  return *this;
}


RHomogeneousMixedTreeLikelihood::RHomogeneousMixedTreeLikelihood(const RHomogeneousMixedTreeLikelihood& lik) :
  RHomogeneousTreeLikelihood(lik),
  treeLikelihoodsContainer_(lik.treeLikelihoodsContainer_.size()),
  probas_(lik.probas_.size())
{
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); ++i)
  {
    treeLikelihoodsContainer_[i] = shared_ptr<RHomogeneousTreeLikelihood>(lik.treeLikelihoodsContainer_[i]->clone());
    probas_.push_back(lik.probas_[i]);
  }
}

RHomogeneousMixedTreeLikelihood::~RHomogeneousMixedTreeLikelihood() {}


void RHomogeneousMixedTreeLikelihood::initialize()
{
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); ++i)
  {
    treeLikelihoodsContainer_[i]->initialize();
  }

  RHomogeneousTreeLikelihood::initialize();
}

void RHomogeneousMixedTreeLikelihood::setData(const AlignmentDataInterface& sites)
{
  RHomogeneousTreeLikelihood::setData(sites);

  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); ++i)
  {
    treeLikelihoodsContainer_[i]->setData(sites);
  }
}


void RHomogeneousMixedTreeLikelihood::fireParameterChanged(const ParameterList& params)
{
  // checks in the model will change
  bool modelC = model_->getParameters().testParametersValues(params);

  applyParameters();
  auto mixedmodel = dynamic_pointer_cast<MixedTransitionModelInterface>(model_);
  size_t s = mixedmodel->getNumberOfModels();

  for (size_t i = 0; i < s; ++i)
  {
    ParameterList pl;
    pl.addParameters(mixedmodel->nModel(i).getParameters());
    pl.includeParameters(getParameters());

    if (modelC)
    {
      treeLikelihoodsContainer_[i]->setParameters(pl);
    }
    else
      treeLikelihoodsContainer_[i]->matchParametersValues(pl);
  }

  probas_ = mixedmodel->getProbabilities();
  minusLogLik_ = -getLogLikelihood();
}

void RHomogeneousMixedTreeLikelihood::computeTreeLikelihood()
{
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTreeLikelihood();
  }
}

/******************************************************************************
*                                   Likelihoods                              *
******************************************************************************/
double RHomogeneousMixedTreeLikelihood::getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  double res = 0;

  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    res += treeLikelihoodsContainer_[i]->getLikelihoodForASiteForARateClass(site, rateClass) * probas_[i];
  }

  return res;
}

double RHomogeneousMixedTreeLikelihood::getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  double x = getLikelihoodForASiteForARateClass(site, rateClass);
  if (x < 0)
    x = 0;
  return log(x);
}

double RHomogeneousMixedTreeLikelihood::getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  double res = 0;

  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    res += treeLikelihoodsContainer_[i]->getLikelihoodForASiteForARateClassForAState(site, rateClass, state) * probas_[i];
  }

  return res;
}

double RHomogeneousMixedTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  double x = getLikelihoodForASiteForARateClassForAState(site, rateClass, state);
  if (x < 0)
    x = 0;
  return log(x);
}


/******************************************************************************
*                           First Order Derivatives                          *
******************************************************************************/
double RHomogeneousMixedTreeLikelihood::getDLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  double res = 0;

  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    res += treeLikelihoodsContainer_[i]->getDLikelihoodForASiteForARateClass(site, rateClass) * probas_[i];
  }

  return res;
}

void RHomogeneousMixedTreeLikelihood::computeTreeDLikelihood(const string& variable)
{
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTreeDLikelihood(variable);
  }
}

/******************************************************************************
*                           Second Order Derivatives                          *
******************************************************************************/
double RHomogeneousMixedTreeLikelihood::getD2LikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  double res = 0;

  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    res += treeLikelihoodsContainer_[i]->getD2LikelihoodForASiteForARateClass(site, rateClass) * probas_[i];
  }

  return res;
}


void RHomogeneousMixedTreeLikelihood::computeTreeD2Likelihood(const string& variable)
{
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTreeD2Likelihood(variable);
  }
}

void RHomogeneousMixedTreeLikelihood::computeSubtreeLikelihood(const Node* node)
{
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeSubtreeLikelihood(node);
  }
}

void RHomogeneousMixedTreeLikelihood::computeDownSubtreeDLikelihood(const Node* node)
{
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeDownSubtreeDLikelihood(node);
  }
}

void RHomogeneousMixedTreeLikelihood::computeDownSubtreeD2Likelihood(const Node* node)
{
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeDownSubtreeD2Likelihood(node);
  }
}


void RHomogeneousMixedTreeLikelihood::computeAllTransitionProbabilities()
{
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeAllTransitionProbabilities();
  }
}


void RHomogeneousMixedTreeLikelihood::computeTransitionProbabilitiesForNode(const Node* node)
{
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTransitionProbabilitiesForNode(node);
  }
}


void RHomogeneousMixedTreeLikelihood::displayLikelihood(const Node* node)
{
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->displayLikelihood(node);
  }
}
