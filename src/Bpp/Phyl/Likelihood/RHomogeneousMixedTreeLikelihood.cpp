//
// File: RHomogeneousMixedTreeLikelihood.h
// Created by: David Fournier, Laurent Gueguen
//

/*
   Copyright or © or Copr. CNRS, (November 16, 2004)

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

#include "RHomogeneousMixedTreeLikelihood.h"


// From the STL:
#include <iostream>

#include <math.h>
#include "PatternTools.h"

#include <NumCalc/VectorTools.h>

#include <Utils/ApplicationTools.h>

using namespace bpp;
using namespace std;

RHomogeneousMixedTreeLikelihood::RHomogeneousMixedTreeLikelihood(
  const Tree& tree,
  SubstitutionModel* model,
  DiscreteDistribution* rDist,
  bool checkRooted,
  bool verbose,
  bool usePatterns) throw (Exception) :
  RHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose, usePatterns),
  treeLikelihoodsContainer_(),
  probas_()
{
  MixedSubstitutionModel* mixedmodel;
  if ((mixedmodel = dynamic_cast<MixedSubstitutionModel*>(model_)) == 0)
    throw Exception("Bad model: RHomogeneousMixedTreeLikelihood needs a MixedSubstitutionModel.");
  unsigned int s = mixedmodel->getNumberOfModels();
  for (unsigned int i = 0; i < s; i++)
  {
   treeLikelihoodsContainer_.push_back(
       new RHomogeneousTreeLikelihood(tree, mixedmodel->getNModel(i), rDist, checkRooted, false, usePatterns));
   probas_.push_back(1.0/s);
  }
}

RHomogeneousMixedTreeLikelihood::RHomogeneousMixedTreeLikelihood(
  const Tree& tree,
  const SiteContainer& data,
  SubstitutionModel* model,
  DiscreteDistribution* rDist,
  bool checkRooted,
  bool verbose,
  bool usePatterns) throw (Exception) :
  RHomogeneousTreeLikelihood(tree, model,rDist, checkRooted, verbose, usePatterns),
  treeLikelihoodsContainer_(),
  probas_()
{
  MixedSubstitutionModel* mixedmodel;

  if ((mixedmodel = dynamic_cast<MixedSubstitutionModel*>(model_)) == 0)
    throw Exception("Bad model: RHomogeneousMixedTreeLikelihood needs a MixedSubstitutionModel.");

  unsigned int s = mixedmodel->getNumberOfModels();

  for (unsigned int i = 0; i < s; i++)
  {
    treeLikelihoodsContainer_.push_back(
        new RHomogeneousTreeLikelihood(tree, mixedmodel->getNModel(i), rDist, checkRooted, false, usePatterns));
    probas_.push_back(1.0/s);
  }
  setData(data);
}

RHomogeneousMixedTreeLikelihood& RHomogeneousMixedTreeLikelihood::operator=(const RHomogeneousMixedTreeLikelihood& lik)
{
  RHomogeneousTreeLikelihood::operator=(lik);

  treeLikelihoodsContainer_.clear();
  probas_.clear();
  
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   treeLikelihoodsContainer_.push_back(lik.treeLikelihoodsContainer_[i]->clone());
   probas_.push_back(lik.probas_[i]);
  }

  return *this;
}


RHomogeneousMixedTreeLikelihood::RHomogeneousMixedTreeLikelihood(const RHomogeneousMixedTreeLikelihood& lik) :
  RHomogeneousTreeLikelihood(lik),
  treeLikelihoodsContainer_(lik.treeLikelihoodsContainer_.size()),
  probas_(lik.probas_.size())
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   treeLikelihoodsContainer_[i] = lik.treeLikelihoodsContainer_[i]->clone();
   probas_.push_back(lik.probas_[i]);
  }
}

RHomogeneousMixedTreeLikelihood::~RHomogeneousMixedTreeLikelihood()
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    delete treeLikelihoodsContainer_[i];
  }
}


void RHomogeneousMixedTreeLikelihood::initialize() throw (Exception)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->initialize();
  }
  RHomogeneousTreeLikelihood::initialize();
}

void RHomogeneousMixedTreeLikelihood::setData(const SiteContainer& sites) throw (Exception)
{
   RHomogeneousTreeLikelihood::setData(sites);

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->setData(sites);
  }
}


double RHomogeneousMixedTreeLikelihood::getLogLikelihood() const
{
   vector<double> reslog;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   reslog.push_back(treeLikelihoodsContainer_[i]->getLogLikelihood());
  }

  return VectorTools::logsumexp(reslog,probas_);
}

double RHomogeneousMixedTreeLikelihood::getLikelihood() const
{
  return exp(getLogLikelihood());
}

double RHomogeneousMixedTreeLikelihood::getLikelihoodForASite(unsigned int site) const
{
  return exp(getLogLikelihoodForASite(site));
}

double RHomogeneousMixedTreeLikelihood::getLogLikelihoodForASite(unsigned int site) const
{
   vector<double> reslog;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   reslog.push_back(treeLikelihoodsContainer_[i]->getLogLikelihoodForASite(site));
  }

  return VectorTools::logsumexp(reslog,probas_);
}

double RHomogeneousMixedTreeLikelihood::getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
  return exp(getLogLikelihoodForASiteForARateClass(site, rateClass));
}

double RHomogeneousMixedTreeLikelihood::getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
   vector<double> reslog;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   reslog.push_back(treeLikelihoodsContainer_[i]->getLogLikelihoodForASiteForARateClass(site, rateClass));
  }

  return VectorTools::logsumexp(reslog,probas_);
}

double RHomogeneousMixedTreeLikelihood::getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
  return exp(getLogLikelihoodForASiteForARateClassForAState(site, rateClass, state));
}

double RHomogeneousMixedTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
   vector<double> reslog;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   reslog.push_back(treeLikelihoodsContainer_[i]->getLogLikelihoodForASiteForARateClassForAState(site, rateClass, state));
  }

  return VectorTools::logsumexp(reslog,probas_);
}

void RHomogeneousMixedTreeLikelihood::fireParameterChanged(const ParameterList& params)
{
  applyParameters();

  MixedSubstitutionModel* mixedmodel = dynamic_cast<MixedSubstitutionModel*>(model_);

  unsigned int s = mixedmodel->getNumberOfModels();

  const SubstitutionModel* pm;
  for (unsigned int i = 0; i < s; i++)
  {
    pm = mixedmodel->getNModel(i);
    treeLikelihoodsContainer_[i]->matchParametersValues(pm->getParameters());
    treeLikelihoodsContainer_[i]->matchParametersValues(getParameters());
  }
  probas_=mixedmodel->getProbabilities();

  minusLogLik_ = -getLogLikelihood();
}

double RHomogeneousMixedTreeLikelihood::getFirstOrderDerivative(const string& variable) const throw (Exception)
{
   vector<double> rescontainer;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   rescontainer.push_back(treeLikelihoodsContainer_[i]->getFirstOrderDerivative(variable));
  }

  return VectorTools::sum<double>(rescontainer,probas_);
}

double RHomogeneousMixedTreeLikelihood::getSecondOrderDerivative(const string& variable) const throw (Exception)
{
   vector<double> rescontainer;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   rescontainer.push_back(treeLikelihoodsContainer_[i]->getSecondOrderDerivative(variable));
  }

  return VectorTools::sum<double>(rescontainer,probas_);
}

void RHomogeneousMixedTreeLikelihood::computeTreeLikelihood()
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTreeLikelihood();
  }
}

double RHomogeneousMixedTreeLikelihood::getDLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
   vector<double> rescontainer;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   rescontainer.push_back(treeLikelihoodsContainer_[i]->getDLikelihoodForASiteForARateClass(site,rateClass));
  }

  return VectorTools::sum<double>(rescontainer,probas_);
}

double RHomogeneousMixedTreeLikelihood::getDLikelihoodForASite(unsigned int site) const
{
   vector<double> rescontainer;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   rescontainer.push_back(treeLikelihoodsContainer_[i]->getDLikelihoodForASite(site));
  }

  return VectorTools::sum<double>(rescontainer,probas_);
}

void RHomogeneousMixedTreeLikelihood::computeTreeDLikelihood(const string& variable)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTreeDLikelihood(variable);
  }
}


double RHomogeneousMixedTreeLikelihood::getD2LikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
   vector<double> rescontainer;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   rescontainer.push_back(treeLikelihoodsContainer_[i]->getD2LikelihoodForASiteForARateClass(site,rateClass));
  }

  return VectorTools::sum<double>(rescontainer,probas_);
}

double RHomogeneousMixedTreeLikelihood::getD2LikelihoodForASite(unsigned int site) const
{
   vector<double> rescontainer;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
   rescontainer.push_back(treeLikelihoodsContainer_[i]->getD2LikelihoodForASite(site));
  }

  return VectorTools::sum<double>(rescontainer,probas_);
}

void RHomogeneousMixedTreeLikelihood::computeTreeD2Likelihood(const string& variable)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTreeD2Likelihood(variable);
  }
}

void RHomogeneousMixedTreeLikelihood::computeSubtreeLikelihood(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeSubtreeLikelihood(node);
  }
}

void RHomogeneousMixedTreeLikelihood::computeDownSubtreeDLikelihood(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeDownSubtreeDLikelihood(node);
  }
}

void RHomogeneousMixedTreeLikelihood::computeDownSubtreeD2Likelihood(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeDownSubtreeD2Likelihood(node);
  }
}


void RHomogeneousMixedTreeLikelihood::computeAllTransitionProbabilities()
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeAllTransitionProbabilities();
  }
}


void RHomogeneousMixedTreeLikelihood::computeTransitionProbabilitiesForNode(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTransitionProbabilitiesForNode(node);
  }
}


void RHomogeneousMixedTreeLikelihood::displayLikelihood(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->displayLikelihood(node);
  }
}
