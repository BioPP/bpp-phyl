//
// File: DRHomogeneousMixedTreeLikelihood.cpp
// Created by: Laurent Gueguen
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

#include "DRHomogeneousMixedTreeLikelihood.h"


// From the STL:
#include <iostream>

#include <math.h>
#include "../PatternTools.h"

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/App/ApplicationTools.h>

using namespace bpp;
using namespace std;

DRHomogeneousMixedTreeLikelihood::DRHomogeneousMixedTreeLikelihood(
  const Tree& tree,
  SubstitutionModel* model,
  DiscreteDistribution* rDist,
  bool checkRooted,
  bool verbose)  throw (Exception) :
  DRHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose), treeLikelihoodsContainer_(),
  probas_()
{
  MixedSubstitutionModel* mixedmodel;

  if ((mixedmodel = dynamic_cast<MixedSubstitutionModel*>(model_)) == NULL)
    throw Exception("Bad model: DRHomogeneousMixedTreeLikelihood needs a MixedSubstitutionModel.");

  unsigned int s = mixedmodel->getNumberOfModels();
  for (unsigned int i = 0; i < s; i++)
  {
    treeLikelihoodsContainer_.push_back(
      new DRHomogeneousTreeLikelihood(tree, mixedmodel->getNModel(i), rDist, checkRooted, false));
    probas_.push_back(mixedmodel->getNProbability(i));
 }
}

DRHomogeneousMixedTreeLikelihood::DRHomogeneousMixedTreeLikelihood(
  const Tree& tree,
  const SiteContainer& data,
  SubstitutionModel* model,
  DiscreteDistribution* rDist,
  bool checkRooted,
  bool verbose)
throw (Exception) :
  DRHomogeneousTreeLikelihood(tree, model,rDist, checkRooted, verbose), treeLikelihoodsContainer_(),
  probas_()
{
  MixedSubstitutionModel* mixedmodel;

  if ((mixedmodel = dynamic_cast<MixedSubstitutionModel*>(model_)) == NULL)
    throw Exception("Bad model: DRHomogeneousMixedTreeLikelihood needs a MixedSubstitutionModel.");

  int s = mixedmodel->getNumberOfModels();

  for (int i = 0; i < s; i++)
  {
   treeLikelihoodsContainer_.push_back(
     new DRHomogeneousTreeLikelihood(tree, mixedmodel->getNModel(i), rDist, checkRooted, false));
   probas_.push_back(mixedmodel->getNProbability(i));
  }
  setData(data);
}


DRHomogeneousMixedTreeLikelihood& DRHomogeneousMixedTreeLikelihood::operator=(const DRHomogeneousMixedTreeLikelihood& lik)
{
  DRHomogeneousTreeLikelihood::operator=(lik);

  treeLikelihoodsContainer_.clear();
  probas_.clear();

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_.push_back(lik.treeLikelihoodsContainer_[i]->clone());
    probas_.push_back(lik.probas_[i]);
  }

  return *this;
}

DRHomogeneousMixedTreeLikelihood::DRHomogeneousMixedTreeLikelihood(const DRHomogeneousMixedTreeLikelihood& lik) :
  DRHomogeneousTreeLikelihood(lik), treeLikelihoodsContainer_(lik.treeLikelihoodsContainer_.size()), probas_(lik.probas_.size())
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_.push_back(lik.treeLikelihoodsContainer_[i]->clone());
    probas_.push_back(lik.probas_[i]);
  }
}

DRHomogeneousMixedTreeLikelihood::~DRHomogeneousMixedTreeLikelihood()
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    delete treeLikelihoodsContainer_[i];
  }
}


void DRHomogeneousMixedTreeLikelihood::initialize() throw (Exception)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->initialize();
  }
  DRHomogeneousTreeLikelihood::initialize();
}

void DRHomogeneousMixedTreeLikelihood::setData(const SiteContainer& sites) throw (Exception)
{
  DRHomogeneousTreeLikelihood::setData(sites);

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->setData(sites);
  }
}


void DRHomogeneousMixedTreeLikelihood::fireParameterChanged(const ParameterList& params)
{
  applyParameters();

  MixedSubstitutionModel* mixedmodel = dynamic_cast<MixedSubstitutionModel*>(model_);

  unsigned int s = mixedmodel->getNumberOfModels();

  const SubstitutionModel* pm;
  for (unsigned int i = 0; i < s; i++)
  {
    ParameterList pl;
    pm = mixedmodel->getNModel(i);
    pl.addParameters(pm->getParameters());
    pl.includeParameters(getParameters());
    treeLikelihoodsContainer_[i]->matchParametersValues(pl);
  }
  probas_=mixedmodel->getProbabilities();

  minusLogLik_ = -getLogLikelihood();
}

void DRHomogeneousMixedTreeLikelihood::resetLikelihoodArrays(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->resetLikelihoodArrays(node);
  }
}

void DRHomogeneousMixedTreeLikelihood::computeTreeLikelihood()
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTreeLikelihood();
  }
}

/******************************************************************************
 *                           Likelihoods                          *
 ******************************************************************************/  

double DRHomogeneousMixedTreeLikelihood::getLikelihood() const
{
  double res=0;
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++) {
    res+=treeLikelihoodsContainer_[i]->getLikelihood() * probas_[i];
  }

  return res;
}

double DRHomogeneousMixedTreeLikelihood::getLogLikelihood() const
{
  double x=getLikelihood();
  if (x<0) x=0;
  return log(x);
}


double DRHomogeneousMixedTreeLikelihood::getLikelihoodForASite(unsigned int site) const
{
  double res=0;
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++) {
    res+=treeLikelihoodsContainer_[i]->getLikelihoodForASite(site) * probas_[i];
  }

  return res;
}

double DRHomogeneousMixedTreeLikelihood::getLogLikelihoodForASite(unsigned int site) const
{
  double x=getLikelihoodForASite(site);
  if (x<0) x=0;
  return log(x);
}

double DRHomogeneousMixedTreeLikelihood::getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
  double res=0;
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++) {
    res+=treeLikelihoodsContainer_[i]->getLikelihoodForASiteForARateClass(site, rateClass) * probas_[i];
  }

  return res;
}

double DRHomogeneousMixedTreeLikelihood::getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
  double x=getLikelihoodForASiteForARateClass(site, rateClass);
  if (x<0) x=0;
  return log(x);
}

double DRHomogeneousMixedTreeLikelihood::getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
  double res=0;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++) {
    res+=treeLikelihoodsContainer_[i]->getLikelihoodForASiteForARateClassForAState(site, rateClass, state) * probas_[i];
  }

  return res;
}

double DRHomogeneousMixedTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
  double x=getLikelihoodForASiteForARateClassForAState(site, rateClass, state);
  if (x<0) x=0;
  return log(x);
}



void DRHomogeneousMixedTreeLikelihood::computeSubtreeLikelihoodPostfix(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeSubtreeLikelihoodPostfix(node);
  }
}

void DRHomogeneousMixedTreeLikelihood::computeSubtreeLikelihoodPrefix(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeSubtreeLikelihoodPostfix(node);
  }
}

void DRHomogeneousMixedTreeLikelihood::computeRootLikelihood()
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeRootLikelihood();
  }
}

void DRHomogeneousMixedTreeLikelihood::computeLikelihoodAtNode_(const Node* node, VVVdouble& likelihoodArray) const
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeLikelihoodAtNode_(node, likelihoodArray);
  }
}

void DRHomogeneousMixedTreeLikelihood::computeLikelihoodFromArrays(
  const vector<const VVVdouble*>& iLik,
  const vector<const VVVdouble*>& tProb,
  VVVdouble& oLik,
  unsigned int nbNodes,
  unsigned int nbDistinctSites,
  unsigned int nbClasses,
  unsigned int nbStates,
  bool reset)
{
  cerr << "Invalid use of static method DRHomogeneousMixedTreeLikelihood::computeLikelihoodFromArrays" << endl;
  exit(0);
}

void DRHomogeneousMixedTreeLikelihood::computeLikelihoodFromArrays(
  const std::vector<const VVVdouble*>& iLik,
  const std::vector<const VVVdouble*>& tProb,
  const VVVdouble* iLikR,
  const VVVdouble* tProbR,
  VVVdouble& oLik,
  unsigned int nbNodes,
  unsigned int nbDistinctSites,
  unsigned int nbClasses,
  unsigned int nbStates,
  bool reset)
{
  cerr << "Invalid use of static method DRHomogeneousMixedTreeLikelihood::computeLikelihoodFromArrays" << endl;
  exit(0);
}

/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/  

void DRHomogeneousMixedTreeLikelihood::computeTreeDLikelihoods()
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
    {
      treeLikelihoodsContainer_[i]->computeTreeDLikelihoods();
    }
}

double DRHomogeneousMixedTreeLikelihood::getFirstOrderDerivative(const string& variable) const throw (Exception)
{
  double res=0;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++) {
    res+=treeLikelihoodsContainer_[i]->getFirstOrderDerivative(variable) * probas_[i];
  }

    return res;
}

/******************************************************************************
 *                           Second Order Derivatives                          *
 ******************************************************************************/  

void DRHomogeneousMixedTreeLikelihood::computeTreeD2LikelihoodAtNode(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
    {
      treeLikelihoodsContainer_[i]->computeTreeD2LikelihoodAtNode(node);
    }
}

void DRHomogeneousMixedTreeLikelihood::computeTreeD2Likelihoods()
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
    {
      treeLikelihoodsContainer_[i]->computeTreeD2Likelihoods();
    }
}

double DRHomogeneousMixedTreeLikelihood::getSecondOrderDerivative(const string& variable) const throw (Exception)
{
  double res=0;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++) {
    res+=treeLikelihoodsContainer_[i]->getSecondOrderDerivative(variable) * probas_[i];
  }
  
  return res;
}

void DRHomogeneousMixedTreeLikelihood::displayLikelihood(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->displayLikelihood(node);
  }
}

