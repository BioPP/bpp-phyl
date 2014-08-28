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
  bool verbose,
  bool rootArray)  throw (Exception) :
  DRHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  treeLikelihoodsContainer_(),
  probas_(),
  rootArray_(rootArray)
{
  MixedSubstitutionModel* mixedmodel;

  if ((mixedmodel = dynamic_cast<MixedSubstitutionModel*>(model_)) == NULL)
    throw Exception("Bad model: DRHomogeneousMixedTreeLikelihood needs a MixedSubstitutionModel.");

  size_t s = mixedmodel->getNumberOfModels();
  for (size_t i = 0; i < s; i++)
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
  bool verbose,
  bool rootArray)
throw (Exception) :
  DRHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  treeLikelihoodsContainer_(),
  probas_(),
  rootArray_(rootArray)
{
  MixedSubstitutionModel* mixedmodel;

  if ((mixedmodel = dynamic_cast<MixedSubstitutionModel*>(model_)) == NULL)
    throw Exception("Bad model: DRHomogeneousMixedTreeLikelihood needs a MixedSubstitutionModel.");

  size_t s = mixedmodel->getNumberOfModels();

  for (size_t i = 0; i < s; i++)
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

  rootArray_=lik.rootArray_;

  return *this;
}

DRHomogeneousMixedTreeLikelihood::DRHomogeneousMixedTreeLikelihood(const DRHomogeneousMixedTreeLikelihood& lik) :
  DRHomogeneousTreeLikelihood(lik),
  treeLikelihoodsContainer_(lik.treeLikelihoodsContainer_.size()),
  probas_(lik.probas_.size()),
  rootArray_(lik.rootArray_)
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
  if(rootArray_)
    computeRootLikelihood();
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

  size_t s = mixedmodel->getNumberOfModels();

  const SubstitutionModel* pm;
  for (size_t i = 0; i < s; i++)
  {
    ParameterList pl;
    pm = mixedmodel->getNModel(i);
    pl.addParameters(pm->getParameters());
    pl.includeParameters(getParameters());
    treeLikelihoodsContainer_[i]->matchParametersValues(pl);
  }
  probas_ = mixedmodel->getProbabilities();

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
  if(rootArray_)
    computeRootLikelihood();
}

/******************************************************************************
*                           Likelihoods                          *
******************************************************************************/

double DRHomogeneousMixedTreeLikelihood::getLikelihood() const
{
  double l = 1.;
  vector<Vdouble*> llik;
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    llik.push_back(&treeLikelihoodsContainer_[i]->likelihoodData_->getRootRateSiteLikelihoodArray());
  }

  double x;
  const vector<unsigned int> * w = &likelihoodData_->getWeights();
  for (unsigned int i = 0; i < nbDistinctSites_; i++)
  {
    x = 0;
    for (unsigned int j = 0; j < treeLikelihoodsContainer_.size(); j++)
    {
      x += (*llik[j])[i] * probas_[j];
    }
    l *= std::pow(x, (int)(*w)[i]);
  }
  return l;
}

double DRHomogeneousMixedTreeLikelihood::getLogLikelihood() const
{
  double ll = 0;

  vector<Vdouble*> llik;
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    llik.push_back(&treeLikelihoodsContainer_[i]->likelihoodData_->getRootRateSiteLikelihoodArray());
  }

  double x;
  const vector<unsigned int> * w = &likelihoodData_->getWeights();
  vector<double> la(nbDistinctSites_);
  for (unsigned int i = 0; i < nbDistinctSites_; i++)
  {
    x = 0;
    for (unsigned int j = 0; j < treeLikelihoodsContainer_.size(); j++)
    {
      x += (*llik[j])[i] * probas_[j];
    }
    la[i] = (*w)[i] * log(x);
  }
  sort(la.begin(), la.end());
  for (size_t i = nbDistinctSites_; i > 0; i--)
  {
    ll += la[i - 1];
  }

  return ll;
}


double DRHomogeneousMixedTreeLikelihood::getLikelihoodForASite(size_t site) const
{
  double res = 0;
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    res += treeLikelihoodsContainer_[i]->getLikelihoodForASite(site) * probas_[i];
  }

  return res;
}

double DRHomogeneousMixedTreeLikelihood::getLogLikelihoodForASite(size_t site) const
{
  double x = getLikelihoodForASite(site);
  if (x < 0) x = 0;
  return log(x);
}

double DRHomogeneousMixedTreeLikelihood::getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  double res = 0;
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    res += treeLikelihoodsContainer_[i]->getLikelihoodForASiteForARateClass(site, rateClass) * probas_[i];
  }

  return res;
}

double DRHomogeneousMixedTreeLikelihood::getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  double x = getLikelihoodForASiteForARateClass(site, rateClass);
  if (x < 0) x = 0;
  return log(x);
}

double DRHomogeneousMixedTreeLikelihood::getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  double res = 0;

  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    res += treeLikelihoodsContainer_[i]->getLikelihoodForASiteForARateClassForAState(site, rateClass, state) * probas_[i];
  }

  return res;
}

double DRHomogeneousMixedTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  double x = getLikelihoodForASiteForARateClassForAState(site, rateClass, state);
  if (x < 0) x = 0;
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
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeRootLikelihood();
  }
}

void DRHomogeneousMixedTreeLikelihood::computeLikelihoodAtNode_(const Node* node, VVVdouble& likelihoodArray, const Node* sonNode) const
{
  likelihoodArray.resize(nbDistinctSites_);
  for (size_t i = 0; i < nbDistinctSites_; i++){
    VVdouble* likelihoodArray_i = &likelihoodArray[i];
    likelihoodArray_i->resize(nbClasses_);
    for (size_t c = 0; c < nbClasses_; c++) {
      Vdouble* likelihoodArray_i_c = &(*likelihoodArray_i)[c];
      likelihoodArray_i_c->resize(nbStates_);
      for (size_t x = 0; x < nbStates_; x++)
        (*likelihoodArray_i_c)[x] = 0;
    }
  }

  VVVdouble lArray;
  for (size_t nm = 0; nm < treeLikelihoodsContainer_.size(); nm++)
  {
    treeLikelihoodsContainer_[nm]->computeLikelihoodAtNode_(node, lArray, sonNode);
    
    for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        VVdouble* likelihoodArray_i = &likelihoodArray[i];
        VVdouble* lArray_i = &lArray[i];
        
        for (size_t c = 0; c < nbClasses_; c++)
          {
            Vdouble* likelihoodArray_i_c = &(*likelihoodArray_i)[c];
            Vdouble* lArray_i_c = &(*lArray_i)[c];
            for (size_t x = 0; x < nbStates_; x++)
              (*likelihoodArray_i_c)[x] += (*lArray_i_c)[x] * probas_[nm];
         }
      }
    
  }
}

/******************************************************************************
*                           First Order Derivatives                          *
******************************************************************************/

void DRHomogeneousMixedTreeLikelihood::computeTreeDLikelihoods()
{
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->computeTreeDLikelihoods();
  }
}

double DRHomogeneousMixedTreeLikelihood::getFirstOrderDerivative(const std::string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("DRHomogeneousTreeLikelihood::getFirstOrderDerivative().", variable);
  if (getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameters are not implemented.");
  }
  if (getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }

  //
  // Computation for branch lengths:
  //

  // Get the node with the branch whose length must be derivated:
  unsigned int brI = TextTools::to<unsigned int>(variable.substr(5));
  const Node* branch = nodes_[brI];
  vector< Vdouble*> _vdLikelihoods_branch;
  for (size_t i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    _vdLikelihoods_branch.push_back(&treeLikelihoodsContainer_[i]->likelihoodData_->getDLikelihoodArray(branch->getId()));
  }

  double d = 0;
  double x;
  const vector<unsigned int> * w = &likelihoodData_->getWeights();
  for (unsigned int i = 0; i < nbDistinctSites_; i++)
  {
    x = 0;
    for (unsigned int j = 0; j < treeLikelihoodsContainer_.size(); j++)
    {
      x += (*_vdLikelihoods_branch[j])[i] * probas_[j];
    }
    d += (*w)[i] * x;
  }

  return -d;
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

double DRHomogeneousMixedTreeLikelihood::getSecondOrderDerivative(const std::string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("DRHomogeneousTreeLikelihood::getFirstOrderDerivative().", variable);
  if (getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameters are not implemented.");
  }
  if (getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }

  //
  // Computation for branch lengths:
  //

  // Get the node with the branch whose length must be derivated:
  unsigned int brI = TextTools::to<unsigned int>(variable.substr(5));
  const Node* branch = nodes_[brI];
  vector< Vdouble*> _vdLikelihoods_branch, _vd2Likelihoods_branch;
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    _vdLikelihoods_branch.push_back(&treeLikelihoodsContainer_[i]->likelihoodData_->getDLikelihoodArray(branch->getId()));
    _vd2Likelihoods_branch.push_back(&treeLikelihoodsContainer_[i]->likelihoodData_->getD2LikelihoodArray(branch->getId()));
  }

  double d = 0;
  double x, x2;
  const vector<unsigned int> * w = &likelihoodData_->getWeights();
  for (unsigned int i = 0; i < nbDistinctSites_; i++)
  {
    x = 0;
    x2 = 0;
    for (unsigned int j = 0; j < treeLikelihoodsContainer_.size(); j++)
    {
      x += (*_vdLikelihoods_branch[j])[i] * probas_[j];
    }
    for (unsigned int j = 0; j < treeLikelihoodsContainer_.size(); j++)
    {
      x2 += (*_vd2Likelihoods_branch[j])[i] * probas_[j];
    }

    d += (*w)[i] * (x2 - pow(x, 2));
  }

  return -d;
}


void DRHomogeneousMixedTreeLikelihood::displayLikelihood(const Node* node)
{
  for (unsigned int i = 0; i < treeLikelihoodsContainer_.size(); i++)
  {
    treeLikelihoodsContainer_[i]->displayLikelihood(node);
  }
}

