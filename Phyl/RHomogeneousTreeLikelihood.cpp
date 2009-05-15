//
// File: RHomogeneousTreeLikelihood.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 17 18:14:51 2003
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

#include "RHomogeneousTreeLikelihood.h"
#include "PatternTools.h"

// From Utils:
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

RHomogeneousTreeLikelihood::RHomogeneousTreeLikelihood(
  const Tree & tree,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose,
  bool usePatterns)
throw (Exception):
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  _likelihoodData(NULL)
{
  _init(usePatterns);
}

/******************************************************************************/

RHomogeneousTreeLikelihood::RHomogeneousTreeLikelihood(
  const Tree & tree,
  const SiteContainer & data,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose,
  bool usePatterns)
throw (Exception):
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  _likelihoodData(NULL)
{
  _init(usePatterns);
  setData(data);
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::_init(bool usePatterns) throw (Exception)
{
  _likelihoodData = new DRASRTreeLikelihoodData(*_tree, _rateDistribution->getNumberOfCategories(), usePatterns);
}

/******************************************************************************/

RHomogeneousTreeLikelihood::RHomogeneousTreeLikelihood(
    const RHomogeneousTreeLikelihood & lik):
  AbstractHomogeneousTreeLikelihood(lik),
  _likelihoodData(NULL)
{
  _likelihoodData = dynamic_cast<DRASRTreeLikelihoodData *>(lik._likelihoodData->clone());
  _likelihoodData->setTree(*_tree);
}

/******************************************************************************/

RHomogeneousTreeLikelihood & RHomogeneousTreeLikelihood::operator=(
    const RHomogeneousTreeLikelihood & lik)
{
  AbstractHomogeneousTreeLikelihood::operator=(lik);
  if(_likelihoodData) delete _likelihoodData;
  _likelihoodData = dynamic_cast<DRASRTreeLikelihoodData *>(lik._likelihoodData->clone());
  _likelihoodData->setTree(*_tree);
  return *this;
}

/******************************************************************************/

RHomogeneousTreeLikelihood::~RHomogeneousTreeLikelihood()
{ 
  delete _likelihoodData;
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::setData(const SiteContainer & sites) throw (Exception)
{
  if(_data) delete _data;
  _data = PatternTools::getSequenceSubset(sites, *_tree->getRootNode());
  if(verbose_) ApplicationTools::displayTask("Initializing data structure");
  _likelihoodData->initLikelihoods(*_data, *model_);
  if(verbose_) ApplicationTools::displayTaskDone();

  nbSites_ = _likelihoodData->getNumberOfSites();
  nbDistinctSites_ = _likelihoodData->getNumberOfDistinctSites();
  nbStates_ = _likelihoodData->getNumberOfStates();
  
  if(verbose_) ApplicationTools::displayResult("Number of distinct sites",
      TextTools::toString(nbDistinctSites_));
  _initialized = false;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLikelihood() const
{
  double l = 1.;
  for(unsigned int i = 0; i < nbSites_; i++)
  {
    l *= getLikelihoodForASite(i);
  }
  return l;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLogLikelihood() const
{
  double ll = 0;
  vector<double> la(nbSites_);
  for(unsigned int i = 0; i < nbSites_; i++)
  {
    la[i] = getLogLikelihoodForASite(i);
  }
  sort(la.begin(), la.end());
  for(unsigned int i = nbSites_; i > 0; i--)
    ll += la[i-1];
  return ll;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLikelihoodForASite(unsigned int site) const
{
  double l = 0;
  for(unsigned int i = 0; i < nbClasses_; i++)
  {
    l += getLikelihoodForASiteForARateClass(site, i) * _rateDistribution->getProbability(i);
  }
  return l;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLogLikelihoodForASite(unsigned int site) const
{
  double l = 0;
  for(unsigned int i = 0; i < nbClasses_; i++)
  {
    l += getLikelihoodForASiteForARateClass(site, i) * _rateDistribution->getProbability(i);
  }
  //if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
  if(l<0) l=0; //May happen because of numerical errors.
  return log(l);
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
  double l = 0;
  Vdouble * la = & _likelihoodData->getLikelihoodArray(_tree->getRootNode()->getId())[_likelihoodData->getRootArrayPosition(site)][rateClass];
  for(unsigned int i = 0; i < nbStates_; i++)
  {
    l += (* la)[i] * rootFreqs_[i];
  }
  return l;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
  double l = 0;
  Vdouble * la = & _likelihoodData->getLikelihoodArray(_tree->getRootNode()->getId())[_likelihoodData->getRootArrayPosition(site)][rateClass];
  for(unsigned int i = 0; i < nbStates_; i++) {
    l += (* la)[i] * rootFreqs_[i];
  }
  //if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
  return log(l);
}

/******************************************************************************/  

double RHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
  return _likelihoodData->getLikelihoodArray(_tree->getRootNode()->getId())[_likelihoodData->getRootArrayPosition(site)][rateClass][state];
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
  return log(_likelihoodData->getLikelihoodArray(_tree->getRootNode()->getId())[_likelihoodData->getRootArrayPosition(site)][rateClass][state]);
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::setParameters(const ParameterList & parameters)
  throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::fireParameterChanged(const ParameterList & params)
{
  applyParameters();

  if(_rateDistribution->getParameters().getCommonParametersWith(params).size() > 0
  || model_->getParameters().getCommonParametersWith(params).size() > 0)
  {
    //Rate parameter changed, need to recompute all probs:
    computeAllTransitionProbabilities();
  }
  else if(params.size() > 0)
  {
    //We may save some computations:
    for(unsigned int i = 0; i < params.size(); i++)
    {
      string s = params[i]->getName();
      if(s.substr(0,5) == "BrLen")
      {
        //Branch length parameter:
        computeTransitionProbabilitiesForNode(nodes_[TextTools::to<unsigned int>(s.substr(5))]);
      }
    }
    rootFreqs_ = model_->getFrequencies();
  }

  computeTreeLikelihood();
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getValue() const
throw (Exception)
{
  if(!isInitialized()) throw Exception("RHomogeneousTreeLikelihood::getValue(). Instance is not initialized.");
  return - getLogLikelihood();
}

/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/  

double RHomogeneousTreeLikelihood::getDLikelihoodForASiteForARateClass(
  unsigned int site,
  unsigned int rateClass) const
{
  double dl = 0;
  Vdouble * dla = & _likelihoodData->getDLikelihoodArray(_tree->getRootNode()->getId())[_likelihoodData->getRootArrayPosition(site)][rateClass];
  for(unsigned int i = 0; i < nbStates_; i++)
  {
    dl += (* dla)[i] * rootFreqs_[i];
  }
  return dl;
}

/******************************************************************************/  

double RHomogeneousTreeLikelihood::getDLikelihoodForASite(unsigned int site) const
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for(unsigned int i = 0; i < nbClasses_; i++)
    dl += getDLikelihoodForASiteForARateClass(site, i) * _rateDistribution->getProbability(i);
  return dl;
}

/******************************************************************************/  

double RHomogeneousTreeLikelihood::getDLogLikelihoodForASite(unsigned int site) const
{
  // d(f(g(x)))/dx = dg(x)/dx . df(g(x))/dg :
  return getDLikelihoodForASite(site) / getLikelihoodForASite(site);
}

/******************************************************************************/  

double RHomogeneousTreeLikelihood::getDLogLikelihood() const
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for(unsigned int i = 0; i < nbSites_; i++)
    dl += getDLogLikelihoodForASite(i);
  return dl;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getFirstOrderDerivative(const string & variable) const
throw (Exception)
{ 
  const Parameter * p = &getParameter(variable);
  if(getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameter are not implemented.");
  }
  if(getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }
  
  const_cast<RHomogeneousTreeLikelihood *>(this)->computeTreeDLikelihood(variable);
  return - getDLogLikelihood();
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::computeTreeDLikelihood(const string & variable)
{  
  // Get the node with the branch whose length must be derivated:
  unsigned int brI = TextTools::to<unsigned int>(variable.substr(5));
  const Node * branch = nodes_[brI];
  const Node * father = branch->getFather();
  VVVdouble * _dLikelihoods_father = & _likelihoodData->getDLikelihoodArray(father->getId());
  
  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  unsigned int nbSites  = _dLikelihoods_father->size();
  for(unsigned int i = 0; i < nbSites; i++)
  {
    VVdouble * _dLikelihoods_father_i = & (* _dLikelihoods_father)[i];
    for(unsigned int c = 0; c < nbClasses_; c++)
    {
      Vdouble * _dLikelihoods_father_i_c = & (* _dLikelihoods_father_i)[c];
      for(unsigned int s = 0; s < nbStates_; s++)
      {
        (* _dLikelihoods_father_i_c)[s] = 1.;  
      }
    }
  }

  unsigned int nbNodes = father->getNumberOfSons();
  for(unsigned int l = 0; l < nbNodes; l++)
  {
    const Node * son = father->getSon(l);

    vector <unsigned int> * _patternLinks_father_son = & _likelihoodData->getArrayPositions(father->getId(), son->getId());
    VVVdouble * _likelihoods_son = & _likelihoodData->getLikelihoodArray(son->getId());

    if(son == branch)
    {
      VVVdouble * dpxy__son = & dpxy_[son->getId()];
      for(unsigned int i = 0; i < nbSites; i++)
      {
        VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_father_son)[i]];
        VVdouble * _dLikelihoods_father_i = & (* _dLikelihoods_father)[i];
        for(unsigned int c = 0; c < nbClasses_; c++)
        {
          Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
          Vdouble * _dLikelihoods_father_i_c = & (* _dLikelihoods_father_i)[c];
          VVdouble * dpxy__son_c = & (* dpxy__son)[c];
          for(unsigned int x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble * dpxy__son_c_x = & (* dpxy__son_c)[x];
            for(unsigned int y = 0; y < nbStates_; y++)
            {
              dl += (* dpxy__son_c_x)[y] * (* _likelihoods_son_i_c)[y];
            }
            (* _dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
    else
    {
      VVVdouble * pxy__son = & pxy_[son->getId()];
      for(unsigned int i = 0; i < nbSites; i++)
      {
        VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_father_son)[i]];
        VVdouble * _dLikelihoods_father_i = & (* _dLikelihoods_father)[i];
        for(unsigned int c = 0; c < nbClasses_; c++)
        {
          Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
          Vdouble * _dLikelihoods_father_i_c = & (* _dLikelihoods_father_i)[c];
          VVdouble * pxy__son_c = & (* pxy__son)[c];
          for(unsigned int x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble * pxy__son_c_x = & (* pxy__son_c)[x];
            for(unsigned int y = 0; y < nbStates_; y++)
            {
              dl += (* pxy__son_c_x)[y] * (* _likelihoods_son_i_c)[y];
            }
            (* _dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
  }

  // Now we go down the tree toward the root node:
  computeDownSubtreeDLikelihood(father);
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::computeDownSubtreeDLikelihood(const Node * node)
{
  const Node * father = node->getFather();
  // We assume that the _dLikelihoods array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if(father == NULL) return; // We reached the root!
    
  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  VVVdouble * _dLikelihoods_father = & _likelihoodData->getDLikelihoodArray(father->getId());
  unsigned int nbSites  = _dLikelihoods_father->size();
  for(unsigned int i = 0; i < nbSites; i++)
  {
    VVdouble * _dLikelihoods_father_i = & (* _dLikelihoods_father)[i];
    for(unsigned int c = 0; c < nbClasses_; c++)
    {
      Vdouble * _dLikelihoods_father_i_c = & (* _dLikelihoods_father_i)[c];
      for(unsigned int s = 0; s < nbStates_; s++)
      {
        (* _dLikelihoods_father_i_c)[s] = 1.;  
      }
    }
  }

  unsigned int nbNodes = father->getNumberOfSons();
  for(unsigned int l = 0; l < nbNodes; l++)
  {
    const Node * son = father->getSon(l);

    VVVdouble * pxy__son = & pxy_[son->getId()];
    vector <unsigned int> * _patternLinks_father_son = & _likelihoodData->getArrayPositions(father->getId(),son->getId());
  
    if(son == node)
    {
      VVVdouble * _dLikelihoods_son = & _likelihoodData->getDLikelihoodArray(son->getId());
      for(unsigned int i = 0; i < nbSites; i++)
      {
        VVdouble * _dLikelihoods_son_i = & (* _dLikelihoods_son)[(* _patternLinks_father_son)[i]];
        VVdouble * _dLikelihoods_father_i = & (* _dLikelihoods_father)[i];
        for(unsigned int c = 0; c < nbClasses_; c++)
        {
          Vdouble * _dLikelihoods_son_i_c = & (* _dLikelihoods_son_i)[c];
          Vdouble * _dLikelihoods_father_i_c = & (* _dLikelihoods_father_i)[c];
          VVdouble * pxy__son_c = & (* pxy__son)[c];
          for(unsigned int x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble * pxy__son_c_x = & (* pxy__son_c)[x];
            for(unsigned int y = 0; y < nbStates_; y++)
            {
              dl +=(* pxy__son_c_x)[y] * (* _dLikelihoods_son_i_c)[y];
            }
            (* _dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
    else
    {
      VVVdouble * _likelihoods_son = & _likelihoodData->getLikelihoodArray(son->getId());
      for(unsigned int i = 0; i < nbSites; i++)
      {
        VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_father_son)[i]];
        VVdouble * _dLikelihoods_father_i = & (* _dLikelihoods_father)[i];
        for(unsigned int c = 0; c < nbClasses_; c++)
        {
          Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
          Vdouble * _dLikelihoods_father_i_c = & (* _dLikelihoods_father_i)[c];
          VVdouble * pxy__son_c = & (* pxy__son)[c];
          for(unsigned int x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble * pxy__son_c_x = & (* pxy__son_c)[x];
            for(unsigned int y = 0; y < nbStates_; y++)
            {
              dl += (* pxy__son_c_x)[y] * (* _likelihoods_son_i_c)[y];
            }
            (* _dLikelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
  }

  //Next step: move toward grand father...
  computeDownSubtreeDLikelihood(father);
}
  
/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/  

double RHomogeneousTreeLikelihood::getD2LikelihoodForASiteForARateClass(
  unsigned int site,
  unsigned int rateClass) const
{
  double d2l = 0;
  Vdouble * d2la = & _likelihoodData->getD2LikelihoodArray(_tree->getRootNode()->getId())[_likelihoodData->getRootArrayPosition(site)][rateClass];
  for(unsigned int i = 0; i < nbStates_; i++)
  {
    d2l += (* d2la)[i] * rootFreqs_[i];
  }
  return d2l;
}

/******************************************************************************/  

double RHomogeneousTreeLikelihood::getD2LikelihoodForASite(unsigned int site) const
{
  // Derivative of the sum is the sum of derivatives:
  double d2l = 0;
  for(unsigned int i = 0; i < nbClasses_; i++)
    d2l += getD2LikelihoodForASiteForARateClass(site, i) * _rateDistribution->getProbability(i);
  return d2l;
}

/******************************************************************************/  

double RHomogeneousTreeLikelihood::getD2LogLikelihoodForASite(unsigned int site) const
{
  return getD2LikelihoodForASite(site) / getLikelihoodForASite(site)
  - pow( getDLikelihoodForASite(site) / getLikelihoodForASite(site), 2);
}

/******************************************************************************/  

double RHomogeneousTreeLikelihood::getD2LogLikelihood() const
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for(unsigned int i = 0; i < nbSites_; i++)
    dl += getD2LogLikelihoodForASite(i);
  return dl;
}

/******************************************************************************/

double RHomogeneousTreeLikelihood::getSecondOrderDerivative(const string & variable) const
throw (Exception)
{
  const Parameter * p = &getParameter(variable);
  if(getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameter are not implemented.");
  }
  if(getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }
  
  const_cast<RHomogeneousTreeLikelihood *>(this)->computeTreeD2Likelihood(variable);
  return - getD2LogLikelihood();
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::computeTreeD2Likelihood(const string & variable)
{
  // Get the node with the branch whose length must be derivated:
  unsigned int brI = TextTools::to<unsigned int>(variable.substr(5));
  const Node * branch = nodes_[brI];
  const Node * father = branch->getFather();
  
  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  VVVdouble * _d2Likelihoods_father = & _likelihoodData->getD2LikelihoodArray(father->getId()); 
  unsigned int nbSites  = _d2Likelihoods_father->size();
  for(unsigned int i = 0; i < nbSites; i++)
  {
    VVdouble * _d2Likelihoods_father_i = & (* _d2Likelihoods_father)[i];
    for(unsigned int c = 0; c < nbClasses_; c++)
    {
      Vdouble * _d2Likelihoods_father_i_c = & (* _d2Likelihoods_father_i)[c];
      for(unsigned int s = 0; s < nbStates_; s++)
      {
        (* _d2Likelihoods_father_i_c)[s] = 1.;  
      }
    }
  }

  unsigned int nbNodes = father->getNumberOfSons();
  for(unsigned int l = 0; l < nbNodes; l++)
  {
    const Node * son = father->getSon(l);
    
    vector <unsigned int> * _patternLinks_father_son = & _likelihoodData->getArrayPositions(father->getId(),son->getId());
    VVVdouble * _likelihoods_son = & _likelihoodData->getLikelihoodArray(son->getId());

    if(son == branch)
    {
      VVVdouble * d2pxy__son = & d2pxy_[son->getId()];
      for(unsigned int i = 0; i < nbSites; i++)
      {
        VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_father_son)[i]];
        VVdouble * _d2Likelihoods_father_i = & (* _d2Likelihoods_father)[i];
        for(unsigned int c = 0; c < nbClasses_; c++)
        {
          Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
          Vdouble * _d2Likelihoods_father_i_c = & (* _d2Likelihoods_father_i)[c];
          VVdouble * d2pxy__son_c = & (* d2pxy__son)[c];
          for(unsigned int x = 0; x < nbStates_; x++)
          {
            double d2l = 0;
            Vdouble * d2pxy__son_c_x = & (* d2pxy__son_c)[x];
            for(unsigned int y = 0; y < nbStates_; y++)
            {
              d2l += (* d2pxy__son_c_x)[y] * (* _likelihoods_son_i_c)[y];
            }
            (* _d2Likelihoods_father_i_c)[x] *= d2l;
          }
        }
      }
    }
    else
    {
      VVVdouble * pxy__son = & pxy_[son->getId()];
      for(unsigned int i = 0; i < nbSites; i++)
      {
        VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_father_son)[i]];
        VVdouble * _d2Likelihoods_father_i = & (* _d2Likelihoods_father)[i];
        for(unsigned int c = 0; c < nbClasses_; c++)
        {
          Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
          Vdouble * _d2Likelihoods_father_i_c = & (* _d2Likelihoods_father_i)[c];
          VVdouble * pxy__son_c = & (* pxy__son)[c];
          for(unsigned int x = 0; x < nbStates_; x++)
          {
            double d2l = 0;
            Vdouble * pxy__son_c_x = & (* pxy__son_c)[x];
            for(unsigned int y = 0; y < nbStates_; y++)
            {
              d2l += (* pxy__son_c_x)[y] * (* _likelihoods_son_i_c)[y];
            }
            (* _d2Likelihoods_father_i_c)[x] *= d2l;
          }
        }
      }
    }  
  }

  // Now we go down the tree toward the root node:
  computeDownSubtreeD2Likelihood(father);
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::computeDownSubtreeD2Likelihood(const Node * node)
{
  const Node * father = node->getFather();
  // We assume that the _dLikelihoods array has been filled for the current node 'node'.
  // We will evaluate the array for the father node.
  if(father == NULL) return; // We reached the root!
    
  // Compute dLikelihoods array for the father node.
  // Fist initialize to 1:
  VVVdouble * _d2Likelihoods_father = & _likelihoodData->getD2LikelihoodArray(father->getId());
  unsigned int nbSites  = _d2Likelihoods_father->size();
  for(unsigned int i = 0; i < nbSites; i++)
  {
    VVdouble * _d2Likelihoods_father_i = & (* _d2Likelihoods_father)[i];
    for(unsigned int c = 0; c < nbClasses_; c++)
    {
      Vdouble * _d2Likelihoods_father_i_c = & (* _d2Likelihoods_father_i)[c];
      for(unsigned int s = 0; s < nbStates_; s++)
      {
        (* _d2Likelihoods_father_i_c)[s] = 1.;  
      }
    }
  }

  unsigned int nbNodes = father->getNumberOfSons();
  for(unsigned int l = 0; l < nbNodes; l++)
  {
    const Node * son = father->getSon(l);

    VVVdouble * pxy__son = & pxy_[son->getId()];
    vector <unsigned int> * _patternLinks_father_son = & _likelihoodData->getArrayPositions(father->getId(),son->getId());
  
    if(son == node)
    {
      VVVdouble * _d2Likelihoods_son = & _likelihoodData->getD2LikelihoodArray(son->getId());
      for(unsigned int i = 0; i < nbSites; i++)
      {
        VVdouble * _d2Likelihoods_son_i = & (* _d2Likelihoods_son)[(* _patternLinks_father_son)[i]];
        VVdouble * _d2Likelihoods_father_i = & (* _d2Likelihoods_father)[i];
        for(unsigned int c = 0; c < nbClasses_; c++)
        {
          Vdouble * _d2Likelihoods_son_i_c = & (* _d2Likelihoods_son_i)[c];
          Vdouble * _d2Likelihoods_father_i_c = & (* _d2Likelihoods_father_i)[c];
          VVdouble * pxy__son_c = & (* pxy__son)[c];
          for(unsigned int x = 0; x < nbStates_; x++)
          {
            double d2l = 0;
            Vdouble * pxy__son_c_x = & (* pxy__son_c)[x];
            for(unsigned int y = 0; y < nbStates_; y++)
            {
              d2l += (* pxy__son_c_x)[y] * (* _d2Likelihoods_son_i_c)[y];
            }
            (* _d2Likelihoods_father_i_c)[x] *= d2l;
          }
        }
      }
    }
    else
    {
      VVVdouble * _likelihoods_son = & _likelihoodData->getLikelihoodArray(son->getId());
      for(unsigned int i = 0; i < nbSites; i++)
      {
        VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_father_son)[i]];
        VVdouble * _d2Likelihoods_father_i = & (* _d2Likelihoods_father)[i];
        for(unsigned int c = 0; c < nbClasses_; c++)
        {
          Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
          Vdouble * _d2Likelihoods_father_i_c = & (* _d2Likelihoods_father_i)[c];
          VVdouble * pxy__son_c = & (* pxy__son)[c];
          for(unsigned int x = 0; x < nbStates_; x++)
          {
            double dl = 0;
            Vdouble * pxy__son_c_x = & (* pxy__son_c)[x];
            for(unsigned int y = 0; y < nbStates_; y++)
            {
              dl += (* pxy__son_c_x)[y] * (* _likelihoods_son_i_c)[y];
            }
            (* _d2Likelihoods_father_i_c)[x] *= dl;
          }
        }
      }
    }
  }

  //Next step: move toward grand father...
  computeDownSubtreeD2Likelihood(father);
}

/******************************************************************************/
  
void RHomogeneousTreeLikelihood::computeTreeLikelihood()
{
  computeSubtreeLikelihood(_tree->getRootNode());
}

/******************************************************************************/  

void RHomogeneousTreeLikelihood::computeSubtreeLikelihood(const Node * node)
{
  if(node->isLeaf()) return;

  unsigned int nbSites = _likelihoodData->getLikelihoodArray(node->getId()).size();
  unsigned int nbNodes = node->getNumberOfSons();
    
  // Must reset the likelihood array first (i.e. set all of them to 1):
  VVVdouble * _likelihoods_node = & _likelihoodData->getLikelihoodArray(node->getId());
  for(unsigned int i = 0; i < nbSites; i++)
  {
    //For each site in the sequence,
    VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i];
    for(unsigned int c = 0; c < nbClasses_; c++)
    {
      //For each rate classe,
      Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c];
      for(unsigned int x = 0; x < nbStates_; x++)
      {
        //For each initial state,
        (* _likelihoods_node_i_c)[x] = 1.;
      }
    }
  }

  for(unsigned int l = 0; l < nbNodes; l++)
  {
    //For each son node,  

    const Node * son = node->getSon(l);
    
    computeSubtreeLikelihood(son); //Recursive method:
    
    VVVdouble * pxy__son = & pxy_[son->getId()];
    vector <unsigned int> * _patternLinks_node_son = & _likelihoodData->getArrayPositions(node->getId(), son->getId());
    VVVdouble * _likelihoods_son = & _likelihoodData->getLikelihoodArray(son->getId());

    for(unsigned int i = 0; i < nbSites; i++)
    {
      //For each site in the sequence,
      VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_node_son)[i]];
      VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i];
      for(unsigned int c = 0; c < nbClasses_; c++)
      {
        //For each rate classe,
        Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
        Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c];
        VVdouble * pxy__son_c = & (* pxy__son)[c];
        for(unsigned int x = 0; x < nbStates_; x++)
        {
          //For each initial state,
          Vdouble * pxy__son_c_x = & (* pxy__son_c)[x];
          double likelihood = 0;
          for(unsigned int y = 0; y < nbStates_; y++)
          {
            likelihood += (* pxy__son_c_x)[y] * (* _likelihoods_son_i_c)[y];
          }
          (* _likelihoods_node_i_c)[x] *= likelihood;
        }
      }
    }
  }
}

/******************************************************************************/

void RHomogeneousTreeLikelihood::displayLikelihood(const Node * node)
{
  cout << "Likelihoods at node " << node->getName() << ": " << endl;
  displayLikelihoodArray(_likelihoodData->getLikelihoodArray(node->getId()));
  cout << "                                         ***" << endl;
}

/*******************************************************************************/

