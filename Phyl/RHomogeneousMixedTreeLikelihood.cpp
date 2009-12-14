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
                                                                 const Tree & tree,
                                                                 SubstitutionModel * model,
                                                                 DiscreteDistribution * rDist,
                                                                 bool checkRooted,
                                                                 bool verbose,
                                                                 bool usePatterns
                                                                 )  throw (Exception):
  RHomogeneousTreeLikelihood(tree, model , rDist, checkRooted, verbose)
{
  if ((_mixedmodel=dynamic_cast<MixedModel*>(model_))==NULL)
    throw Exception("Bad model: RHomogeneousMixedTreeLikelihood needs a MixedModel.");
  
  unsigned int s=_mixedmodel->getNumberOfModels();
  for (unsigned int i=0;i<s;i++)
    _treelikelihoodscontainer.push_back(
					new RHomogeneousTreeLikelihood(tree, _mixedmodel->getNModel(i), rDist, checkRooted, false)
					);
  
}

RHomogeneousMixedTreeLikelihood::RHomogeneousMixedTreeLikelihood(
  const Tree & tree,
  const SiteContainer & data,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose,
  bool usePatterns
)

  throw (Exception):
  RHomogeneousTreeLikelihood(tree, model,rDist, checkRooted, verbose)
{
  if ((_mixedmodel=dynamic_cast<MixedModel*>(model_))==NULL)
    throw Exception("Bad model: RHomogeneousMixedTreeLikelihood needs a MixedModel.");
  
  int s=_mixedmodel->getNumberOfModels();
  
  for (int i=0;i<s;i++)
    _treelikelihoodscontainer.push_back(
                                        new RHomogeneousTreeLikelihood(tree, _mixedmodel->getNModel(i), rDist, checkRooted, false)
					);
  setData(data);
}


RHomogeneousMixedTreeLikelihood & RHomogeneousMixedTreeLikelihood:: operator=(const RHomogeneousMixedTreeLikelihood & lik)
{
  RHomogeneousTreeLikelihood::operator=(lik);

  _mixedmodel=lik._mixedmodel->clone();
  
  int s=lik._mixedmodel->getNumberOfModels();
  
  for (int i=0;i<s;i++)
    _treelikelihoodscontainer.push_back(lik._treelikelihoodscontainer[i]->clone());

  return *this;
}


RHomogeneousMixedTreeLikelihood::RHomogeneousMixedTreeLikelihood(const RHomogeneousMixedTreeLikelihood & lik) : RHomogeneousTreeLikelihood(lik), _mixedmodel(lik._mixedmodel->clone())
{
  int s=lik._mixedmodel->getNumberOfModels();
  
  for (int i=0;i<s;i++)
    _treelikelihoodscontainer.push_back(lik._treelikelihoodscontainer[i]->clone());
  
}

RHomogeneousMixedTreeLikelihood::~RHomogeneousMixedTreeLikelihood()
{
  for(unsigned int i=0; i<_treelikelihoodscontainer.size(); i++)
    delete _treelikelihoodscontainer[i];
}


void RHomogeneousMixedTreeLikelihood::initialize() throw(Exception)
{
  for(unsigned int i=0; i<_treelikelihoodscontainer.size(); i++){
    _treelikelihoodscontainer[i]->initialize();
  }
  RHomogeneousTreeLikelihood::initialize();
}

void RHomogeneousMixedTreeLikelihood::setData(const SiteContainer & sites) throw (Exception)
{
  RHomogeneousTreeLikelihood::setData(sites);

  for(unsigned int i=0; i<_treelikelihoodscontainer.size(); i++)
    _treelikelihoodscontainer[i]->setData(sites);
}


double RHomogeneousMixedTreeLikelihood::getLogLikelihood() const
{
  vector<double> reslog;

  for(unsigned int i=0; i<_treelikelihoodscontainer.size(); i++)
    reslog.push_back(_treelikelihoodscontainer[i]->getLogLikelihood());

  return VectorTools::logmeanexp(reslog);
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
  unsigned int i;

  for(i=0; i<_treelikelihoodscontainer.size(); i++)
    reslog.push_back(_treelikelihoodscontainer[i]->getLogLikelihoodForASite(site));

  return VectorTools::logmeanexp(reslog);
}

double RHomogeneousMixedTreeLikelihood::getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
  return exp(getLogLikelihoodForASiteForARateClass(site, rateClass));
}

double RHomogeneousMixedTreeLikelihood::getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const 
{

  vector<double> reslog;
  unsigned int i;

  for(i=0; i<_treelikelihoodscontainer.size(); i++)
    reslog.push_back(_treelikelihoodscontainer[i]->getLogLikelihoodForASiteForARateClass(site, rateClass));

  return VectorTools::logmeanexp(reslog);
}

double RHomogeneousMixedTreeLikelihood::getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
  return exp(getLogLikelihoodForASiteForARateClassForAState(site, rateClass, state));
}

double RHomogeneousMixedTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const 
{
  vector<double> reslog;
  unsigned int i;

  for(i=0; i<_treelikelihoodscontainer.size(); i++)
    reslog.push_back(_treelikelihoodscontainer[i]->getLogLikelihoodForASiteForARateClassForAState(site, rateClass, state));

  return VectorTools::logmeanexp(reslog);
}

void RHomogeneousMixedTreeLikelihood::fireParameterChanged(const ParameterList& params)
{
  RHomogeneousTreeLikelihood::fireParameterChanged(params);
  
  AbstractSubstitutionModel *pm;
  for (unsigned int i=0; i<_treelikelihoodscontainer.size(); i++){
    pm=_mixedmodel->getNModel(i);
    _treelikelihoodscontainer[i]->matchParametersValues(pm->getParameters());
    _treelikelihoodscontainer[i]->matchParametersValues(getParameters());
  }
}

double RHomogeneousMixedTreeLikelihood::getValue() const
throw (Exception)
{
  for (unsigned int i=0; i<_treelikelihoodscontainer.size(); i++)
    if(!_treelikelihoodscontainer[i]->isInitialized())
      throw Exception("RHomogeneousMixedTreeLikelihood::getValue(). Instance is not initialized.");

  return - getLogLikelihood();
}

double RHomogeneousMixedTreeLikelihood::getFirstOrderDerivative(const string & variable) const throw (Exception)
{
  vector<double> rescontainer;
  unsigned int i;

  for(i=0; i<_treelikelihoodscontainer.size(); i++)
    rescontainer.push_back(_treelikelihoodscontainer[i]->getFirstOrderDerivative(variable));

  return ((double)VectorTools::sum(rescontainer))/rescontainer.size();
}

double RHomogeneousMixedTreeLikelihood::getSecondOrderDerivative(const string & variable) const throw (Exception)
{

  vector<double> rescontainer;
  unsigned int i;

  for(i=0; i<_treelikelihoodscontainer.size(); i++)
    rescontainer.push_back(_treelikelihoodscontainer[i]->getSecondOrderDerivative(variable));

  return ((double)VectorTools::sum(rescontainer))/rescontainer.size();
}

void RHomogeneousMixedTreeLikelihood::computeTreeLikelihood()
{

  for(unsigned int i=0; i<_treelikelihoodscontainer.size(); i++)
    _treelikelihoodscontainer[i]->computeTreeLikelihood();

}

double RHomogeneousMixedTreeLikelihood::getDLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass)const
{
  vector<double> rescontainer;
  unsigned int i;

  for(i=0; i<_treelikelihoodscontainer.size(); i++)
    rescontainer.push_back(_treelikelihoodscontainer[i]->getDLikelihoodForASiteForARateClass(site,rateClass));

  return ((double)VectorTools::sum(rescontainer))/rescontainer.size();
}

double RHomogeneousMixedTreeLikelihood::getDLikelihoodForASite(unsigned int site) const
{
  vector<double> rescontainer;
  unsigned int i;

  for(i=0; i<_treelikelihoodscontainer.size(); i++)
    rescontainer.push_back(_treelikelihoodscontainer[i]->getDLikelihoodForASite(site));

  return ((double)VectorTools::sum(rescontainer))/rescontainer.size();
}

void RHomogeneousMixedTreeLikelihood::computeTreeDLikelihood(const string & variable)
{
  for(unsigned int i=0; i<_treelikelihoodscontainer.size(); i++)
    _treelikelihoodscontainer[i]->computeTreeDLikelihood(variable);
}


double RHomogeneousMixedTreeLikelihood::getD2LikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass)const
{
  vector<double> rescontainer;
  unsigned int i;

  for(i=0; i<_treelikelihoodscontainer.size(); i++)
    rescontainer.push_back(_treelikelihoodscontainer[i]->getD2LikelihoodForASiteForARateClass(site,rateClass));

  return ((double)VectorTools::sum(rescontainer))/rescontainer.size();
}

double RHomogeneousMixedTreeLikelihood::getD2LikelihoodForASite(unsigned int site) const
{
  vector<double> rescontainer;
  unsigned int i;

  for(i=0; i<_treelikelihoodscontainer.size(); i++)
    rescontainer.push_back(_treelikelihoodscontainer[i]->getD2LikelihoodForASite(site));

  return ((double)VectorTools::sum(rescontainer))/rescontainer.size();
}

void RHomogeneousMixedTreeLikelihood::computeTreeD2Likelihood(const string & variable)
{
  for(unsigned int i=0; i<_treelikelihoodscontainer.size(); i++)
    _treelikelihoodscontainer[i]->computeTreeD2Likelihood(variable);
}

void RHomogeneousMixedTreeLikelihood::computeSubtreeLikelihood(const Node * node)
{
  for(unsigned int i=0; i<_treelikelihoodscontainer.size(); i++)
    _treelikelihoodscontainer[i]->computeSubtreeLikelihood(node);
}

void RHomogeneousMixedTreeLikelihood::computeDownSubtreeDLikelihood(const Node * node)
{
  for(unsigned int i=0; i<_treelikelihoodscontainer.size(); i++)
    _treelikelihoodscontainer[i]->computeDownSubtreeDLikelihood(node);
}

void RHomogeneousMixedTreeLikelihood::computeDownSubtreeD2Likelihood(const Node * node)
{
  for(unsigned int i=0; i<_treelikelihoodscontainer.size(); i++)
    _treelikelihoodscontainer[i]->computeDownSubtreeD2Likelihood(node);
}


void RHomogeneousMixedTreeLikelihood::computeAllTransitionProbabilities()
{
  for(unsigned int i=0; i<_treelikelihoodscontainer.size(); i++)
    _treelikelihoodscontainer[i]->computeAllTransitionProbabilities();
}


void RHomogeneousMixedTreeLikelihood::computeTransitionProbabilitiesForNode(const Node * node)
{
  for(unsigned int i=0; i<_treelikelihoodscontainer.size(); i++)
    _treelikelihoodscontainer[i]->computeTransitionProbabilitiesForNode(node);
}


void RHomogeneousMixedTreeLikelihood::displayLikelihood(const Node * node)
{
  for(unsigned int i=0; i<_treelikelihoodscontainer.size(); i++)
    _treelikelihoodscontainer[i]->displayLikelihood(node);
}
