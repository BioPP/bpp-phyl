//
// File: AbstractHomogeneousTreeLikelihood.cpp
// Created by: Julien Dutheil
// Created on: Thr Dec 23 12:03 2004
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

#include "AbstractHomogeneousTreeLikelihood.h"
#include "PatternTools.h"

//From SeqLib:
#include <Seq/SiteTools.h>
#include <Seq/SequenceContainerTools.h>

// From Utils:
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>

// From the STL:
#include <iostream>
using namespace std;

/******************************************************************************/

AbstractHomogeneousTreeLikelihood::AbstractHomogeneousTreeLikelihood(
  const Tree & tree,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose)
  throw (Exception):
  AbstractDiscreteRatesAcrossSitesTreeLikelihood(rDist, verbose)
{
  init(tree, model, rDist, checkRooted, verbose);
}

/******************************************************************************/

AbstractHomogeneousTreeLikelihood::AbstractHomogeneousTreeLikelihood(
    const AbstractHomogeneousTreeLikelihood & lik) :
  //AbstractParametrizable(lik),
  //AbstractTreeLikelihood(lik),
  AbstractDiscreteRatesAcrossSitesTreeLikelihood(lik)
{
  _model           = lik._model;
  _pxy             = lik._pxy;
  _dpxy            = lik._dpxy;
  _d2pxy           = lik._d2pxy;
  _nodes = _tree->getNodes();
  _nodes.pop_back(); //Remove the root node (the last added!).  
	_nbSites         = lik._nbSites;
  _nbDistinctSites = lik._nbDistinctSites;
	_nbClasses       = lik._nbClasses;
	_nbStates        = lik._nbStates;
	_nbNodes         = lik._nbNodes;
  _verbose         = lik._verbose;
  _minimumBrLen    = lik._minimumBrLen;
  _brLenParameters = lik._brLenParameters;
  _brLenConstraint = lik._brLenConstraint->clone();
  //Update the constraints on branch lengths:
  for(unsigned int i = 0; i < _brLenParameters.size(); i++)
  {
    _brLenParameters[i]->setConstraint(_brLenConstraint);
  }
}

/******************************************************************************/

AbstractHomogeneousTreeLikelihood & AbstractHomogeneousTreeLikelihood::operator=(
    const AbstractHomogeneousTreeLikelihood & lik)
{
  AbstractDiscreteRatesAcrossSitesTreeLikelihood::operator=(lik);
  _model           = lik._model;
  _pxy             = lik._pxy;
  _dpxy            = lik._dpxy;
  _d2pxy           = lik._d2pxy;
  _nodes = _tree->getNodes();
  _nodes.pop_back(); //Remove the root node (the last added!).  
	_nbSites         = lik._nbSites;
  _nbDistinctSites = lik._nbDistinctSites;
	_nbClasses       = lik._nbClasses;
	_nbStates        = lik._nbStates;
	_nbNodes         = lik._nbNodes;
  _verbose         = lik._verbose;
  _minimumBrLen    = lik._minimumBrLen;
  _brLenParameters = lik._brLenParameters;
  _brLenConstraint = lik._brLenConstraint->clone();
  //Update the constraints on branch lengths:
  for(unsigned int i = 0; i < _brLenParameters.size(); i++)
  {
    _brLenParameters[i]->setConstraint(_brLenConstraint);
  }
  return *this;
}

/******************************************************************************/

AbstractHomogeneousTreeLikelihood::~AbstractHomogeneousTreeLikelihood()
{
  delete _brLenConstraint;
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::init(const Tree & tree,
			SubstitutionModel * model,
			DiscreteDistribution * rDist,
      bool checkRooted,
			bool verbose) throw (Exception)
{
  _tree = new TreeTemplate<Node>(tree);
  if(checkRooted && _tree->isRooted())
  {
    if(verbose) ApplicationTools::displayWarning("Tree has been unrooted.");
    _tree->unroot();
    _tree->resetNodesId();
  }
  _nodes = _tree->getNodes();
  _nodes.pop_back(); //Remove the root node (the last added!).  
  _nbNodes = _nodes.size();
  
  _model = model;

  _nbStates = model->getNumberOfStates();
  
  _nbClasses = _rateDistribution->getNumberOfCategories();

  //Allocate transition probabilities arrays:
  for(unsigned int l = 0; l < _nbNodes; l++)
  {
    //For each son node,
    Node * son = _nodes[l];

    VVVdouble * _pxy_son = & _pxy[son->getId()];
    _pxy_son->resize(_nbClasses);
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      VVdouble * _pxy_son_c = & (* _pxy_son)[c];
      _pxy_son_c->resize(_nbStates);
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        (*_pxy_son_c)[x].resize(_nbStates);
      }
    }
  
    VVVdouble * _dpxy_son = & _dpxy[son->getId()];
    _dpxy_son->resize(_nbClasses);
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      VVdouble * _dpxy_son_c = & (* _dpxy_son)[c];
      _dpxy_son_c->resize(_nbStates);
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        (* _dpxy_son_c)[x].resize(_nbStates);
      }
    }
      
    VVVdouble * _d2pxy_son = & _d2pxy[son->getId()];
    _d2pxy_son->resize(_nbClasses);
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      VVdouble * _d2pxy_son_c = & (* _d2pxy_son)[c];
      _d2pxy_son_c->resize(_nbStates);
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        (* _d2pxy_son_c)[x].resize(_nbStates);
      }
    }
  }

  _verbose = verbose;

  _minimumBrLen = 0.000001;
  _brLenConstraint = new IncludingPositiveReal(_minimumBrLen);
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::initialize() throw (Exception)
{
  if(_initialized) throw Exception("AbstractHomogeneousTreeLikelihood::initialize(). Object is already initialized.");
  if(_data == NULL) throw Exception("AbstractHomogeneousTreeLikelihood::initialize(). Data are no set.");
  initParameters();
  _initialized = true;
  fireParameterChanged(_parameters);
}

/******************************************************************************/

ParameterList AbstractHomogeneousTreeLikelihood::getBranchLengthsParameters() const
{
  if(!_initialized) throw Exception("AbstractHomogeneousTreeLikelihood::getBranchLengthsParameters(). Object is not initialized.");
  return _brLenParameters.getCommonParametersWith(_parameters);
}

/******************************************************************************/

ParameterList AbstractHomogeneousTreeLikelihood::getSubstitutionModelParameters() const
{
  if(!_initialized) throw Exception("AbstractHomogeneousTreeLikelihood::getSubstitutionModelParameters(). Object is not initialized.");
  return _model->getParameters().getCommonParametersWith(_parameters);
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::initParameters()
{
  // Reset parameters:
  _parameters.reset();
  
  // Branch lengths:
  initBranchLengthsParameters();
  _parameters.addParameters(_brLenParameters);
  
  // Substitution model:
  _parameters.addParameters(_model->getParameters());
  
  // Rate distribution:
  _parameters.addParameters(_rateDistribution->getParameters());
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::ignoreParameter(const string & name)
throw (ParameterNotFoundException)
{
  _parameters.deleteParameter(name);
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::applyParameters() throw (Exception)
{
  if(!_initialized) throw Exception("AbstractHomogeneousTreeLikelihood::applyParameters(). Object not initialized.");
  //Apply branch lengths:
  for(unsigned int i = 0; i < _nbNodes; i++)
  {
    const Parameter * brLen = _parameters.getParameter(string("BrLen") + TextTools::toString(_nodes[i]->getId()));
    if(brLen != NULL) _nodes[i]->setDistanceToFather(brLen->getValue());
  }
  //Apply substitution model parameters:
  _model->matchParametersValues(_parameters);
  //Apply rate distribution parameters:
  _rateDistribution->matchParametersValues(_parameters);
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::initBranchLengthsParameters()
{
  _brLenParameters.reset();
  for(unsigned int i = 0; i < _nbNodes; i++)
  {
    double d = _minimumBrLen;
    if(!_nodes[i]->hasDistanceToFather())
    {
      cout << "WARNING!!! Missing branch length " << i << ". Value is set to " << _minimumBrLen << endl;
      _nodes[i]->setDistanceToFather(_minimumBrLen);
    }
    else
    {
      d = _nodes[i]->getDistanceToFather();
      if (d < _minimumBrLen)
      {
        cout << "WARNING!!! Branch length " << i << " is too small: " << d << ". Value is set to " << _minimumBrLen << endl;
        _nodes[i]->setDistanceToFather(_minimumBrLen);
        d = _minimumBrLen;
      }
    }
    _brLenParameters.addParameter(Parameter("BrLen" + TextTools::toString(_nodes[i]->getId()), d, _brLenConstraint));
  }
}

/*******************************************************************************/

void AbstractHomogeneousTreeLikelihood::computeAllTransitionProbabilities()
{
  for(unsigned int l = 0; l < _nbNodes; l++)
  {
    //For each node,
    Node * node = _nodes[l];
    computeTransitionProbabilitiesForNode(node);
  }
}

/*******************************************************************************/

void AbstractHomogeneousTreeLikelihood::computeTransitionProbabilitiesForNode(const Node * node)
{
  double l = node->getDistanceToFather(); 

  //Computes all pxy and pyx once for all:
  VVVdouble * _pxy_node = & _pxy[node->getId()];
  for(unsigned int c = 0; c < _nbClasses; c++)
  {
    VVdouble * _pxy_node_c = & (* _pxy_node)[c];
    RowMatrix<double> Q = _model->getPij_t(l * _rateDistribution->getCategory(c));
    for(unsigned int x = 0; x < _nbStates; x++)
    {
      Vdouble * _pxy_node_c_x = & (* _pxy_node_c)[x];
      for(unsigned int y = 0; y < _nbStates; y++)
      {
        (* _pxy_node_c_x)[y] = Q(x, y);
      }
    }
  }
  
  if(_computeFirstOrderDerivatives)
  {
    //Computes all dpxy/dt once for all:
    VVVdouble * _dpxy_node = & _dpxy[node->getId()];
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      VVdouble * _dpxy_node_c = & (* _dpxy_node)[c];
      double rc = _rateDistribution->getCategory(c);
      RowMatrix<double> dQ = _model->getdPij_dt(l * rc);  
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        Vdouble * _dpxy_node_c_x = & (* _dpxy_node_c)[x];
        for(unsigned int y = 0; y < _nbStates; y++)
        {
          (* _dpxy_node_c_x)[y] = rc * dQ(x, y); 
        }
      }
    }
  }
      
  if(_computeSecondOrderDerivatives)
  {
    //Computes all d2pxy/dt2 once for all:
    VVVdouble * _d2pxy_node = & _d2pxy[node->getId()];
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      VVdouble * _d2pxy_node_c = & (* _d2pxy_node)[c];
      double rc =  _rateDistribution->getCategory(c);
      RowMatrix<double> d2Q = _model->getd2Pij_dt2(l * rc);
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        Vdouble * _d2pxy_node_c_x = & (* _d2pxy_node_c)[x];
        for(unsigned int y = 0; y < _nbStates; y++)
        {
          (* _d2pxy_node_c_x)[y] = rc * rc * d2Q(x, y);
        }
      }
    }
  }
}

/*******************************************************************************/

