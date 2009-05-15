//
// File: AbstractNonHomogeneousTreeLikelihood.cpp
// Created by: Julien Dutheil
// Created on: Tue Oct 09 16:07 2007
// From file: AbstractHomogeneousTreeLikelihood.cpp
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
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "AbstractNonHomogeneousTreeLikelihood.h"
#include "PatternTools.h"

//From SeqLib:
#include <Seq/SiteTools.h>
#include <Seq/SequenceContainerTools.h>

// From Utils:
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

AbstractNonHomogeneousTreeLikelihood::AbstractNonHomogeneousTreeLikelihood(
  const Tree & tree,
  SubstitutionModelSet * modelSet,
  DiscreteDistribution * rDist,
  bool verbose)
  throw (Exception):
  AbstractDiscreteRatesAcrossSitesTreeLikelihood(rDist, verbose)
{
  _init(tree, modelSet, rDist, verbose);
}

/******************************************************************************/

AbstractNonHomogeneousTreeLikelihood::AbstractNonHomogeneousTreeLikelihood(
    const AbstractNonHomogeneousTreeLikelihood & lik) :
  AbstractDiscreteRatesAcrossSitesTreeLikelihood(lik)
{
  _modelSet        = lik._modelSet;
  _pxy             = lik._pxy;
  _dpxy            = lik._dpxy;
  _d2pxy           = lik._d2pxy;
  _nodes           = _tree->getNodes();
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
  _rootFreqs       = lik._rootFreqs;
  _root1           = lik._root1;
  _root2           = lik._root2;
  //Rebuild nodes index:
  for(unsigned int i = 0; i < _nodes.size(); i++)
  {
    const Node * node = _nodes[i];
    _idToNode[node->getId()] = node;
  }
}

/******************************************************************************/

AbstractNonHomogeneousTreeLikelihood & AbstractNonHomogeneousTreeLikelihood::operator=(
    const AbstractNonHomogeneousTreeLikelihood & lik)
{
  AbstractDiscreteRatesAcrossSitesTreeLikelihood::operator=(lik);
  _modelSet        = lik._modelSet;
  _pxy             = lik._pxy;
  _dpxy            = lik._dpxy;
  _d2pxy           = lik._d2pxy;
  _nodes           = _tree->getNodes();
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
  _rootFreqs       = lik._rootFreqs;
  _root1           = lik._root1;
  _root2           = lik._root2;
  //Rebuild nodes index:
  for(unsigned int i = 0; i < _nodes.size(); i++)
  {
    const Node * node = _nodes[i];
    _idToNode[node->getId()] = node;
  }
  return *this;
}

/******************************************************************************/

AbstractNonHomogeneousTreeLikelihood::~AbstractNonHomogeneousTreeLikelihood()
{
  delete _brLenConstraint;
}

/******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::_init(const Tree & tree,
			SubstitutionModelSet * modelSet,
			DiscreteDistribution * rDist,
			bool verbose) throw (Exception)
{
  TreeTools::checkIds(tree, true);
  _tree = new TreeTemplate<Node>(tree);
  _root1 = _tree->getRootNode()->getSon(0)->getId();
  _root2 = _tree->getRootNode()->getSon(1)->getId();
  _nodes = _tree->getNodes();
  _nodes.pop_back(); //Remove the root node (the last added!).  
  _nbNodes = _nodes.size();
  //Build nodes index:
  for(unsigned int i = 0; i < _nodes.size(); i++)
  {
    const Node * node = _nodes[i];
    _idToNode[node->getId()] = node;
  }
  _nbClasses = _rateDistribution->getNumberOfCategories();

  _verbose = verbose;

  _minimumBrLen = 0.000001;
  _brLenConstraint = new IncludingPositiveReal(_minimumBrLen);
  setSubstitutionModelSet(modelSet);
}

/******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::setSubstitutionModelSet(SubstitutionModelSet * modelSet) throw (Exception)
{
  //Check:
  if(_data)
  {
    if(modelSet->getAlphabet()->getAlphabetType() != _data->getAlphabet()->getAlphabetType())
      throw Exception("AbstractNonHomogeneousTreeLikelihood::setSubstitutionModelSet(). Model alphabet do not match existing data.");
  }

  _modelSet = modelSet;
  
  if(_data)
  {
    if(modelSet->getNumberOfStates() != _modelSet->getNumberOfStates())
      setData(*_data); //Have to reinitialize the whole data structure.
  }
  
  _nbStates = modelSet->getNumberOfStates();

  //Allocate transition probabilities arrays:
  for(unsigned int l = 0; l < _nbNodes; l++)
  {
    //For each son node,i
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

  //We have to reset parameters. If the instance is not initialized, this will be done by the initialize method.
  if(_initialized) initParameters();
}

/******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::initialize() throw (Exception)
{
  if(_initialized) throw Exception("AbstractNonHomogeneousTreeLikelihood::initialize(). Object is already initialized.");
  if(_data == NULL) throw Exception("AbstractNonHomogeneousTreeLikelihood::initialize(). Data are no set.");
  initParameters();
  _initialized = true;
  computeAllTransitionProbabilities();
  fireParameterChanged(getParameters());
}

/******************************************************************************/

ParameterList AbstractNonHomogeneousTreeLikelihood::getBranchLengthsParameters() const
{
  if(!_initialized) throw Exception("AbstractNonHomogeneousTreeLikelihood::getBranchLengthsParameters(). Object is not initialized.");
  return _brLenParameters.getCommonParametersWith(getParameters());
}

/******************************************************************************/

ParameterList AbstractNonHomogeneousTreeLikelihood::getSubstitutionModelParameters() const
{
  if(!_initialized) throw Exception("AbstractNonHomogeneousTreeLikelihood::getSubstitutionModelParameters(). Object is not initialized.");
  return _modelSet->getParameters().getCommonParametersWith(getParameters());
}

/******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::initParameters()
{
  // Reset parameters:
  resetParameters_();
  
  // Branch lengths:
  initBranchLengthsParameters();
  addParameters_(_brLenParameters);
  
  // Substitution model:
  addParameters_(_modelSet->getParameters());
  
  // Rate distribution:
  addParameters_(_rateDistribution->getIndependentParameters());
}

/******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::applyParameters() throw (Exception)
{
  if(!_initialized) throw Exception("AbstractNonHomogeneousTreeLikelihood::applyParameters(). Object not initialized.");
  //Apply branch lengths:
  for(unsigned int i = 0; i < _nbNodes; i++)
  {
    int id = _nodes[i]->getId();
    if(id == _root1)
    {
      const Parameter * rootBrLen = &getParameter("BrLenRoot");
      const Parameter * rootPos = &getParameter("RootPosition");
      _nodes[i]->setDistanceToFather(rootBrLen->getValue() * rootPos->getValue());
    }
    else if(id == _root2)
    {
      const Parameter * rootBrLen = &getParameter("BrLenRoot");
      const Parameter * rootPos = &getParameter("RootPosition");
      _nodes[i]->setDistanceToFather(rootBrLen->getValue() * (1. - rootPos->getValue()));
    }
    else
    {
      const Parameter * brLen = &getParameter(string("BrLen") + TextTools::toString(i));
      if(brLen != NULL) _nodes[i]->setDistanceToFather(brLen->getValue());
    }
  }
  //Apply substitution model parameters:
  _modelSet->matchParametersValues(getParameters());
  //Apply rate distribution parameters:
  _rateDistribution->matchParametersValues(getParameters());
}

/******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::initBranchLengthsParameters()
{
  _brLenParameters.reset();
  double l1 = 0, l2 = 0;
  for(unsigned int i = 0; i < _nbNodes; i++)
  {
    double d = _minimumBrLen;
    if(!_nodes[i]->hasDistanceToFather())
    {
      ApplicationTools::displayWarning("Missing branch length " + TextTools::toString(i) + ". Value is set to " + TextTools::toString(_minimumBrLen));
      _nodes[i]->setDistanceToFather(_minimumBrLen);
    }
    else
    {
      d = _nodes[i]->getDistanceToFather();
      if (d < _minimumBrLen)
      {
        ApplicationTools::displayWarning("Branch length " + TextTools::toString(i) + " is too small: " + TextTools::toString(d) + ". Value is set to " + TextTools::toString(_minimumBrLen));
        _nodes[i]->setDistanceToFather(_minimumBrLen);
        d = _minimumBrLen;
      }
    }
    if(_nodes[i]->getId() == _root1)
      l1 = d;
    else if(_nodes[i]->getId() == _root2)
      l2 = d;
    else
    {
      _brLenParameters.addParameter(Parameter("BrLen" + TextTools::toString(i), d, _brLenConstraint->clone(), true)); //Attach constraint to avoid clonage problems!
    }
  }
  _brLenParameters.addParameter(Parameter("BrLenRoot", l1 + l2, _brLenConstraint->clone(), true));
  _brLenParameters.addParameter(Parameter("RootPosition", l1 / (l1 + l2), &Parameter::PROP_CONSTRAINT_IN));
}

/*******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::computeAllTransitionProbabilities()
{
  for(unsigned int l = 0; l < _nbNodes; l++)
  {
    //For each node,
    Node * node = _nodes[l];
    computeTransitionProbabilitiesForNode(node);
  }
  _rootFreqs = _modelSet->getRootFrequencies();
}

/*******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::computeTransitionProbabilitiesForNode(const Node * node)
{
  const SubstitutionModel * model = _modelSet->getModelForNode(node->getId());
  double l = node->getDistanceToFather(); 

  //Computes all pxy and pyx once for all:
  VVVdouble * _pxy_node = & _pxy[node->getId()];
  for(unsigned int c = 0; c < _nbClasses; c++)
  {
    VVdouble * _pxy_node_c = & (* _pxy_node)[c];
    RowMatrix<double> Q = model->getPij_t(l * _rateDistribution->getCategory(c));
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
      RowMatrix<double> dQ = model->getdPij_dt(l * rc);  
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
      RowMatrix<double> d2Q = model->getd2Pij_dt2(l * rc);
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

