//
// File: AbstractHomogeneousTreeLikelihood.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Thr Dec 23 12:03 2004
//

#include "AbstractHomogeneousTreeLikelihood.h"
#include "PatternTools.h"
#include "ApplicationTools.h"

//From SeqLib:
#include <Seq/SiteTools.h>
#include <Seq/SequenceContainerTools.h>

// From Utils:
#include <Utils/TextTools.h>

// From the STL:
#include <iostream>
using namespace std;

/******************************************************************************/

AbstractHomogeneousTreeLikelihood::AbstractHomogeneousTreeLikelihood(
	Tree<Node> & tree,
	const SiteContainer & data,
	SubstitutionModel * model,
	DiscreteDistribution * rDist,
	bool verbose)
	throw (Exception):
	AbstractTreeLikelihood()
{
	_tree = &tree;
	if(_tree -> isRooted()) {
		if(verbose) ApplicationTools::displayWarning("Tree has been unrooted.");
		_tree -> unroot();
	}
	//Sequences will be in the same order than in the tree:
	_data = PatternTools::getSequenceSubset(data, * _tree -> getRootNode());
	_model = model;
	if(_data -> getAlphabet() -> getAlphabetType()
			!= _model -> getAlphabet() -> getAlphabetType())
		throw AlphabetMismatchException("AbstractHomogeneousTreeLikelihood::AbstractHomogeneousTreeLikelihood. Data and model must have the same alphabet type.",
				_data -> getAlphabet(),
				_model -> getAlphabet());

	_rateDistribution = rDist;
	
	_nodes = _tree -> getNodes();
	
	_nodes.pop_back(); //Remove the root node (the last added!).
	
	_nbSites   = _data -> getNumberOfSites();
	_nbClasses = _rateDistribution -> getNumberOfCategories();
	_nbStates  = _model -> getAlphabet() -> getSize();
	_nbNodes   = _nodes.size();
	
}

/******************************************************************************/

AbstractHomogeneousTreeLikelihood::~AbstractHomogeneousTreeLikelihood()
{
	delete _data; 
}

/******************************************************************************/

ParameterList AbstractHomogeneousTreeLikelihood::getBranchLengthsParameters() const {
	return _brLenParameters.getCommonParametersWith(_parameters);
}

/******************************************************************************/

ParameterList AbstractHomogeneousTreeLikelihood::getSubstitutionModelParameters() const {
	return _model -> getParameters().getCommonParametersWith(_parameters);
}

/******************************************************************************/

ParameterList AbstractHomogeneousTreeLikelihood::getRateDistributionParameters() const {
	return _rateDistribution -> getParameters().getCommonParametersWith(_parameters);
}

/******************************************************************************/

VVdouble AbstractHomogeneousTreeLikelihood::getLikelihoodForEachSiteForEachRateClass() const
{
	VVdouble l(_nbSites);
	for(unsigned int i = 0; i < _nbSites; i++) {
		l[i].resize(_nbClasses);
		for(unsigned int j = 0; j < _nbClasses; j++)
			l[i][j] = getLikelihoodForASiteForARateClass(i, j);
	}
	return l;
}

/******************************************************************************/

VVdouble AbstractHomogeneousTreeLikelihood::getLogLikelihoodForEachSiteForEachRateClass() const
{
	VVdouble l(_nbSites);
	for(unsigned int i = 0; i < _nbSites; i++) {
		l[i] = Vdouble(_nbClasses);
		for(unsigned int j = 0; j < _nbClasses; j++)
			l[i][j] = getLogLikelihoodForASiteForARateClass(i, j);
	}
	return l;
}

/******************************************************************************/

VVVdouble AbstractHomogeneousTreeLikelihood::getLikelihoodForEachSiteForEachRateClassForEachState() const
{
	VVVdouble l(_nbSites);
	for(unsigned int i = 0; i < _nbSites; i++) {
		l[i].resize(_nbClasses);
		for(unsigned int j = 0; j < _nbClasses; j++) {
			l[i][j].resize(_nbStates);
			for(int x = 0; x < _nbStates; x++) {
				l[i][j][x] = getLikelihoodForASiteForARateClassForAState(i, j, x);
			}
		}
	}
	return l;
}

/******************************************************************************/

VVVdouble AbstractHomogeneousTreeLikelihood::getLogLikelihoodForEachSiteForEachRateClassForEachState() const
{
	VVVdouble l(_nbSites);
	for(unsigned int i = 0; i < _nbSites; i++) {
		l[i].resize(_nbClasses);
		for(unsigned int j = 0; j < _nbClasses; j++) {
			l[i][j].resize(_nbStates);
			for(int x = 0; x < _nbStates; x++) {
				l[i][j][x] = getLogLikelihoodForASiteForARateClassForAState(i, j, x);
			}
		}
	}
	return l;
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
	_parameters.addParameters(_model -> getParameters());
	
	// Rate distribution:
	_parameters.addParameters(_rateDistribution -> getParameters());
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
	//Apply branch lengths:
	for(unsigned int i = 0; i < _nbNodes; i++) {
		const Parameter * brLen = _parameters.getParameter(string("BrLen") + TextTools::toString(i));
		_nodes[i] -> setDistanceToFather(brLen -> getValue());
	}
	//Apply substitution model parameters:
	_model -> matchParametersValues(_parameters);
	//Apply rate distribution parameters:
	_rateDistribution -> matchParametersValues(_parameters);
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::initBranchLengthsParameters()
{
	for(unsigned int i = 0; i < _nbNodes; i++) {
		double d = _nodes[i] -> getDistanceToFather();
		if (d <= 0) {
			cout << "WARNING!!! Branch length " << i << " is <=0. Value is set to 0.000001." << endl;
			_nodes[i] -> setDistanceToFather(0.000001);
			d = 0.000001;
		}
		_brLenParameters.addParameter(Parameter("BrLen" + TextTools::toString(i), d, & Parameter::R_PLUS_STAR));
	}
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::resetLikelihoodArray(VVVdouble & likelihoodArray)
{
	unsigned int nbSites   = likelihoodArray.size();
	unsigned int nbClasses = likelihoodArray[0].size();
	unsigned int nbStates  = likelihoodArray[0][0].size();
	for(unsigned int i = 0; i < nbSites; i++) {
		for(unsigned int c = 0; c < nbClasses; c++) {
			for(unsigned int s = 0; s < nbStates; s++) {
				likelihoodArray[i][c][s] = 1.;
			}
		}
	}
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::displayLikelihoodArray(const VVVdouble & likelihoodArray)
{
	unsigned int nbSites   = likelihoodArray.size();
	unsigned int nbClasses = likelihoodArray[0].size();
	unsigned int nbStates  = likelihoodArray[0][0].size();
	for(unsigned int i = 0; i < nbSites; i++) {
		cout << "Site " << i << ":" << endl;
		for(unsigned int c = 0; c < nbClasses; c++) {
			cout << "Rate class " << c;
			for(unsigned int s = 0; s < nbStates; s++) {
				cout << "\t" << likelihoodArray[i][c][s];
			}
			cout << endl;
		}
		cout << endl;
	}
}

/*******************************************************************************/
