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
	TreeTemplate<Node> & tree,
	const SiteContainer & data,
	SubstitutionModel * model,
	DiscreteDistribution * rDist,
  bool checkRooted,
	bool verbose)
	throw (Exception):
	AbstractDiscreteRatesAcrossSitesTreeLikelihood(rDist, verbose)
{
	_tree = &tree;
	if(checkRooted && _tree -> isRooted()) {
		if(verbose) ApplicationTools::displayWarning("Tree has been unrooted.");
		_tree -> unroot();
	}
	//Sequences will be in the same order than in the tree:
	_data = PatternTools::getSequenceSubset(data, * _tree -> getRootNode());
	if(_data -> getNumberOfSequences() == 1) throw Exception("Error, only 1 sequence!");
	if(_data -> getNumberOfSequences() == 0) throw Exception("Error, no sequence!");
	_model = model;
	if(_data -> getAlphabet() -> getAlphabetType()
			!= _model -> getAlphabet() -> getAlphabetType())
		throw AlphabetMismatchException("AbstractHomogeneousTreeLikelihood::AbstractHomogeneousTreeLikelihood. Data and model must have the same alphabet type.",
				_data -> getAlphabet(),
				_model -> getAlphabet());
	_alphabet = _data -> getAlphabet();

	_nodes = _tree -> getNodes();
	
	_nodes.pop_back(); //Remove the root node (the last added!).
	
	_nbSites   = _data -> getNumberOfSites();
	_nbClasses = _rateDistribution -> getNumberOfCategories();
	_nbStates  = _alphabet -> getSize();
	_nbNodes   = _nodes.size();
	
}

/******************************************************************************/

AbstractHomogeneousTreeLikelihood::~AbstractHomogeneousTreeLikelihood()
{
	delete _data; 
}

/******************************************************************************/

ParameterList AbstractHomogeneousTreeLikelihood::getBranchLengthsParameters() const
{
	return _brLenParameters.getCommonParametersWith(_parameters);
}

/******************************************************************************/

ParameterList AbstractHomogeneousTreeLikelihood::getSubstitutionModelParameters() const
{
	return _model -> getParameters().getCommonParametersWith(_parameters);
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
    double d = 0;
    if(!_nodes[i] -> hasDistanceToFather()) {
			cout << "WARNING!!! Missing branch length " << i << ". Value is set to 0." << endl;
			_nodes[i] -> setDistanceToFather(0.);
    } else {
  		d = _nodes[i] -> getDistanceToFather();
	  	if (d < 0) {
		  	cout << "WARNING!!! Branch length " << i << " is <0. Value is set to 0." << endl;
			  _nodes[i] -> setDistanceToFather(0.);
			  d = 0.;
		  }
    }
		_brLenParameters.addParameter(Parameter("BrLen" + TextTools::toString(i), d, & Parameter::R_PLUS));
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
