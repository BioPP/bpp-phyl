//
// File: HomogeneousTreeLikelihood.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Fri Oct 17 18:14:51 2003
//

#include "DRHomogeneousTreeLikelihood.h"
#include "PatternTools.h"
#include "ApplicationTools.h"

//From SeqLib:
#include <Seq/SiteTools.h>
#include <Seq/AlignedSequenceContainer.h>
#include <Seq/SequenceContainerTools.h>

// From Utils:
#include <Utils/TextTools.h>

// From NumCalc:
using namespace VectorFunctions;

// From the STL:
#include <iostream>
using namespace std;

/******************************************************************************/

DRHomogeneousTreeLikelihood::DRHomogeneousTreeLikelihood(
	Tree & tree,
	const SiteContainer & data,
	SubstitutionModel * model,
	DiscreteDistribution * rDist,
	bool verbose)
	throw (Exception):
	AbstractHomogeneousTreeLikelihood(tree, data, model, rDist, verbose)
{
	if(verbose) ApplicationTools::message << "Double-Recursive Homogeneous Tree Likelihood" << endl;	
	
	//Initialize root patterns:
	_shrunkData = PatternTools::shrinkSiteSet(* _data);
	_nbDistinctSites = _shrunkData -> getNumberOfSites();
	if(verbose) ApplicationTools::displayResult("Number of distinct sites",
			TextTools::toString(_nbDistinctSites));
	
	if(verbose) ApplicationTools::displayTask("Init root patterns");
	_rootPatternLinks.resize(_nbSites);
	for(unsigned int i = 0; i < _nbSites; i++) {
		const Site * site1 =  _data -> getSite(i);
		for(unsigned int ii = 0; ii < _nbDistinctSites; ii++) {
			if(SiteTools::areSitesIdentical(* _shrunkData -> getSite(ii), * site1)) {
				_rootPatternLinks[i] = ii;
				break;
			}
		}
	}
	if(verbose) ApplicationTools::displayTaskDone();
	
	//Init _likelihoods:
	if(verbose) ApplicationTools::displayTask("Init likelihoods arrays recursively");
	// Clone data for more efficiency on sequences access:
	const SiteContainer * sequences = new AlignedSequenceContainer(* _shrunkData);
	initTreeLikelihoods(_tree -> getRootNode(), * sequences);
	delete sequences;

	// Now initialize root likelihoods and derivatives:
	_rootLikelihoods.resize(_nbDistinctSites);
	for(unsigned int i = 0; i < _nbDistinctSites; i++) {
		VVdouble * _rootLikelihoods_i = & _rootLikelihoods[i];
		_rootLikelihoods_i -> resize(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			Vdouble * _rootLikelihoods_i_c = & (* _rootLikelihoods_i)[c];
			_rootLikelihoods_i_c -> resize(_nbStates);
			for(unsigned int x = 0; x < _nbStates; x++) {
				(* _rootLikelihoods_i_c)[x] = 1.;
			}
		}
	}

	if(verbose) ApplicationTools::displayTaskDone();
	
	// Now initializes all parameters:
	initParameters();
	fireParameterChanged(_parameters);
}

/******************************************************************************/

DRHomogeneousTreeLikelihood::~DRHomogeneousTreeLikelihood() {
	delete _data; 
	delete _shrunkData;
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLikelihood() const
{
	double l = 1.;
	for(unsigned int i = 0; i < _nbSites; i++) {
		l *= getLikelihoodForASite(i);
	}
	return l;
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihood() const
{
	double ll = 0;
	for(unsigned int i = 0; i < _nbSites; i++) {
		ll += getLogLikelihoodForASite(i);
	}
	return ll;
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLikelihoodForASite(unsigned int site) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbClasses; i++) {
		l += getLikelihoodForASiteForARateClass(site, i) * _rateDistribution -> getProbability(i);
	}
	return l;
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihoodForASite(unsigned int site) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbClasses; i++) {
		l += getLikelihoodForASiteForARateClass(site, i) * _rateDistribution -> getProbability(i);
	}
	//if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
	return log(l);
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLikelihoodForASiteForAState(unsigned int site, int state) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbClasses; i++) {
		l += getLikelihoodForASiteForARateClassForAState(site, i, state) * _rateDistribution -> getProbability(i);
	}
	return l;
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihoodForASiteForAState(unsigned int site, int state) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbClasses; i++) {
		l += getLikelihoodForASiteForARateClassForAState(site, i, state) * _rateDistribution -> getProbability(i);
	}
	//if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
	return log(l);
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbStates; i++) {
		l += _rootLikelihoods[_rootPatternLinks[site]][rateClass][i] * _model -> freq(i);
	}
	return l;
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbStates; i++) {
		l += _rootLikelihoods[_rootPatternLinks[site]][rateClass][i] * _model -> freq(i);
	}
	//if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
	return log(l);
}

/******************************************************************************/	

double DRHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
	return _rootLikelihoods[_rootPatternLinks[site]][rateClass][state];
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
	return log(_rootLikelihoods[_rootPatternLinks[site]][rateClass][state]);
}

/******************************************************************************/	


VVdouble DRHomogeneousTreeLikelihood::getPosteriorProbabilitiesOfEachRate() const
{
	VVdouble pb = getLikelihoodForEachSiteForEachRateClass();
	Vdouble  l  = getLikelihoodForEachSite();
	for(unsigned int i = 0; i < _nbSites; i++) {
		for(unsigned int j = 0; j < _nbClasses; j++) pb[i][j] = pb[i][j] * _rateDistribution -> getProbability(j) / l[i]; 
	}
	return pb;
}
	
/******************************************************************************/

Vint DRHomogeneousTreeLikelihood::getRateClassWithMaxPostProbOfEachSite() const
{
	VVdouble l = getLikelihoodForEachSiteForEachRateClass();
	Vint classes(_nbSites);
	for(unsigned int i = 0; i < _nbSites; i++) classes[i] = posmax<double>(l[i]);
	return classes;
}

/******************************************************************************/

Vdouble DRHomogeneousTreeLikelihood::getRateWithMaxPostProbOfEachSite() const
{
	VVdouble l = getLikelihoodForEachSiteForEachRateClass();
	Vdouble rates(_nbSites);
	for(unsigned int i = 0; i < _nbSites; i++) {
		rates[i] = _rateDistribution -> getCategory(posmax<double>(l[i]));
	}
	return rates;
}

/******************************************************************************/

Vdouble DRHomogeneousTreeLikelihood::getPosteriorRateOfEachSite() const
{
	VVdouble lr = getLikelihoodForEachSiteForEachRateClass();
	Vdouble  l  = getLikelihoodForEachSite();
	Vdouble rates(_nbSites, 0.);
	for(unsigned int i = 0; i < _nbSites; i++) {
		for(unsigned int j = 0; j < _nbClasses; j++) {
			rates[i] += (lr[i][j] / l[i]) * _rateDistribution -> getProbability(j) *  _rateDistribution -> getCategory(j);
		}
	}
	return rates;
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::setParameters(const ParameterList & parameters)
	throw (ParameterNotFoundException, ConstraintException)
{
	setParametersValues(parameters);
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::fireParameterChanged(const ParameterList & params)
{
	applyParameters();

	// For now we ignore the parameter that changed and we recompute all arrays...
	for(unsigned int l = 0; l < _nbNodes; l++) {
		//For each son node,
		Node * son = _nodes[l];
		double l = son -> getDistanceToFather(); 
	
		//Computes all pxy and pyx once for all:
		VVVdouble * _pxy_son = & _pxy[son];
		_pxy_son -> resize(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			VVdouble * _pxy_son_c = & (* _pxy_son)[c];
			_pxy_son_c -> resize(_nbStates);
			Matrix Q = _model -> getPij_t(l * _rateDistribution -> getCategory(c));
			for(unsigned int x = 0; x < _nbStates; x++) {
				Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
				_pxy_son_c_x -> resize(_nbStates);
				for(unsigned int y = 0; y < _nbStates; y++) {
					(* _pxy_son_c_x)[y] = Q(x, y);
				}
			}
		}

		//Computes all dpxy/dt once for all:
		VVVdouble * _dpxy_son = & _dpxy[son];
		_dpxy_son -> resize(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			VVdouble * _dpxy_son_c = & (* _dpxy_son)[c];
			_dpxy_son_c -> resize(_nbStates);
			double rc = _rateDistribution -> getCategory(c);
			Matrix dQ = _model -> getdPij_dt(l * rc);  
			for(unsigned int x = 0; x < _nbStates; x++) {
				Vdouble * _dpxy_son_c_x = & (* _dpxy_son_c)[x];
				_dpxy_son_c_x -> resize(_nbStates);
				for(unsigned int y = 0; y < _nbStates; y++) {
					(* _dpxy_son_c_x)[y] =  rc * dQ(x, y); 
				}
			}
		}
		
		//Computes all d2pxy/dt2 once for all:
		VVVdouble * _d2pxy_son = & _d2pxy[son];
		_d2pxy_son -> resize(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			VVdouble * _d2pxy_son_c = & (* _d2pxy_son)[c];
			_d2pxy_son_c -> resize(_nbStates);
			double rc =  _rateDistribution -> getCategory(c);
			Matrix d2Q = _model -> getd2Pij_dt2(l * rc);
			for(unsigned int x = 0; x < _nbStates; x++) {
				Vdouble * _d2pxy_son_c_x = & (* _d2pxy_son_c)[x];
				_d2pxy_son_c_x -> resize(_nbStates);
				for(unsigned int y = 0; y < _nbStates; y++) {
					(* _d2pxy_son_c_x)[y] =  rc * rc * d2Q(x, y);
				}
			}
		}
	}

	computeTreeLikelihood();
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getValue() const
throw (Exception)
{
	//double f = - getLogLikelihood(); // For minimization.
	//if(isnan(f)) f = -log(0.); // (+inf if unlikely!)
	//return f;
	return - getLogLikelihood();
}

/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/	

double DRHomogeneousTreeLikelihood::getDLikelihoodForASiteForARateClass(
	unsigned int site,
	unsigned int rateClass) const
{
	double dl = 0;
	for(unsigned int i = 0; i < _nbStates; i++) {
		dl += _dLikelihoods[_tree -> getRootNode()][_rootPatternLinks[site]][rateClass][i] * _model -> freq(i);
	}
	return dl;
}

/******************************************************************************/	

double DRHomogeneousTreeLikelihood::getDLikelihoodForASite(unsigned int site) const
{
	// Derivative of the sum is the sum of derivatives:
	double dl = 0;
	for(unsigned int i = 0; i < _nbClasses; i++)
		dl += getDLikelihoodForASiteForARateClass(site, i) * _rateDistribution -> getProbability(i);
	return dl;
}

/******************************************************************************/	

double DRHomogeneousTreeLikelihood::getDLogLikelihoodForASite(unsigned int site) const
{
	// d(f(g(x)))/dx = dg(x)/dx . df(g(x))/dg :
	return getDLikelihoodForASite(site) / getLikelihoodForASite(site);
}

/******************************************************************************/	

double DRHomogeneousTreeLikelihood::getDLogLikelihood() const
{
	// Derivative of the sum is the sum of derivatives:
	double dl = 0;
	for(unsigned int i = 0; i < _nbSites; i++)
		dl += getDLogLikelihoodForASite(i);
	return dl;
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getFirstOrderDerivative(const string & variable) const
throw (Exception)
{ 
	Parameter * p = _parameters.getParameter(variable);
	if(p == NULL) throw ParameterNotFoundException("HomogeneousTreeLikelihood::df", variable);
	if(getRateDistributionParameters().getParameter(variable) != NULL) {
		cout << "DEBUB: WARNING!!! Derivatives respective to rate distribution parameter are not implemented." << endl;
		return log(-1.);
	}
	if(getSubstitutionModelParameters().getParameter(variable) != NULL) {
		cout << "DEBUB: WARNING!!! Derivatives respective to substitution model parameters are not implemented." << endl;
		return log(-1.);
	}
	
	const_cast<DRHomogeneousTreeLikelihood *>(this) -> computeTreeDLikelihood(variable);
	return - getDLogLikelihood();
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeTreeDLikelihood(const string & variable)
{
	// Get the node with the branch whose length must be derivated:
	int brI = TextTools::toInt(variable.substr(5));
	Node * branch = _nodes[brI];
	Node * father = branch -> getFather();
	VVVdouble * _dLikelihoods_father = & _dLikelihoods[father];
	map<const Node *, VVVdouble> * _likelihoods_father = & _likelihoods[father];
	
	// Compute dLikelihoods array for the father node.
	// Fist initialize to 1:
	resetLikelihoodArray(* _dLikelihoods_father);

	unsigned int nbNodes = father -> getNumberOfSons();
	for(unsigned int l = 0; l < nbNodes; l++) {
		
		Node * son = father -> getSon(l);

		VVVdouble * _likelihoods_son = & (* _likelihoods_father)[son];

		if(son == branch) {
			VVVdouble * _dpxy_son = & _dpxy[son];
			for(unsigned int i = 0; i < _nbDistinctSites; i++) {
				VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[i];
				VVdouble * _dLikelihoods_father_i = & (* _dLikelihoods_father)[i];
				for(unsigned int c = 0; c < _nbClasses; c++) {
					Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
					Vdouble * _dLikelihoods_father_i_c = & (* _dLikelihoods_father_i)[c];
					VVdouble * _dpxy_son_c = & (* _dpxy_son)[c];
					for(unsigned int x = 0; x < _nbStates; x++) {
						double dl = 0;
						Vdouble * _dpxy_son_c_x = & (* _dpxy_son_c)[x];
						for(unsigned int y = 0; y < _nbStates; y++) {
							dl += (* _dpxy_son_c_x)[y] * (* _likelihoods_son_i_c)[y];
						}
						(* _dLikelihoods_father_i_c)[x] *= dl;
					}
				}
			}
		} else {
			VVVdouble * _pxy_son = & _pxy[son];
			for(unsigned int i = 0; i < _nbDistinctSites; i++) {
				VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[i];
				VVdouble * _dLikelihoods_father_i = & (* _dLikelihoods_father)[i];
				for(unsigned int c = 0; c < _nbClasses; c++) {
					Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
					Vdouble * _dLikelihoods_father_i_c = & (* _dLikelihoods_father_i)[c];
					VVdouble * _pxy_son_c = & (* _pxy_son)[c];
					for(unsigned int x = 0; x < _nbStates; x++) {
						double dl = 0;
						Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
						for(unsigned int y = 0; y < _nbStates; y++) {
							dl += (* _pxy_son_c_x)[y] * (* _likelihoods_son_i_c)[y];
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

void DRHomogeneousTreeLikelihood::computeDownSubtreeDLikelihood(const Node * node)
{
	const Node * father = node -> getFather();
	// We assume that the _dLikelihoods array has been filled for the current node 'node'.
	// We will evaluate the array for the father node.
	if(father == NULL) return; // We reached the root!
		
	// Compute dLikelihoods array for the father node.
	// Fist initialize to 1:
	VVVdouble * _dLikelihoods_father = & _dLikelihoods[father];
	map<const Node *, VVVdouble> * _likelihoods_father = & _likelihoods[father];
	resetLikelihoodArray(* _dLikelihoods_father);

	unsigned int nbNodes = father -> getNumberOfSons();
	for(unsigned int l = 0; l < nbNodes; l++) {
		const Node * son = father -> getSon(l);

		VVVdouble * _pxy_son = & _pxy[son];
	
		if(son == node) {
			VVVdouble * _dLikelihoods_son = & _dLikelihoods[son];
			for(unsigned int i = 0; i < _nbDistinctSites; i++) {
				VVdouble * _dLikelihoods_son_i = & (* _dLikelihoods_son)[i];
				VVdouble * _dLikelihoods_father_i = & (* _dLikelihoods_father)[i];
				for(unsigned int c = 0; c < _nbClasses; c++) {
					Vdouble * _dLikelihoods_son_i_c = & (* _dLikelihoods_son_i)[c];
					Vdouble * _dLikelihoods_father_i_c = & (* _dLikelihoods_father_i)[c];
					VVdouble * _pxy_son_c = & (* _pxy_son)[c];
					for(unsigned int x = 0; x < _nbStates; x++) {
						double dl = 0;
						Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
						for(unsigned int y = 0; y < _nbStates; y++) {
							dl +=(* _pxy_son_c_x)[y] * (* _dLikelihoods_son_i_c)[y];
						}
						(* _dLikelihoods_father_i_c)[x] *= dl;
					}
				}
			}
		} else {
			VVVdouble * _likelihoods_son = & (* _likelihoods_father)[son];
			for(unsigned int i = 0; i < _nbDistinctSites; i++) {
				VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[i];
				VVdouble * _dLikelihoods_father_i = & (* _dLikelihoods_father)[i];
				for(unsigned int c = 0; c < _nbClasses; c++) {
					Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
					Vdouble * _dLikelihoods_father_i_c = & (* _dLikelihoods_father_i)[c];
					VVdouble * _pxy_son_c = & (* _pxy_son)[c];
					for(unsigned int x = 0; x < _nbStates; x++) {
						double dl = 0;
						Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
						for(unsigned int y = 0; y < _nbStates; y++) {
							dl += (* _pxy_son_c_x)[y] * (* _likelihoods_son_i_c)[y];
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

double DRHomogeneousTreeLikelihood::getD2LikelihoodForASiteForARateClass(
	unsigned int site,
	unsigned int rateClass) const
{
	double d2l = 0;
	for(unsigned int i = 0; i < _nbStates; i++) {
		d2l += _d2Likelihoods[_tree -> getRootNode()][_rootPatternLinks[site]][rateClass][i] * _model -> freq(i);
	}
	return d2l;
}

/******************************************************************************/	

double DRHomogeneousTreeLikelihood::getD2LikelihoodForASite(unsigned int site) const
{
	// Derivative of the sum is the sum of derivatives:
	double d2l = 0;
	for(unsigned int i = 0; i < _nbClasses; i++)
		d2l += getD2LikelihoodForASiteForARateClass(site, i) * _rateDistribution -> getProbability(i);
	return d2l;
}

/******************************************************************************/	

double DRHomogeneousTreeLikelihood::getD2LogLikelihoodForASite(unsigned int site) const
{
	return getD2LikelihoodForASite(site) / getLikelihoodForASite(site)
	- pow( getDLikelihoodForASite(site) / getLikelihoodForASite(site), 2);
}

/******************************************************************************/	

double DRHomogeneousTreeLikelihood::getD2LogLikelihood() const
{
	// Derivative of the sum is the sum of derivatives:
	double dl = 0;
	for(unsigned int i = 0; i < _nbSites; i++)
		dl += getD2LogLikelihoodForASite(i);
	return dl;
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getSecondOrderDerivative(const string & variable) const
throw (Exception)
{
	Parameter * p = _parameters.getParameter(variable);
	if(p == NULL) throw ParameterNotFoundException("HomogeneousTreeLikelihood::df", variable);
	if(getRateDistributionParameters().getParameter(variable) != NULL) {
		cout << "DEBUB: WARNING!!! Derivatives respective to rate distribution parameter are not implemented." << endl;
		return log(-1.);
	}
	if(getSubstitutionModelParameters().getParameter(variable) != NULL) {
		cout << "DEBUB: WARNING!!! Derivatives respective to substitution model parameters are not implemented." << endl;
		return log(-1.);
	}
	
	const_cast<DRHomogeneousTreeLikelihood *>(this) -> computeTreeD2Likelihood(variable);
	return - getD2LogLikelihood();
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeTreeD2Likelihood(const string & variable)
{
	// Get the node with the branch whose length must be derivated:
	int brI = TextTools::toInt(variable.substr(5));
	const Node * branch = _nodes[brI];
	const Node * father = branch -> getFather();
	
	// Compute dLikelihoods array for the father node.
	// Fist initialize to 1:
	VVVdouble * _d2Likelihoods_father = & _d2Likelihoods[father]; 
	map<const Node *, VVVdouble> * _likelihoods_father = & _likelihoods[father];
	resetLikelihoodArray(* _d2Likelihoods_father);

	unsigned int nbNodes = father -> getNumberOfSons();
	for(unsigned int l = 0; l < nbNodes; l++) {
		
		const Node * son = father -> getSon(l);
		
		VVVdouble * _likelihoods_son = & (* _likelihoods_father)[son];

		if(son == branch) {
			VVVdouble * _d2pxy_son = & _d2pxy[son];
			for(unsigned int i = 0; i < _nbDistinctSites; i++) {
				VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[i];
				VVdouble * _d2Likelihoods_father_i = & (* _d2Likelihoods_father)[i];
				for(unsigned int c = 0; c < _nbClasses; c++) {
					Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
					Vdouble * _d2Likelihoods_father_i_c = & (* _d2Likelihoods_father_i)[c];
					VVdouble * _d2pxy_son_c = & (* _d2pxy_son)[c];
					for(unsigned int x = 0; x < _nbStates; x++) {
						double d2l = 0;
						Vdouble * _d2pxy_son_c_x = & (* _d2pxy_son_c)[x];
						for(unsigned int y = 0; y < _nbStates; y++) {
							d2l += (* _d2pxy_son_c_x)[y] * (* _likelihoods_son_i_c)[y];
						}
						(* _d2Likelihoods_father_i_c)[x] *= d2l;
					}
				}
			}
		} else {
			VVVdouble * _pxy_son = & _pxy[son];
			for(unsigned int i = 0; i < _nbDistinctSites; i++) {
				VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[i];
				VVdouble * _d2Likelihoods_father_i = & (* _d2Likelihoods_father)[i];
				for(unsigned int c = 0; c < _nbClasses; c++) {
					Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
					Vdouble * _d2Likelihoods_father_i_c = & (* _d2Likelihoods_father_i)[c];
					VVdouble * _pxy_son_c = & (* _pxy_son)[c];
					for(unsigned int x = 0; x < _nbStates; x++) {
						double d2l = 0;
						Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
						for(unsigned int y = 0; y < _nbStates; y++) {
							d2l += (* _pxy_son_c_x)[y] * (* _likelihoods_son_i_c)[y];
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

void DRHomogeneousTreeLikelihood::computeDownSubtreeD2Likelihood(const Node * node)
{
	const Node * father = node -> getFather();
	// We assume that the _dLikelihoods array has been filled for the current node 'node'.
	// We will evaluate the array for the father node.
	if(father == NULL) return; // We reached the root!
		
	// Compute dLikelihoods array for the father node.
	// Fist initialize to 1:
	VVVdouble * _d2Likelihoods_father = & _d2Likelihoods[father];
	map<const Node *, VVVdouble> * _likelihoods_father = & _likelihoods[father];
	resetLikelihoodArray(* _d2Likelihoods_father);

	unsigned int nbNodes = father -> getNumberOfSons();
	for(unsigned int l = 0; l < nbNodes; l++) {
		const Node * son = father -> getSon(l);

		VVVdouble * _pxy_son = & _pxy[son];
	
		if(son == node) {
			VVVdouble * _d2Likelihoods_son = & _d2Likelihoods[son];
			for(unsigned int i = 0; i < _nbDistinctSites; i++) {
				VVdouble * _d2Likelihoods_son_i = & (* _d2Likelihoods_son)[i];
				VVdouble * _d2Likelihoods_father_i = & (* _d2Likelihoods_father)[i];
				for(unsigned int c = 0; c < _nbClasses; c++) {
					Vdouble * _d2Likelihoods_son_i_c = & (* _d2Likelihoods_son_i)[c];
					Vdouble * _d2Likelihoods_father_i_c = & (* _d2Likelihoods_father_i)[c];
					VVdouble * _pxy_son_c = & (* _pxy_son)[c];
					for(unsigned int x = 0; x < _nbStates; x++) {
						double d2l = 0;
						Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
						for(unsigned int y = 0; y < _nbStates; y++) {
							d2l += (* _pxy_son_c_x)[y] * (* _d2Likelihoods_son_i_c)[y];
						}
						(* _d2Likelihoods_father_i_c)[x] *= d2l;
					}
				}
			}
		} else {
			VVVdouble * _likelihoods_son = & (* _likelihoods_father)[son];
			for(unsigned int i = 0; i < _nbDistinctSites; i++) {
				VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[i];
				VVdouble * _d2Likelihoods_father_i = & (* _d2Likelihoods_father)[i];
				for(unsigned int c = 0; c < _nbClasses; c++) {
					Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
					Vdouble * _d2Likelihoods_father_i_c = & (* _d2Likelihoods_father_i)[c];
					VVdouble * _pxy_son_c = & (* _pxy_son)[c];
					for(unsigned int x = 0; x < _nbStates; x++) {
						double dl = 0;
						Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
						for(unsigned int y = 0; y < _nbStates; y++) {
							dl += (* _pxy_son_c_x)[y] * (* _likelihoods_son_i_c)[y];
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

void DRHomogeneousTreeLikelihood::resetLikelihoodArrays(const Node * node)
{
	for(unsigned int n = 0; n < node -> getNumberOfSons(); n++) {
		const Node * subNode = node -> getSon(n);
		resetLikelihoodArray(_likelihoods[node][subNode]);
	}
	if(node -> hasFather()) {
		const Node * father = node -> getFather();
		resetLikelihoodArray(_likelihoods[node][father]);
	}
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::initTreeLikelihoods(const Node * node, const SequenceContainer & sequences) throw (Exception)
{
	if(node -> isLeaf()) {
		// Init leaves likelihoods:
		//cout << "Leaf:\t" << node -> getName() << endl;
		VVdouble * _leavesLikelihoods_leaf = & _leavesLikelihoods[node];
		_leavesLikelihoods_leaf -> resize(_nbDistinctSites);
		for(unsigned int i = 0; i < _nbDistinctSites; i++) {
			Vdouble * _leavesLikelihoods_leaf_i = & (* _leavesLikelihoods_leaf)[i];
			_leavesLikelihoods_leaf_i -> resize(_nbStates);
			for(unsigned int s = 0; s < _nbStates; s++) {
				//Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
				//otherwise value set to 0:
				try {
					int state = sequences.getSequence(node -> getName()) -> getValue(i);
					(* _leavesLikelihoods_leaf_i)[s] = _model -> getInitValue(s, state);
				} catch (SequenceNotFoundException snfe) {
					throw SequenceNotFoundException("DRHomogeneousTreeLikelihood::initTreelikelihoods. Leaf name in tree not found in site conainer: ", (node -> getName()));
				}
			}
		}

		//Initialize likelihood vector:
		map<const Node *, VVVdouble> * _likelihoods_node = & _likelihoods[node];

		//Only one neighbor: the father node!
		const Node * neighbor = node -> getFather();
		VVVdouble * _likelihoods_node_neighbor = & (* _likelihoods_node)[neighbor];
		
		_likelihoods_node_neighbor -> resize(_nbDistinctSites);

		//The neighbor can't be a leaf: 
		for(unsigned int i = 0; i < _nbDistinctSites; i++) {
			VVdouble * _likelihoods_node_neighbor_i = & (* _likelihoods_node_neighbor)[i];
			_likelihoods_node_neighbor_i -> resize(_nbClasses);
			for(unsigned int c = 0; c < _nbClasses; c++) {
				Vdouble * _likelihoods_node_neighbor_i_c = & (* _likelihoods_node_neighbor_i)[c];
				_likelihoods_node_neighbor_i_c -> resize(_nbStates);
				for(unsigned int s = 0; s < _nbStates; s++) {
					(* _likelihoods_node_neighbor_i_c)[s] = 1.; //All likelihoods are initialized to 1.
				}
			}
		}
	
	
	} else {
		//cout << "Node:\t" << node -> getId() << endl;
		// We initialize each son node first:
		unsigned int nbSonNodes = node -> getNumberOfSons();
		for(unsigned int l = 0; l < nbSonNodes; l++) {
			//For each son node,
			initTreeLikelihoods(node -> getSon(l), sequences);
		}

		//Initialize likelihood vector:
		map<const Node *, VVVdouble> * _likelihoods_node = & _likelihoods[node];
	
		int nbSons = node -> getNumberOfSons();
	
		for(int n = (node -> hasFather() ? -1 : 0); n < nbSons; n++) {
			const Node * neighbor = (* node)[n];
			VVVdouble * _likelihoods_node_neighbor = & (* _likelihoods_node)[neighbor];
		
			_likelihoods_node_neighbor -> resize(_nbDistinctSites);

			if(neighbor -> isLeaf()) {
				VVdouble * _leavesLikelihoods_leaf = & _leavesLikelihoods[neighbor];
				for(unsigned int i = 0; i < _nbDistinctSites; i++) {
					Vdouble  * _leavesLikelihoods_leaf_i = & (* _leavesLikelihoods_leaf)[i];
					VVdouble * _likelihoods_node_neighbor_i = & (* _likelihoods_node_neighbor)[i];
					_likelihoods_node_neighbor_i -> resize(_nbClasses);
					for(unsigned int c = 0; c < _nbClasses; c++) {
						Vdouble * _likelihoods_node_neighbor_i_c = & (* _likelihoods_node_neighbor_i)[c];
						_likelihoods_node_neighbor_i_c -> resize(_nbStates);
						for(unsigned int s = 0; s < _nbStates; s++) {
							(* _likelihoods_node_neighbor_i_c)[s] = (* _leavesLikelihoods_leaf_i)[s];
						}
					}
				}
			} else {
				for(unsigned int i = 0; i < _nbDistinctSites; i++) {
					VVdouble * _likelihoods_node_neighbor_i = & (* _likelihoods_node_neighbor)[i];
					_likelihoods_node_neighbor_i -> resize(_nbClasses);
					for(unsigned int c = 0; c < _nbClasses; c++) {
						Vdouble * _likelihoods_node_neighbor_i_c = & (* _likelihoods_node_neighbor_i)[c];
						_likelihoods_node_neighbor_i_c -> resize(_nbStates);
						for(unsigned int s = 0; s < _nbStates; s++) {
							(* _likelihoods_node_neighbor_i_c)[s] = 1.; //All likelihoods are initialized to 1.
						}
					}
				}
			}
		}	
	}
	
	// Initialize d and d2 likelihoods:
	VVVdouble * _dLikelihoods_node = & _dLikelihoods[node];
	VVVdouble * _d2Likelihoods_node = & _d2Likelihoods[node];
	_dLikelihoods_node -> resize(_nbDistinctSites);
	_d2Likelihoods_node -> resize(_nbDistinctSites);
	for(unsigned int i = 0; i < _nbDistinctSites; i++) {
		VVdouble * _dLikelihoods_node_i = & (* _dLikelihoods_node)[i];
		VVdouble * _d2Likelihoods_node_i = & (* _d2Likelihoods_node)[i];
		_dLikelihoods_node_i -> resize(_nbClasses);
		_d2Likelihoods_node_i -> resize(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			Vdouble * _dLikelihoods_node_i_c = & (* _dLikelihoods_node_i)[c];
			Vdouble * _d2Likelihoods_node_i_c = & (* _d2Likelihoods_node_i)[c];
			_dLikelihoods_node_i_c -> resize(_nbStates);
			_d2Likelihoods_node_i_c -> resize(_nbStates);
			for(unsigned int s = 0; s < _nbStates; s++) {
				(* _dLikelihoods_node_i_c)[s] = 0.;
				(* _d2Likelihoods_node_i_c)[s] = 0.;
			}
		}
	}
}

/******************************************************************************/
	
void DRHomogeneousTreeLikelihood::computeTreeLikelihood()
{
	computeSubtreeLikelihoodPostfix(_tree -> getRootNode());
	computeSubtreeLikelihoodPrefix(_tree -> getRootNode());
	computeRootLikelihood();
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeSubtreeLikelihoodPostfix(const Node * node)
{
	if(node -> isLeaf()) return;

	// Set all likelihood arrays to 1 for a start:
	resetLikelihoodArrays(node);
	
	map<const Node *, VVVdouble> * _likelihoods_node = & _likelihoods[node];
	unsigned int nbNodes = node -> getNumberOfSons();
	for(unsigned int l = 0; l < nbNodes; l++) {
		//For each son node...	

		const Node * son = node -> getSon(l);
		VVVdouble * _likelihoods_node_son = & (* _likelihoods_node)[son];
		
		if(son -> isLeaf()) {
			VVdouble * _likelihoods_leaf = & _leavesLikelihoods[son];
			for(unsigned int i = 0; i < _nbDistinctSites; i++) {
				//For each site in the sequence,
				Vdouble * _likelihoods_leaf_i     = & (* _likelihoods_leaf)[i];
				VVdouble * _likelihoods_node_son_i = & (* _likelihoods_node_son)[i];
				for(unsigned int c = 0; c < _nbClasses; c++) {
					//For each rate classe,
					Vdouble * _likelihoods_node_son_i_c = & (* _likelihoods_node_son_i)[c];
					for(unsigned int x = 0; x < _nbStates; x++) {
						//For each initial state,
						(* _likelihoods_node_son_i_c)[x] = (* _likelihoods_leaf_i)[x];
					}
				}
			}
		} else {
			
			computeSubtreeLikelihoodPostfix(son); //Recursive method:
			unsigned int nbSons = son -> getNumberOfSons();
			map<const Node *, VVVdouble> * _likelihoods_son = & _likelihoods[son];

			for(unsigned int n = 0; n < nbSons; n++) {

				const Node * sonson = son -> getSon(n);
				VVVdouble * _pxy_sonson = & _pxy[sonson];
				VVVdouble * _likelihoods_son_son = & (* _likelihoods_son)[sonson];

				for(unsigned int i = 0; i < _nbDistinctSites; i++) {
					//For each site in the sequence,
					VVdouble * _likelihoods_son_son_i  = & (* _likelihoods_son_son)[i];
					VVdouble * _likelihoods_node_son_i = & (* _likelihoods_node_son)[i];
					for(unsigned int c = 0; c < _nbClasses; c++) {
						//For each rate classe,
						Vdouble * _likelihoods_son_son_i_c  = & (* _likelihoods_son_son_i)[c];
						Vdouble * _likelihoods_node_son_i_c = & (* _likelihoods_node_son_i)[c];
						VVdouble * _pxy_sonson_c = & (* _pxy_sonson)[c];
						for(unsigned int x = 0; x < _nbStates; x++) {
							//For each initial state,
							Vdouble * _pxy_sonson_c_x = & (* _pxy_sonson_c)[x];
							double likelihood = 0;
							for(unsigned int y = 0; y < _nbStates; y++) {
								likelihood += (* _pxy_sonson_c_x)[y] * (* _likelihoods_son_son_i_c)[y];
							}
							// We store this conditionnal likelihood into the corresponding array:
							(* _likelihoods_node_son_i_c)[x] *= likelihood;
						}
					}
				}
			}
		}
	}
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeSubtreeLikelihoodPrefix(const Node * node)
{
	//if(node -> isLeaf()) {
	//	return;
	//}
	if(! node -> hasFather()) { // 'node' is the root of the tree.
		// Just call the method on each son node:
		unsigned int nbSons = node -> getNumberOfSons();
		for(unsigned int n = 0; n < nbSons; n++) computeSubtreeLikelihoodPrefix(node -> getSon(n));
		return;
	} else {
		const Node * father = node -> getFather();
		vector<const Node *> nodes;
	  // Add brothers:
		unsigned int nbFatherSons = father -> getNumberOfSons();
		for(unsigned int n = 0; n < nbFatherSons; n++) {
			const Node * son = father -> getSon(n);
			if(son -> getId() != node -> getId()) nodes.push_back(son); //This is a real brother, not current node!
		}
		// Now the real stuff... We've got to compute the likelihoods for the
		// subtree defined by node 'father'.
		// This is the same as postfix method, but with different subnodes.
		map<const Node *, VVVdouble> * _likelihoods_node = & _likelihoods[node];
		map<const Node *, VVVdouble> * _likelihoods_father = & _likelihoods[father];
		VVVdouble * _likelihoods_node_father = & (* _likelihoods_node)[father];
		
		unsigned int nbSons = nodes.size(); // In case of a bifurcating tree, this is equal to 1, excepted for the root.
		for(unsigned int n = 0; n < nbSons; n++) {

			const Node * fatherSon = nodes[n];
			VVVdouble * _pxy_fatherSon = & _pxy[fatherSon];
			VVVdouble * _likelihoods_father_son = & (* _likelihoods_father)[fatherSon];

			for(unsigned int i = 0; i < _nbDistinctSites; i++) {
				//For each site in the sequence,
				VVdouble * _likelihoods_father_son_i  = & (* _likelihoods_father_son)[i];
				VVdouble * _likelihoods_node_father_i = & (* _likelihoods_node_father)[i];
				for(unsigned int c = 0; c < _nbClasses; c++) {
					//For each rate classe,
					Vdouble * _likelihoods_father_son_i_c  = & (* _likelihoods_father_son_i)[c];
					Vdouble * _likelihoods_node_father_i_c = & (* _likelihoods_node_father_i)[c];
					VVdouble * _pxy_fatherSon_c = & (* _pxy_fatherSon)[c];
					for(unsigned int x = 0; x < _nbStates; x++) {
						//For each initial state,
						Vdouble * _pxy_fatherSon_c_x = & (* _pxy_fatherSon_c)[x];
						double likelihood = 0;
						for(unsigned int y = 0; y < _nbStates; y++) {
							likelihood += (* _pxy_fatherSon_c_x)[y] * (* _likelihoods_father_son_i_c)[y];
						}
						// We store this conditionnal likelihood into the corresponding array:
						(* _likelihoods_node_father_i_c)[x] *= likelihood;
					}
				}
			}
		}
		
		if(father -> hasFather()) {
			// Also take grand-father into account:
			const Node * fatherFather = father -> getFather();
			VVVdouble * _pxy_fatherFather = & _pxy[father]; //!!! the difference here is that we use
			                                                //!!! _pxy[father] instead of _pxy[fatherFather].
			VVVdouble * _likelihoods_father_father = & (* _likelihoods_father)[fatherFather];

			for(unsigned int i = 0; i < _nbDistinctSites; i++) {
				//For each site in the sequence,
				VVdouble * _likelihoods_father_father_i = & (* _likelihoods_father_father)[i];
				VVdouble * _likelihoods_node_father_i   = & (* _likelihoods_node_father)[i];
				for(unsigned int c = 0; c < _nbClasses; c++) {
					//For each rate classe,
					Vdouble * _likelihoods_father_father_i_c = & (* _likelihoods_father_father_i)[c];
					Vdouble * _likelihoods_node_father_i_c   = & (* _likelihoods_node_father_i)[c];
					VVdouble * _pxy_fatherFather_c = & (* _pxy_fatherFather)[c];
					for(unsigned int x = 0; x < _nbStates; x++) {
						//For each initial state,
						Vdouble * _pxy_fatherFather_c_x = & (* _pxy_fatherFather_c)[x];
						double likelihood = 0;
						for(unsigned int y = 0; y < _nbStates; y++) {
							likelihood += (* _pxy_fatherFather_c_x)[y] * (* _likelihoods_father_father_i_c)[y];
						}
						// We store this conditionnal likelihood into the corresponding array:
						(* _likelihoods_node_father_i_c)[x] *= likelihood;
					}
				}
			}
			
		}


		// Call the method on each son node:
		unsigned int nbNodeSons = node -> getNumberOfSons();
		for(unsigned int i = 0; i < nbNodeSons; i++)
			computeSubtreeLikelihoodPrefix(node -> getSon(i)); //Recursive method.
	}
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeRootLikelihood()
{
	const Node * root = _tree -> getRootNode();
	// Set all likelihoods to 1 for a start:
	resetLikelihoodArray(_rootLikelihoods);

	map<const Node *, VVVdouble> * _likelihoods_root = & _likelihoods[root];
	unsigned int nbNodes = root -> getNumberOfSons();
	for(unsigned int n = 0; n < nbNodes; n++) {

		const Node * son = root -> getSon(n);
		VVVdouble * _pxy_son = & _pxy[son];
		VVVdouble * _likelihoods_root_son = & (* _likelihoods_root)[son];

		for(unsigned int i = 0; i < _nbDistinctSites; i++) {
			//For each site in the sequence,
			VVdouble * _likelihoods_root_son_i  = & (* _likelihoods_root_son)[i];
			VVdouble * _rootLikelihoods_i = & _rootLikelihoods[i];
			for(unsigned int c = 0; c < _nbClasses; c++) {
				//For each rate classe,
				Vdouble * _likelihoods_root_son_i_c  = & (* _likelihoods_root_son_i)[c];
				Vdouble * _rootLikelihoods_i_c = & (* _rootLikelihoods_i)[c];
				VVdouble * _pxy_son_c = & (* _pxy_son)[c];
				for(unsigned int x = 0; x < _nbStates; x++) {
					//For each initial state,
					Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
					double likelihood = 0;
					for(unsigned int y = 0; y < _nbStates; y++) {
						likelihood += (* _pxy_son_c_x)[y] * (* _likelihoods_root_son_i_c)[y];
					}
					// We store this conditionnal likelihood into the corresponding array:
					(* _rootLikelihoods_i_c)[x] *= likelihood;
				}
			}
		}
	}
	//displayLikelihood(root);
	//cout << "Result : " << endl;
	//displayLikelihood(_rootLikelihoods);
}

/******************************************************************************/

VVVdouble DRHomogeneousTreeLikelihood::computeLikelihoodAtNode(const Node * node)
{
	VVVdouble likelihoodArray(_nbDistinctSites);
	map<const Node *, VVVdouble> * _likelihoods_node = & _likelihoods[node];
	
	if(node -> isLeaf()) {
		VVdouble * _leavesLikelihoods_node = & _leavesLikelihoods[node];
		for(unsigned int i = 0; i < _nbDistinctSites; i++) {
			VVdouble * likelihoodArray_i = & likelihoodArray[i];
			Vdouble * _leavesLikelihoods_node_i = & (* _leavesLikelihoods_node)[i];
			likelihoodArray_i -> resize(_nbClasses);
			for(unsigned int c = 0; c < _nbClasses; c++) {
				Vdouble * likelihoodArray_i_c = & (* likelihoodArray_i)[c];
				likelihoodArray_i_c -> resize(_nbStates);
				for(unsigned int x = 0; x < _nbStates; x++) {
					(* likelihoodArray_i_c)[x] = (* _leavesLikelihoods_node_i)[x];
				}
			}
		}
		
	} else {
	
		// Otherwise:
		// Set all likelihoods to 1 for a start:
		for(unsigned int i = 0; i < _nbDistinctSites; i++) {
			VVdouble * likelihoodArray_i = & likelihoodArray[i];
			likelihoodArray_i -> resize(_nbClasses);
			for(unsigned int c = 0; c < _nbClasses; c++) {
				Vdouble * likelihoodArray_i_c = & (* likelihoodArray_i)[c];
				likelihoodArray_i_c -> resize(_nbStates);
				for(unsigned int s = 0; s < _nbStates; s++) {
					(* likelihoodArray_i_c)[s] = 1.;
				}
			}
		}

		unsigned int nbNodes = node -> getNumberOfSons();
		for(unsigned int n = 0; n < nbNodes; n++) {

			const Node * son = node -> getSon(n);
			VVVdouble * _pxy_son = & _pxy[son];
			VVVdouble * _likelihoods_node_son = & (* _likelihoods_node)[son];

			for(unsigned int i = 0; i < _nbDistinctSites; i++) {
				//For each site in the sequence,
				VVdouble * _likelihoods_node_son_i  = & (* _likelihoods_node_son)[i];
				VVdouble * likelihoodArray_i = & likelihoodArray[i];
				for(unsigned int c = 0; c < _nbClasses; c++) {
					//For each rate classe,
					Vdouble * _likelihoods_node_son_i_c  = & (* _likelihoods_node_son_i)[c];
					Vdouble * likelihoodArray_i_c = & (* likelihoodArray_i)[c];
					VVdouble * _pxy_son_c = & (* _pxy_son)[c];
					for(unsigned int x = 0; x < _nbStates; x++) {
						//For each initial state,
						Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
						double likelihood = 0;
						for(unsigned int y = 0; y < _nbStates; y++) {
							likelihood += (* _pxy_son_c_x)[y] * (* _likelihoods_node_son_i_c)[y];
						}
						// We store this conditionnal likelihood into the corresponding array:
						(* likelihoodArray_i_c)[x] *= likelihood;
					}
				}
			}
		}
	}
	
	if(node -> hasFather()) {
		const Node * son = node -> getFather();
		VVVdouble * _pxy_son = & _pxy[node]; // and not son!!!
		VVVdouble * _likelihoods_node_son = & (* _likelihoods_node)[son];

		for(unsigned int i = 0; i < _nbDistinctSites; i++) {
			//For each site in the sequence,
			VVdouble * _likelihoods_node_son_i  = & (* _likelihoods_node_son)[i];
			VVdouble * likelihoodArray_i = & likelihoodArray[i];
			for(unsigned int c = 0; c < _nbClasses; c++) {
				//For each rate classe,
				Vdouble * _likelihoods_node_son_i_c  = & (* _likelihoods_node_son_i)[c];
				Vdouble * likelihoodArray_i_c = & (* likelihoodArray_i)[c];
				VVdouble * _pxy_son_c = & (* _pxy_son)[c];
				for(unsigned int x = 0; x < _nbStates; x++) {
					//For each initial state,
					Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
					double likelihood = 0;
					for(unsigned int y = 0; y < _nbStates; y++) {
						likelihood += (* _pxy_son_c_x)[y] * (* _likelihoods_node_son_i_c)[y];
					}
					// We store this conditionnal likelihood into the corresponding array:
					(* likelihoodArray_i_c)[x] *= likelihood;
				}
			}
		}
	}
	return likelihoodArray;
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::displayLikelihood(const Node * node)
{
	cout << "Likelihoods at node " << node -> getId() << ": " << endl;
	for(unsigned int n = 0; n < node -> getNumberOfSons(); n++) {
		const Node * subNode = node -> getSon(n);
		cout << "Array for sub-node " << subNode -> getId() << endl;
		displayLikelihoodArray(_likelihoods[node][subNode]);
	}
	if(node -> hasFather()) {
		const Node * father = node -> getFather();
		cout << "Array for father node " << father -> getId() << endl;
		displayLikelihoodArray(_likelihoods[node][father]);
	}
	cout << "                                         ***" << endl;
}

/*******************************************************************************/

