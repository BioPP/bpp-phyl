//
// File: HomogeneousTreeLikelihood.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Fri Oct 17 18:14:51 2003
//

#include "HomogeneousTreeLikelihood.h"
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

HomogeneousTreeLikelihood::HomogeneousTreeLikelihood(
	Tree & tree,
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
		throw AlphabetMismatchException("HomogeneousTreeLikelihood::HomogeneousTreeLikelihood. Data and model must have the same alphabet type.",
				_data -> getAlphabet(),
				_model -> getAlphabet());

	_rateDistribution = rDist;
	
	_nodes = _tree -> getNodes();
	
	_nodes.pop_back(); //Remove the root node (the last added!).
	
	_nbSites   = _data -> getNumberOfSites();
	_nbClasses = _rateDistribution -> getNumberOfCategories();
	_nbStates  = _model -> getAlphabet() -> getSize();
	_nbNodes   = _nodes.size();
	
	//Init _likelihoods:
	if(verbose) ApplicationTools::message << "Homogeneous Tree Likelihood" << endl;	
	if(verbose) ApplicationTools::displayTask("Init likelihoods arrays recursively");
	const SiteContainer * subSubSequencesShrunk =
		initTreeLikelihoodsWithPatterns(_tree -> getRootNode(), *_data);
	if(verbose) ApplicationTools::displayTaskDone();
	
	if(verbose) ApplicationTools::displayResult("Number of distinct sites",
			TextTools::toString(subSubSequencesShrunk -> getNumberOfSites()));
	
	//Initialize root patterns:
	if(verbose) ApplicationTools::displayTask("Init root patterns");
	_rootPatternLinks.resize(_nbSites);
	for(unsigned int i = 0; i < _nbSites; i++) {
		const Site * site1 =  _data -> getSite(i);
		for(unsigned int ii = 0; ii < subSubSequencesShrunk -> getNumberOfSites(); ii++) {
			if(SiteTools::areSitesIdentical(* subSubSequencesShrunk -> getSite(ii), * site1)) {
				//_rootPatternLinks[i] = & _likelihoods[_tree.getRootNode()][ii];
				_rootPatternLinks[i] = ii;
				break;
			}
		}
		//if(_rootPatternLinks[i] == NULL) cerr << "ERROR while initializing HomogeneousTreeLikelihood at site " << i << endl;
	}
	delete subSubSequencesShrunk;
	if(verbose) ApplicationTools::displayTaskDone();
	
	// Now initializes all parameters:
	initParameters();
	fireParameterChanged(_parameters);
}

/******************************************************************************/

HomogeneousTreeLikelihood::~HomogeneousTreeLikelihood() { delete _data; }

/******************************************************************************/

inline double HomogeneousTreeLikelihood::getLikelihood() const
{
	double l = 1.;
	for(unsigned int i = 0; i < _nbSites; i++) {
		l *= getLikelihoodForASite(i);
	}
	return l;
}

/******************************************************************************/

inline double HomogeneousTreeLikelihood::getLogLikelihood() const
{
	double ll = 0;
	for(unsigned int i = 0; i < _nbSites; i++) {
		ll += getLogLikelihoodForASite(i);
	}
	return ll;
}

/******************************************************************************/

inline double HomogeneousTreeLikelihood::getLikelihoodForASite(unsigned int site) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbClasses; i++) {
		l += getLikelihoodForASiteForARateClass(site, i) * _rateDistribution -> getProbability(i);
	}
	return l;
}

/******************************************************************************/

inline double HomogeneousTreeLikelihood::getLogLikelihoodForASite(unsigned int site) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbClasses; i++) {
		l += getLikelihoodForASiteForARateClass(site, i) * _rateDistribution -> getProbability(i);
	}
	//if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
	return log(l);
}

/******************************************************************************/

inline double HomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbStates; i++) {
		//l += _rootPatternLinks[site] -> operator[](rateClass)[i] * _model -> freq(i);
		l += _likelihoods[_tree -> getRootNode()][_rootPatternLinks[site]][rateClass][i] * _model -> freq(i);
	}
	return l;
}

/******************************************************************************/

inline double HomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbStates; i++) {
		//l += _rootPatternLinks[site] -> operator[](rateClass)[i] * _model -> freq(i);
		l += _likelihoods[_tree -> getRootNode()][_rootPatternLinks[site]][rateClass][i] * _model -> freq(i);
	}
	//if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
	return log(l);
}

/******************************************************************************/	

VVdouble HomogeneousTreeLikelihood::getLikelihoodForEachSiteForEachRate() const
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

VVdouble HomogeneousTreeLikelihood::getLogLikelihoodForEachSiteForEachRate() const
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

VVdouble HomogeneousTreeLikelihood::getPosteriorProbabilitiesOfEachRate() const
{
	VVdouble pb = getLikelihoodForEachSiteForEachRate();
	Vdouble  l  = getLikelihoodForEachSite();
	for(unsigned int i = 0; i < _nbSites; i++) {
		for(unsigned int j = 0; j < _nbClasses; j++) pb[i][j] = pb[i][j] * _rateDistribution -> getProbability(j) / l[i]; 
	}
	return pb;
}
	
/******************************************************************************/

Vint HomogeneousTreeLikelihood::getRateClassWithMaxPostProbOfEachSite() const
{
	VVdouble l = getLikelihoodForEachSiteForEachRate();
	Vint classes(_nbSites);
	for(unsigned int i = 0; i < _nbSites; i++) classes[i] = posmax<double>(l[i]);
	return classes;
}

/******************************************************************************/

Vdouble HomogeneousTreeLikelihood::getRateWithMaxPostProbOfEachSite() const
{
	VVdouble l = getLikelihoodForEachSiteForEachRate();
	Vdouble rates(_nbSites);
	for(unsigned int i = 0; i < _nbSites; i++) {
		rates[i] = _rateDistribution -> getCategory(posmax<double>(l[i]));
	}
	return rates;
}

/******************************************************************************/

Vdouble HomogeneousTreeLikelihood::getPosteriorRateOfEachSite() const
{
	VVdouble lr = getLikelihoodForEachSiteForEachRate();
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

void HomogeneousTreeLikelihood::setParameters(const ParameterList & parameters)
	throw (ParameterNotFoundException, ConstraintException)
{
	setParametersValues(parameters);
}

/******************************************************************************/

void HomogeneousTreeLikelihood::fireParameterChanged(const ParameterList & params)
{
	applyParameters();

	// For now we ignore the parameter that changed and we recompute all arrays...
	for(unsigned int l = 0; l < _nbNodes; l++) {
		//For each son node,
		Node * son = _nodes[l];
		double l = son -> getDistanceToFather(); 
	
		//Computes all pxy once for all:
		VVVdouble * _pxy_son = & _pxy[son];
		(* _pxy_son) = VVVdouble(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			(* _pxy_son)[c] = VVdouble(_nbStates);
			Matrix Q = _model -> getPij_t(l * _rateDistribution -> getCategory(c));
			for(unsigned int x = 0; x < _nbStates; x++) {
				(* _pxy_son)[c][x] = Vdouble(_nbStates);
				for(unsigned int y = 0; y < _nbStates; y++) {
					(* _pxy_son)[c][x][y] = Q(x, y);
				}
			}
		}

		//Computes all dpxy/dt once for all:
		VVVdouble * _dpxy_son = & _dpxy[son];
		(* _dpxy_son) = VVVdouble(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			(* _dpxy_son)[c] = VVdouble(_nbStates);
			double rc = _rateDistribution -> getCategory(c);
			Matrix dQ = _model -> getdPij_dt(l * rc);  
			for(unsigned int x = 0; x < _nbStates; x++) {
				(* _dpxy_son)[c][x] = Vdouble(_nbStates);
				for(unsigned int y = 0; y < _nbStates; y++) {
					(* _dpxy_son)[c][x][y] =  rc * dQ(x, y); 
				}
			}
		}
		
		//Computes all d2pxy/dt2 once for all:
		VVVdouble * _d2pxy_son = & _d2pxy[son];
		(* _d2pxy_son) = VVVdouble(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			(* _d2pxy_son)[c] = VVdouble(_nbStates);
			double rc =  _rateDistribution -> getCategory(c);
			Matrix d2Q = _model -> getd2Pij_dt2(l * rc);
			for(unsigned int x = 0; x < _nbStates; x++) {
				(* _d2pxy_son)[c][x] = Vdouble(_nbStates);
				for(unsigned int y = 0; y < _nbStates; y++) {
					(* _d2pxy_son)[c][x][y] =  rc * rc * d2Q(x, y);
				}
			}
		}
	}

	computeTreeLikelihood();
}

/******************************************************************************/

double HomogeneousTreeLikelihood::getValue() const
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

inline double HomogeneousTreeLikelihood::getDLikelihoodForASiteForARateClass(
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

inline double HomogeneousTreeLikelihood::getDLikelihoodForASite(unsigned int site) const
{
	// Derivative of the sum is the sum of derivatives:
	double dl = 0;
	for(unsigned int i = 0; i < _nbClasses; i++)
		dl += getDLikelihoodForASiteForARateClass(site, i) * _rateDistribution -> getProbability(i);
	return dl;
}

/******************************************************************************/	

inline double HomogeneousTreeLikelihood::getDLogLikelihoodForASite(unsigned int site) const
{
	// d(f(g(x)))/dx = dg(x)/dx . df(g(x))/dg :
	return getDLikelihoodForASite(site) / getLikelihoodForASite(site);
}

/******************************************************************************/	

inline double HomogeneousTreeLikelihood::getDLogLikelihood() const
{
	// Derivative of the sum is the sum of derivatives:
	double dl = 0;
	for(unsigned int i = 0; i < _nbSites; i++)
		dl += getDLogLikelihoodForASite(i);
	return dl;
}

/******************************************************************************/

double HomogeneousTreeLikelihood::getFirstOrderDerivative(const string & variable) const
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
	
	const_cast<HomogeneousTreeLikelihood *>(this) -> computeTreeDLikelihood(variable);
	return - getDLogLikelihood();
}

/******************************************************************************/

void HomogeneousTreeLikelihood::computeTreeDLikelihood(const string & variable)
{
	
	// Get the node with the branch whose length must be derivated:
	int brI = TextTools::toInt(variable.substr(5));
	Node * branch = _nodes[brI];
	Node * father = branch -> getFather();
	VVVdouble * _dLikelihoods_father = & _dLikelihoods[father];
	
	// Compute dLikelihoods array for the father node.
	// Fist initialize to 1:
	unsigned int nbSites = _likelihoods[father].size();
	for(unsigned int i = 0; i < nbSites; i++) {
		VVdouble * _dLikelihoods_father_i = & (* _dLikelihoods_father)[i];
		for(unsigned int c = 0; c < _nbClasses; c++) {
			Vdouble * _dLikelihoods_father_i_c = & (* _dLikelihoods_father_i)[c];
			for(unsigned int s = 0; s < _nbStates; s++) {
				(* _dLikelihoods_father_i_c)[s] = 1.;	
			}
		}
	}

	unsigned int nbNodes = father -> getNumberOfSons();
	for(unsigned int l = 0; l < nbNodes; l++) {
		
		Node * son = father -> getSon(l);

		vector <unsigned int> * _patternLinks_father_son = & _patternLinks[father][son];
		VVVdouble * _likelihoods_son = & _likelihoods[son];

		if(son == branch) {
			VVVdouble * _dpxy_son = & _dpxy[son];
			for(unsigned int i = 0; i < nbSites; i++) {
				VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_father_son)[i]];
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
			for(unsigned int i = 0; i < nbSites; i++) {
				VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_father_son)[i]];
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

void HomogeneousTreeLikelihood::computeDownSubtreeDLikelihood(Node * node)
{
	Node * father = node -> getFather();
	// We assume that the _dLikelihoods array has been filled for the current node 'node'.
	// We will evaluate the array for the father node.
	if(father == NULL) return; // We reached the root!
		
	// Compute dLikelihoods array for the father node.
	// Fist initialize to 1:
	unsigned int nbSites = _likelihoods[father].size();
	VVVdouble * _dLikelihoods_father = & _dLikelihoods[father];
	for(unsigned int i = 0; i < nbSites; i++) {
		VVdouble * _dLikelihoods_father_i = & (* _dLikelihoods_father)[i];
		for(unsigned int c = 0; c < _nbClasses; c++) {
			Vdouble * _dLikelihoods_father_i_c = & (* _dLikelihoods_father_i)[c];
			for(unsigned int s = 0; s < _nbStates; s++) {
				(* _dLikelihoods_father_i_c)[s] = 1.;	
			}
		}
	}

	unsigned int nbNodes = father -> getNumberOfSons();
	for(unsigned int l = 0; l < nbNodes; l++) {
		Node * son = father -> getSon(l);

		VVVdouble * _pxy_son = & _pxy[son];
		vector <unsigned int> * _patternLinks_father_son = & _patternLinks[father][son];
	
		if(son == node) {
			VVVdouble * _dLikelihoods_son = & _dLikelihoods[son];
			for(unsigned int i = 0; i < nbSites; i++) {
				VVdouble * _dLikelihoods_son_i = & (* _dLikelihoods_son)[(* _patternLinks_father_son)[i]];
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
			VVVdouble * _likelihoods_son = & _likelihoods[son];
			for(unsigned int i = 0; i < nbSites; i++) {
				VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_father_son)[i]];
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

inline double HomogeneousTreeLikelihood::getD2LikelihoodForASiteForARateClass(
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

inline double HomogeneousTreeLikelihood::getD2LikelihoodForASite(unsigned int site) const
{
	// Derivative of the sum is the sum of derivatives:
	double d2l = 0;
	for(unsigned int i = 0; i < _nbClasses; i++)
		d2l += getD2LikelihoodForASiteForARateClass(site, i) * _rateDistribution -> getProbability(i);
	return d2l;
}

/******************************************************************************/	

inline double HomogeneousTreeLikelihood::getD2LogLikelihoodForASite(unsigned int site) const
{
	return getD2LikelihoodForASite(site) / getLikelihoodForASite(site)
	- pow( getDLikelihoodForASite(site) / getLikelihoodForASite(site), 2);
}

/******************************************************************************/	

inline double HomogeneousTreeLikelihood::getD2LogLikelihood() const
{
	// Derivative of the sum is the sum of derivatives:
	double dl = 0;
	for(unsigned int i = 0; i < _nbSites; i++)
		dl += getD2LogLikelihoodForASite(i);
	return dl;
}

/******************************************************************************/

double HomogeneousTreeLikelihood::getSecondOrderDerivative(const string & variable) const
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
	
	const_cast<HomogeneousTreeLikelihood *>(this) -> computeTreeD2Likelihood(variable);
	return - getD2LogLikelihood();
}

/******************************************************************************/

void HomogeneousTreeLikelihood::computeTreeD2Likelihood(const string & variable)
{
	
	// Get the node with the branch whose length must be derivated:
	int brI = TextTools::toInt(variable.substr(5));
	Node * branch = _nodes[brI];
	Node * father = branch -> getFather();
	
	// Compute dLikelihoods array for the father node.
	// Fist initialize to 1:
	unsigned int nbSites = _likelihoods[father].size();
	VVVdouble * _d2Likelihoods_father = & _d2Likelihoods[father]; 
	for(unsigned int i = 0; i < nbSites; i++) {
		VVdouble * _d2Likelihoods_father_i = & (* _d2Likelihoods_father)[i];
		for(unsigned int c = 0; c < _nbClasses; c++) {
			Vdouble * _d2Likelihoods_father_i_c = & (* _d2Likelihoods_father_i)[c];
			for(unsigned int s = 0; s < _nbStates; s++) {
				(* _d2Likelihoods_father_i_c)[s] = 1.;	
			}
		}
	}

	unsigned int nbNodes = father -> getNumberOfSons();
	for(unsigned int l = 0; l < nbNodes; l++) {
		
		Node * son = father -> getSon(l);
		
		vector <unsigned int> * _patternLinks_father_son = & _patternLinks[father][son];
		VVVdouble * _likelihoods_son = & _likelihoods[son];

		if(son == branch) {
			VVVdouble * _d2pxy_son = & _d2pxy[son];
			for(unsigned int i = 0; i < nbSites; i++) {
				VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_father_son)[i]];
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
			for(unsigned int i = 0; i < nbSites; i++) {
				VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_father_son)[i]];
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

void HomogeneousTreeLikelihood::computeDownSubtreeD2Likelihood(Node * node)
{
	Node * father = node -> getFather();
	// We assume that the _dLikelihoods array has been filled for the current node 'node'.
	// We will evaluate the array for the father node.
	if(father == NULL) return; // We reached the root!
		
	// Compute dLikelihoods array for the father node.
	// Fist initialize to 1:
	unsigned int nbSites = _likelihoods[father].size();
	VVVdouble * _d2Likelihoods_father = & _d2Likelihoods[father];
	for(unsigned int i = 0; i < nbSites; i++) {
		VVdouble * _d2Likelihoods_father_i = & (* _d2Likelihoods_father)[i];
		for(unsigned int c = 0; c < _nbClasses; c++) {
			Vdouble * _d2Likelihoods_father_i_c = & (* _d2Likelihoods_father_i)[c];
			for(unsigned int s = 0; s < _nbStates; s++) {
				(* _d2Likelihoods_father_i_c)[s] = 1.;	
			}
		}
	}

	unsigned int nbNodes = father -> getNumberOfSons();
	for(unsigned int l = 0; l < nbNodes; l++) {
		Node * son = father -> getSon(l);

		VVVdouble * _pxy_son = & _pxy[son];
		vector <unsigned int> * _patternLinks_father_son = & _patternLinks[father][son];
	
		if(son == node) {
			VVVdouble * _d2Likelihoods_son = & _d2Likelihoods[son];
			for(unsigned int i = 0; i < nbSites; i++) {
				VVdouble * _d2Likelihoods_son_i = & (* _d2Likelihoods_son)[(* _patternLinks_father_son)[i]];
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
			VVVdouble * _likelihoods_son = & _likelihoods[son];
			for(unsigned int i = 0; i < nbSites; i++) {
				VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_father_son)[i]];
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

void HomogeneousTreeLikelihood::applyParameters() throw (Exception)
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

ParameterList HomogeneousTreeLikelihood::getBranchLengthsParameters() const {
	return _brLenParameters.getCommonParametersWith(_parameters);
}

/******************************************************************************/

ParameterList HomogeneousTreeLikelihood::getSubstitutionModelParameters() const {
	return _model -> getParameters().getCommonParametersWith(_parameters);
}

/******************************************************************************/

ParameterList HomogeneousTreeLikelihood::getRateDistributionParameters() const {
	return _rateDistribution -> getParameters().getCommonParametersWith(_parameters);
}

/******************************************************************************/

void HomogeneousTreeLikelihood::initParameters()
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

void HomogeneousTreeLikelihood::initBranchLengthsParameters()
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

void HomogeneousTreeLikelihood::ignoreParameter(const string & name)
throw (ParameterNotFoundException)
{
	_parameters.deleteParameter(name);
}

/******************************************************************************/

const SubstitutionModel * HomogeneousTreeLikelihood::getSubstitutionModel() const {
	return _model;
}

/******************************************************************************/

SubstitutionModel * HomogeneousTreeLikelihood::getSubstitutionModel() {
	return _model;
}

/******************************************************************************/

const DiscreteDistribution * HomogeneousTreeLikelihood::getRateDistribution() const {
	return _rateDistribution;
}

/******************************************************************************/

DiscreteDistribution * HomogeneousTreeLikelihood::getRateDistribution() {
	return _rateDistribution;
}

/******************************************************************************/

void HomogeneousTreeLikelihood::initTreeLikelihoods(Node * node, const SiteContainer & sequences) throw (Exception)
{
	unsigned int nbSites = _nbSites;

	//Initialize likelihood vector:
	VVVdouble * _likelihoods_node = & _likelihoods[node];
	VVVdouble * _dLikelihoods_node = & _dLikelihoods[node];
	VVVdouble * _d2Likelihoods_node = & _d2Likelihoods[node];
	
	_likelihoods_node -> resize(nbSites);
	_dLikelihoods_node -> resize(nbSites);
	_d2Likelihoods_node -> resize(nbSites);

	for(unsigned int i = 0; i < nbSites; i++) {
		VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i];
		VVdouble * _dLikelihoods_node_i = & (* _dLikelihoods_node)[i];
		VVdouble * _d2Likelihoods_node_i = & (* _d2Likelihoods_node)[i];
		_likelihoods_node_i -> resize(_nbClasses);
		_dLikelihoods_node_i -> resize(_nbClasses);
		_d2Likelihoods_node_i -> resize(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c];
			Vdouble * _dLikelihoods_node_i_c = & (* _dLikelihoods_node_i)[c];
			Vdouble * _d2Likelihoods_node_i_c = & (* _d2Likelihoods_node_i)[c];
			_likelihoods_node_i_c -> resize(_nbStates);
			_dLikelihoods_node_i_c -> resize(_nbStates);
			_d2Likelihoods_node_i_c -> resize(_nbStates);
			for(unsigned int s = 0; s < _nbStates; s++) {
				(* _likelihoods_node_i_c)[s] = 1; //All likelihoods are initialized to 1.
				(* _dLikelihoods_node_i_c)[s] = 0; //All dLikelihoods are initialized to 0.
				(* _d2Likelihoods_node_i_c)[s] = 0; //All d2Likelihoods are initialized to 0.
			}
		}
	}

	//Now initialize likelihood values and pointers:
	
	if(node -> isLeaf()) {
		for(unsigned int i = 0; i < nbSites; i++) {
			VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i]; 
			for(unsigned int c = 0; c < _nbClasses; c++) {
				Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c]; 
				for(unsigned int s = 0; s < _nbStates; s++) {
					//Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
					//otherwise value set to 0:
					try {
						//cout << "i=" << i << "\tc=" << c << "\ts=" << s << endl;
						int state = sequences.getSequence(node -> getName()) -> getValue(i);
						(* _likelihoods_node_i_c)[s] = _model -> getInitValue(s, state);
					} catch (SequenceNotFoundException snfe) {
						throw SequenceNotFoundException("HomogeneousTreeLikelihood::initTreelikelihoods. Leaf name in tree not found in site conainer: ", (node -> getName()));
					}
				}
			}
		}
	} else {
		//'node' is an internal node.
		map<Node *, vector<unsigned int> > * _patternLinks_node = & _patternLinks[node];
		unsigned int nbSonNodes = node -> getNumberOfSons();
		for(unsigned int l = 0; l < nbSonNodes; l++) {
			//For each son node,
			Node * son = (* node)[l];
			initTreeLikelihoods(son, sequences);
			vector<unsigned int> * _patternLinks_node_son = & _patternLinks[node][son];

			//Init map:
			_patternLinks_node_son -> resize(nbSites);

			for(unsigned int i = 0; i < nbSites; i++) {
				(* _patternLinks_node_son)[i] = i;
			}
		}
	}
}

/******************************************************************************/

SiteContainer * HomogeneousTreeLikelihood::initTreeLikelihoodsWithPatterns( Node * node, const SiteContainer & sequences) throw (Exception)
{
	SiteContainer * tmp = PatternTools::getSequenceSubset(sequences, * node);
	SiteContainer * subSequences = PatternTools::shrinkSiteSet(* tmp);
	delete tmp;

	unsigned int nbSites = subSequences -> getNumberOfSites();

	//Initialize likelihood vector:
	VVVdouble * _likelihoods_node = & _likelihoods[node];
	VVVdouble * _dLikelihoods_node = & _dLikelihoods[node];
	VVVdouble * _d2Likelihoods_node = & _d2Likelihoods[node];
	_likelihoods_node -> resize(nbSites);
	_dLikelihoods_node -> resize(nbSites);
	_d2Likelihoods_node -> resize(nbSites);

	for(unsigned int i = 0; i < nbSites; i++) {
		VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i];
		VVdouble * _dLikelihoods_node_i = & (* _dLikelihoods_node)[i];
		VVdouble * _d2Likelihoods_node_i = & (* _d2Likelihoods_node)[i];
		_likelihoods_node_i -> resize(_nbClasses);
		_dLikelihoods_node_i -> resize(_nbClasses);
		_d2Likelihoods_node_i -> resize(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c];
			Vdouble * _dLikelihoods_node_i_c = & (* _dLikelihoods_node_i)[c];
			Vdouble * _d2Likelihoods_node_i_c = & (* _d2Likelihoods_node_i)[c];
			_likelihoods_node_i_c -> resize(_nbStates);
			_dLikelihoods_node_i_c -> resize(_nbStates);
			_d2Likelihoods_node_i_c -> resize(_nbStates);
			for(unsigned int s = 0; s < _nbStates; s++) {
				(* _likelihoods_node_i_c)[s] = 1; //All likelihoods are initialized to 1.
				(* _dLikelihoods_node_i_c)[s] = 0; //All dLikelihoods are initialized to 0.
				(* _d2Likelihoods_node_i_c)[s] = 0; //All d2Likelihoods are initialized to 0.
			}
		}
	}

	//Now initialize likelihood values and pointers:
	
	//For efficiency, we make a copy of the data for direct accesss to sequences:
	const AlignedSequenceContainer asc(* subSequences);
	
	if(node -> isLeaf()) {
		for(unsigned int i = 0; i < nbSites; i++) {
			VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i];
			for(unsigned int c = 0; c < _nbClasses; c++) {
				Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c];
				for(unsigned int s = 0; s < _nbStates; s++) {
					//Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
					//otherwise value set to 0:
					try {
						//cout << "i=" << i << "\tc=" << c << "\ts=" << s << endl;
						int state = asc.getSequence(node -> getName()) -> getValue(i);
						(* _likelihoods_node_i_c)[s] = _model -> getInitValue(s, state);
					} catch (SequenceNotFoundException snfe) {
						throw SequenceNotFoundException("HomogeneousTreeLikelihood::initTreelikelihoodsWithPatterns. Leaf name in tree not found in site conainer: ", (node -> getName()));
					}
				}
			}
		}
	} else {
		//'node' is an internal node.
		map<Node *, vector<unsigned int> > * _patternLinks_node = & _patternLinks[node];
		
		//Now initialize pattern links:
		unsigned int nbSonNodes = node -> getNumberOfSons();
		for(unsigned int l = 0; l < nbSonNodes; l++) {
			//For each son node,
			Node * son = (* node)[l];

			vector<unsigned int> * _patternLinks_node_son = & _patternLinks[node][son];
			
			//Init map:
			_patternLinks_node_son -> resize(nbSites);

			//Initialize subtree 'l' and retrieves corresponding subSequences:
			SiteContainer * subSubSequencesShrunk = initTreeLikelihoodsWithPatterns(son, sequences);
			SiteContainer * subSubSequencesExpanded = PatternTools::getSequenceSubset(* subSequences, * son);

			for(unsigned int i = 0; i < nbSites; i++) {
				for(unsigned int ii = 0; ii < subSubSequencesShrunk -> getNumberOfSites(); ii++) {
					if(SiteTools::areSitesIdentical(* subSubSequencesShrunk -> getSite(ii), (* subSubSequencesExpanded -> getSite(i)))) {
						(* _patternLinks_node_son)[i] = ii;
					}
				}
			}
			delete subSubSequencesShrunk;
			delete subSubSequencesExpanded;
		}
	}
	//displayLikelihood(node);
	return subSequences;
}

/******************************************************************************/

void HomogeneousTreeLikelihood::computeSubtreeLikelihood(Node * node)
{
	if(node -> isLeaf()) return;

	unsigned int nbSites = _likelihoods[node].size();
	unsigned int nbNodes = node -> getNumberOfSons();
		
	// Must reset the likelihood array first (i.e. set all of them to 1):
	VVVdouble * _likelihoods_node = & _likelihoods[node];
	for(unsigned int i = 0; i < nbSites; i++) {
		//For each site in the sequence,
		VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i];
		for(unsigned int c = 0; c < _nbClasses; c++) {
			//For each rate classe,
			Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c];
			for(unsigned int x = 0; x < _nbStates; x++) {
				//For each initial state,
				(* _likelihoods_node_i_c)[x] = 1.;
			}
		}
	}

	for(unsigned int l = 0; l < nbNodes; l++) {
		//For each son node,	

		Node * son = node -> getSon(l);
		
		computeSubtreeLikelihood(son); //Recursive method:
		
		VVVdouble * _pxy_son = & _pxy[son];
		vector <unsigned int> * _patternLinks_node_son = & _patternLinks[node][son];
		VVVdouble * _likelihoods_son = & _likelihoods[son];

		for(unsigned int i = 0; i < nbSites; i++) {
			//For each site in the sequence,
			VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[(* _patternLinks_node_son)[i]];
			VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i];
			for(unsigned int c = 0; c < _nbClasses; c++) {
				//For each rate classe,
				Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
				Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c];
				VVdouble * _pxy_son_c = & (* _pxy_son)[c];
				for(unsigned int x = 0; x < _nbStates; x++) {
					//For each initial state,
					Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
					double likelihood = 0;
					for(unsigned int y = 0; y < _nbStates; y++) {
						likelihood += (* _pxy_son_c_x)[y] * (* _likelihoods_son_i_c)[y];
					}
					(* _likelihoods_node_i_c)[x] *= likelihood;
				}
			}
		}
	}
}

/******************************************************************************/

void HomogeneousTreeLikelihood::displayLikelihood(Node * node)
{
	cout << "Likelihoods at node " << node -> getName() << ": " << endl;
	for(unsigned int i = 0; i < _likelihoods[node].size(); i++) {
		cout << "Site " << i << ":" << endl;
		for(unsigned int c = 0; c < _likelihoods[node][i].size(); c++) {
			cout << "Rate class " << c;
			for(unsigned int s = 0; s < _likelihoods[node][i][c].size(); s++) {
				cout << "\t" << _likelihoods[node][i][c][s];
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << "                                         ***" << endl;
}

/*******************************************************************************/
