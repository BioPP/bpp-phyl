//
// File: HomogeneousTreeLikelihood.cpp
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

#include "HomogeneousTreeLikelihood.h"
#include "PatternTools.h"

//From SeqLib:
#include <Seq/SiteTools.h>
#include <Seq/AlignedSequenceContainer.h>
#include <Seq/SequenceContainerTools.h>

// From Utils:
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>

// From NumCalc:
using namespace VectorFunctions;

// From the STL:
#include <iostream>
using namespace std;

/******************************************************************************/

void DRASRTreeLikelihoodData::init(const SiteContainer & sites, const SubstitutionModel & model) throw (Exception)
{
	_alphabet = sites.getAlphabet();
	_nbStates = sites.getAlphabet() -> getSize();
 	_nbSites  = sites.getNumberOfSites();
	_shrunkData =	initTreeLikelihoodsWithPatterns(_tree -> getRootNode(), sites, model);
	_nbDistinctSites = _shrunkData -> getNumberOfSites();
	//Initialize root patterns:
	_rootPatternLinks.resize(_nbSites);
	for(unsigned int i = 0; i < _nbSites; i++) {
		const Site * site1 =  sites.getSite(i);
		for(unsigned int ii = 0; ii < _nbDistinctSites; ii++) {
			if(SiteTools::areSitesIdentical(* _shrunkData -> getSite(ii), * site1)) {
				//_rootPatternLinks[i] = & _likelihoods[_tree.getRootNode()][ii];
				_rootPatternLinks[i] = ii;
				break;
			}
		}
	}
}

/******************************************************************************/

void DRASRTreeLikelihoodData::initTreeLikelihoods(const Node * node, const SiteContainer & sequences, const SubstitutionModel & model) throw (Exception)
{
	//Initialize likelihood vector:
	VVVdouble * _likelihoods_node = & _nodeData[node].getLikelihoodArray();
	VVVdouble * _dLikelihoods_node = & _nodeData[node].getDLikelihoodArray();
	VVVdouble * _d2Likelihoods_node = & _nodeData[node].getD2LikelihoodArray();
	
	_likelihoods_node -> resize(_nbSites);
	_dLikelihoods_node -> resize(_nbSites);
	_d2Likelihoods_node -> resize(_nbSites);

	for(unsigned int i = 0; i < _nbSites; i++) {
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
		const Sequence * seq;
		try {
			seq = sequences.getSequence(node -> getName());
		} catch (SequenceNotFoundException snfe) {
			throw SequenceNotFoundException("HomogeneousTreeLikelihood::initTreelikelihoods. Leaf name in tree not found in site conainer: ", (node -> getName()));
		}	
		for(unsigned int i = 0; i < _nbSites; i++) {
			VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i]; 
			for(unsigned int c = 0; c < _nbClasses; c++) {
				Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c]; 
				for(unsigned int s = 0; s < _nbStates; s++) {
					//Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
					//otherwise value set to 0:
					//cout << "i=" << i << "\tc=" << c << "\ts=" << s << endl;
					int state = seq -> getValue(i);
					(* _likelihoods_node_i_c)[s] = model.getInitValue(s, state);
				}
			}
		}
	} else {
		//'node' is an internal node.
		map<const Node *, vector<unsigned int> > * _patternLinks_node = & _patternLinks[node];
		unsigned int nbSonNodes = node -> getNumberOfSons();
		for(unsigned int l = 0; l < nbSonNodes; l++) {
			//For each son node,
			const Node * son = (* node)[l];
			initTreeLikelihoods(son, sequences, model);
			vector<unsigned int> * _patternLinks_node_son = & (* _patternLinks_node)[son];

			//Init map:
			_patternLinks_node_son -> resize(_nbSites);

			for(unsigned int i = 0; i < _nbSites; i++) {
				(* _patternLinks_node_son)[i] = i;
			}
		}
	}
}

/******************************************************************************/

SiteContainer * DRASRTreeLikelihoodData::initTreeLikelihoodsWithPatterns(const Node * node, const SiteContainer & sequences, const SubstitutionModel & model) throw (Exception)
{
	SiteContainer * tmp = PatternTools::getSequenceSubset(sequences, * node);
	SiteContainer * subSequences = PatternTools::shrinkSiteSet(* tmp);
	delete tmp;

	unsigned int nbSites = subSequences -> getNumberOfSites();
	//_likelihoodData._tree = _tree;
	
	//Initialize likelihood vector:
	VVVdouble * _likelihoods_node = & _nodeData[node].getLikelihoodArray();
	VVVdouble * _dLikelihoods_node = & _nodeData[node].getDLikelihoodArray();
	VVVdouble * _d2Likelihoods_node = & _nodeData[node].getD2LikelihoodArray();
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
		const Sequence * seq;
		try {
			seq = subSequences -> getSequence(node -> getName());
		} catch (SequenceNotFoundException snfe) {
			throw SequenceNotFoundException("HomogeneousTreeLikelihood::initTreelikelihoodsWithPatterns. Leaf name in tree not found in site conainer: ", (node -> getName()));
		}	
		for(unsigned int i = 0; i < nbSites; i++) {
			VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i];
			for(unsigned int c = 0; c < _nbClasses; c++) {
				Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c];
				for(unsigned int s = 0; s < _nbStates; s++) {
					//Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
					//otherwise value set to 0:
					//cout << "i=" << i << "\tc=" << c << "\ts=" << s << endl;
					int state = seq -> getValue(i);
					(* _likelihoods_node_i_c)[s] = model.getInitValue(s, state);
				}
			}
		}
	} else {
		//'node' is an internal node.
		map<const Node *, vector<unsigned int> > * _patternLinks_node = & _patternLinks[node];
		
		//Now initialize pattern links:
		unsigned int nbSonNodes = node -> getNumberOfSons();
		for(unsigned int l = 0; l < nbSonNodes; l++) {
			//For each son node,
			const Node * son = (* node)[l];

			vector<unsigned int> * _patternLinks_node_son = & (* _patternLinks_node)[son];
			
			//Init map:
			_patternLinks_node_son -> resize(nbSites);

			//Initialize subtree 'l' and retrieves corresponding subSequences:
			SiteContainer * subSubSequencesShrunk = initTreeLikelihoodsWithPatterns(son, sequences, model);
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

HomogeneousTreeLikelihood::HomogeneousTreeLikelihood(
	TreeTemplate<Node> & tree,
	const SiteContainer & data,
	SubstitutionModel * model,
	DiscreteDistribution * rDist,
	bool verbose
)	throw (Exception):
	AbstractDiscreteRatesAcrossSitesTreeLikelihood(rDist, verbose), // We must do this since AbstractTreeLikelihood is virtual
	AbstractHomogeneousTreeLikelihood(tree, data, model, rDist, verbose)
{
	//Init _likelihoods:
	if(verbose) ApplicationTools::message << "Homogeneous Tree Likelihood" << endl;	
	_likelihoodData = new DRASRTreeLikelihoodData(* _tree, rDist -> getNumberOfCategories());
	
	if(verbose) ApplicationTools::displayTask("Init likelihoods data structure");
	_likelihoodData -> init(*_data, *_model);
	if(verbose) ApplicationTools::displayTaskDone();
	
	if(verbose) ApplicationTools::displayResult("Number of distinct sites",
			TextTools::toString(_likelihoodData -> getNumberOfDistinctSites()));
		
	// Now initializes all parameters:
	initParameters();
	fireParameterChanged(_parameters);
}

/******************************************************************************/

HomogeneousTreeLikelihood::~HomogeneousTreeLikelihood() { delete _likelihoodData; }

/******************************************************************************/

double HomogeneousTreeLikelihood::getLikelihood() const
{
	double l = 1.;
	for(unsigned int i = 0; i < _nbSites; i++) {
		l *= getLikelihoodForASite(i);
	}
	return l;
}

/******************************************************************************/

double HomogeneousTreeLikelihood::getLogLikelihood() const
{
	double ll = 0;
	for(unsigned int i = 0; i < _nbSites; i++) {
		ll += getLogLikelihoodForASite(i);
	}
	return ll;
}

/******************************************************************************/

double HomogeneousTreeLikelihood::getLikelihoodForASite(unsigned int site) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbClasses; i++) {
		l += getLikelihoodForASiteForARateClass(site, i) * _rateDistribution -> getProbability(i);
	}
	return l;
}

/******************************************************************************/

double HomogeneousTreeLikelihood::getLogLikelihoodForASite(unsigned int site) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbClasses; i++) {
		l += getLikelihoodForASiteForARateClass(site, i) * _rateDistribution -> getProbability(i);
	}
	//if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
	if(l<0) l=0; //May happen because of numerical errors.
	return log(l);
}

/******************************************************************************/

double HomogeneousTreeLikelihood::getLikelihoodForASiteForAState(unsigned int site, int state) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbClasses; i++) {
		l += getLikelihoodForASiteForARateClassForAState(site, i, state) * _rateDistribution -> getProbability(i);
	}
	return l;
}

/******************************************************************************/

double HomogeneousTreeLikelihood::getLogLikelihoodForASiteForAState(unsigned int site, int state) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbClasses; i++) {
		l += getLikelihoodForASiteForARateClassForAState(site, i, state) * _rateDistribution -> getProbability(i);
	}
	//if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
	return log(l);
}

/******************************************************************************/

double HomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
	double l = 0;
	Vdouble * la = & _likelihoodData -> getLikelihoodArray(_tree -> getRootNode())[_likelihoodData->getRootArrayPosition(site)][rateClass];
	for(unsigned int i = 0; i < _nbStates; i++) {
		l += (* la)[i] * _model -> freq(i);
	}
	return l;
}

/******************************************************************************/

double HomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
	double l = 0;
	Vdouble * la = & _likelihoodData -> getLikelihoodArray(_tree -> getRootNode())[_likelihoodData->getRootArrayPosition(site)][rateClass];
	for(unsigned int i = 0; i < _nbStates; i++) {
		l += (* la)[i] * _model -> freq(i);
	}
	//if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
	return log(l);
}

/******************************************************************************/	

double HomogeneousTreeLikelihood::getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
	//return _likelihoods[_tree -> getRootNode()][_rootPatternLinks[site]][rateClass][state];
	return _likelihoodData -> getLikelihoodArray(_tree -> getRootNode())[_likelihoodData->getRootArrayPosition(site)][rateClass][state];
}

/******************************************************************************/

double HomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
	//return log(_likelihoods[_tree -> getRootNode()][_rootPatternLinks[site]][rateClass][state]);
	return log(_likelihoodData->getLikelihoodArray(_tree -> getRootNode())[_likelihoodData->getRootArrayPosition(site)][rateClass][state]);
}

/******************************************************************************/	


VVdouble HomogeneousTreeLikelihood::getPosteriorProbabilitiesOfEachRate() const
{
	VVdouble pb = getLikelihoodForEachSiteForEachRateClass();
	Vdouble  l  = getLikelihoodForEachSite();
	for(unsigned int i = 0; i < _nbSites; i++) {
		for(unsigned int j = 0; j < _nbClasses; j++) pb[i][j] = pb[i][j] * _rateDistribution -> getProbability(j) / l[i]; 
	}
	return pb;
}
	
/******************************************************************************/

Vdouble HomogeneousTreeLikelihood::getRateWithMaxPostProbOfEachSite() const
{
	VVdouble l = getLikelihoodForEachSiteForEachRateClass();
	Vdouble rates(_nbSites);
	for(unsigned int i = 0; i < _nbSites; i++) {
		rates[i] = _rateDistribution -> getCategory(whichmax<double>(l[i]));
	}
	return rates;
}

/******************************************************************************/

Vdouble HomogeneousTreeLikelihood::getPosteriorRateOfEachSite() const
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
	
		//Computes all pxy and pyx once for all:
		VVVdouble * _pxy_son = & _pxy[son];
		_pxy_son -> resize(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			VVdouble * _pxy_son_c = & (* _pxy_son)[c];
			_pxy_son_c -> resize(_nbStates);
			RowMatrix<double> Q = _model -> getPij_t(l * _rateDistribution -> getCategory(c));
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
			RowMatrix<double> dQ = _model -> getdPij_dt(l * rc);  
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
			RowMatrix<double> d2Q = _model -> getd2Pij_dt2(l * rc);
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

double HomogeneousTreeLikelihood::getDLikelihoodForASiteForARateClass(
	unsigned int site,
	unsigned int rateClass) const
{
	double dl = 0;
	Vdouble * dla = & _likelihoodData->getDLikelihoodArray(_tree -> getRootNode())[_likelihoodData->getRootArrayPosition(site)][rateClass];
	for(unsigned int i = 0; i < _nbStates; i++) {
		dl += (* dla)[i] * _model -> freq(i);
	}
	return dl;
}

/******************************************************************************/	

double HomogeneousTreeLikelihood::getDLikelihoodForASite(unsigned int site) const
{
	// Derivative of the sum is the sum of derivatives:
	double dl = 0;
	for(unsigned int i = 0; i < _nbClasses; i++)
		dl += getDLikelihoodForASiteForARateClass(site, i) * _rateDistribution -> getProbability(i);
	return dl;
}

/******************************************************************************/	

double HomogeneousTreeLikelihood::getDLogLikelihoodForASite(unsigned int site) const
{
	// d(f(g(x)))/dx = dg(x)/dx . df(g(x))/dg :
	return getDLikelihoodForASite(site) / getLikelihoodForASite(site);
}

/******************************************************************************/	

double HomogeneousTreeLikelihood::getDLogLikelihood() const
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
	const Node * branch = _nodes[brI];
	const Node * father = branch -> getFather();
	//VVVdouble * _dLikelihoods_father = & _dLikelihoods[father];
	VVVdouble * _dLikelihoods_father = & _likelihoodData->getDLikelihoodArray(father);
	
	// Compute dLikelihoods array for the father node.
	// Fist initialize to 1:
	unsigned int nbSites = _dLikelihoods_father -> size();
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
		
		const Node * son = father -> getSon(l);

		vector <unsigned int> * _patternLinks_father_son = & _likelihoodData -> getArrayPositions(father, son);
		VVVdouble * _likelihoods_son = & _likelihoodData -> getLikelihoodArray(son);

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

void HomogeneousTreeLikelihood::computeDownSubtreeDLikelihood(const Node * node)
{
	const Node * father = node -> getFather();
	// We assume that the _dLikelihoods array has been filled for the current node 'node'.
	// We will evaluate the array for the father node.
	if(father == NULL) return; // We reached the root!
		
	// Compute dLikelihoods array for the father node.
	// Fist initialize to 1:
	VVVdouble * _dLikelihoods_father = & _likelihoodData->getDLikelihoodArray(father);
	unsigned int nbSites = _dLikelihoods_father -> size();
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
		const Node * son = father -> getSon(l);

		VVVdouble * _pxy_son = & _pxy[son];
		vector <unsigned int> * _patternLinks_father_son = & _likelihoodData->getArrayPositions(father,son);
	
		if(son == node) {
			VVVdouble * _dLikelihoods_son = & _likelihoodData->getDLikelihoodArray(son);
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
			VVVdouble * _likelihoods_son = & _likelihoodData->getLikelihoodArray(son);
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

double HomogeneousTreeLikelihood::getD2LikelihoodForASiteForARateClass(
	unsigned int site,
	unsigned int rateClass) const
{
	double d2l = 0;
	Vdouble * d2la = & _likelihoodData->getD2LikelihoodArray(_tree -> getRootNode())[_likelihoodData->getRootArrayPosition(site)][rateClass];
	for(unsigned int i = 0; i < _nbStates; i++) {
		d2l += (* d2la)[i] * _model -> freq(i);
	}
	return d2l;
}

/******************************************************************************/	

double HomogeneousTreeLikelihood::getD2LikelihoodForASite(unsigned int site) const
{
	// Derivative of the sum is the sum of derivatives:
	double d2l = 0;
	for(unsigned int i = 0; i < _nbClasses; i++)
		d2l += getD2LikelihoodForASiteForARateClass(site, i) * _rateDistribution -> getProbability(i);
	return d2l;
}

/******************************************************************************/	

double HomogeneousTreeLikelihood::getD2LogLikelihoodForASite(unsigned int site) const
{
	return getD2LikelihoodForASite(site) / getLikelihoodForASite(site)
	- pow( getDLikelihoodForASite(site) / getLikelihoodForASite(site), 2);
}

/******************************************************************************/	

double HomogeneousTreeLikelihood::getD2LogLikelihood() const
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
	VVVdouble * _d2Likelihoods_father = & _likelihoodData->getD2LikelihoodArray(father); 
	unsigned int nbSites = _d2Likelihoods_father -> size();
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
		
		const Node * son = father -> getSon(l);
		
		vector <unsigned int> * _patternLinks_father_son = & _likelihoodData->getArrayPositions(father,son);
		VVVdouble * _likelihoods_son = & _likelihoodData->getLikelihoodArray(son);

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

void HomogeneousTreeLikelihood::computeDownSubtreeD2Likelihood(const Node * node)
{
	const Node * father = node -> getFather();
	// We assume that the _dLikelihoods array has been filled for the current node 'node'.
	// We will evaluate the array for the father node.
	if(father == NULL) return; // We reached the root!
		
	// Compute dLikelihoods array for the father node.
	// Fist initialize to 1:
	VVVdouble * _d2Likelihoods_father = & _likelihoodData -> getD2LikelihoodArray(father);
	unsigned int nbSites = _d2Likelihoods_father -> size();
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
		const Node * son = father -> getSon(l);

		VVVdouble * _pxy_son = & _pxy[son];
		vector <unsigned int> * _patternLinks_father_son = & _likelihoodData->getArrayPositions(father,son);
	
		if(son == node) {
			VVVdouble * _d2Likelihoods_son = & _likelihoodData->getD2LikelihoodArray(son);
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
			VVVdouble * _likelihoods_son = & _likelihoodData->getLikelihoodArray(son);
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
	
void HomogeneousTreeLikelihood::computeTreeLikelihood()
{
	computeSubtreeLikelihood(_tree -> getRootNode());
}

/******************************************************************************/	

void HomogeneousTreeLikelihood::computeSubtreeLikelihood(const Node * node)
{
	if(node -> isLeaf()) return;

	unsigned int nbSites = _likelihoodData -> getLikelihoodArray(node).size();
	unsigned int nbNodes = node -> getNumberOfSons();
		
	// Must reset the likelihood array first (i.e. set all of them to 1):
	VVVdouble * _likelihoods_node = & _likelihoodData -> getLikelihoodArray(node);
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

		const Node * son = node -> getSon(l);
		
		computeSubtreeLikelihood(son); //Recursive method:
		
		VVVdouble * _pxy_son = & _pxy[son];
		vector <unsigned int> * _patternLinks_node_son = & _likelihoodData -> getArrayPositions(node,son);
		VVVdouble * _likelihoods_son = & _likelihoodData -> getLikelihoodArray(son);

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

void HomogeneousTreeLikelihood::displayLikelihood(const Node * node)
{
	cout << "Likelihoods at node " << node -> getName() << ": " << endl;
	displayLikelihoodArray(_likelihoodData -> getLikelihoodArray(node));
	cout << "                                         ***" << endl;
}

/*******************************************************************************/

