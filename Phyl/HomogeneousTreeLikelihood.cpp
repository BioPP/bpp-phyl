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
}

/******************************************************************************/

HomogeneousTreeLikelihood::~HomogeneousTreeLikelihood() { delete _data; }

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
	return log(l);
}

/******************************************************************************/

double HomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
	double l = 0;
	for(unsigned int i = 0; i < _nbStates; i++) {
		//l += _rootPatternLinks[site] -> operator[](rateClass)[i] * _model -> freq(i);
		l += _likelihoods[_tree -> getRootNode()][_rootPatternLinks[site]][rateClass][i] * _model -> freq(i);
	}
	return l;
}

/******************************************************************************/

double HomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
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

double HomogeneousTreeLikelihood::f(const ParameterList & parameters) const
{
	try { 
		const_cast<HomogeneousTreeLikelihood *>(this) -> setParametersValues(parameters);
	} catch(ConstraintException ce) {
		return -log(0.); // (+inf if unlikely!)
	}
	const_cast<HomogeneousTreeLikelihood *>(this) -> applyParameters();
	const_cast<HomogeneousTreeLikelihood *>(this) -> computeTreeLikelihood();
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
	for(unsigned int i = 0; i < _nbStates; i++) {
		//dl += _rootPatternLinks[site] -> operator[](rateClass)[i] * _model -> freq(i);
		dl += _likelihoods[_tree -> getRootNode()][_rootPatternLinks[site]][rateClass][i] * _model -> freq(i);
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
	return 1. / getLikelihoodForASite(site) * getDLikelihoodForASite(site);
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

double HomogeneousTreeLikelihood::df (
	const string & variable,
	const ParameterList & parameters) const
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

void HomogeneousTreeLikelihood::computeTreeDLikelihood(const string & variable) {

	// First reset all dLikelihoods arrays to 0:
	for(unsigned int l = 0; l < _nbNodes; l++) {
		unsigned int nbSites = _likelihoods[_nodes[l]].size();
		for(unsigned int i = 0; i < nbSites; i++) {
			for(unsigned int c = 0; c < _nbClasses; c++) {
				for(unsigned int s = 0; s < _nbStates; s++) {
					_dLikelihoods[_nodes[l]][i][c][s] = 0;
				}
			}
		}
	}
	
	// Then get the node with the branch whose length must be derivated:
	int brI = TextTools::toInt(variable.substr(5));
	Node * branch = _nodes[brI];
	Node * father = branch -> getFather();
	
	// Compute dLikelihoods array for the father node.
	// Fist initialize to 1:
	unsigned int nbSites = _likelihoods[father].size();
	for(unsigned int i = 0; i < nbSites; i++) {
		for(unsigned int c = 0; c < _nbClasses; c++) {
			for(unsigned int s = 0; s < _nbStates; s++) {
				_dLikelihoods[father][i][c][s] = 1.;	
			}
		}
	}

	unsigned int nbNodes = father -> getNumberOfSons();
	for(unsigned int l = 0; l < nbNodes; l++) {
		Node * son = father -> getSon(l);

		//Computes all pxy once for all:
		VVVdouble pxy = VVVdouble(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			pxy[c] = VVdouble(_nbStates);
			for(unsigned int x = 0; x < _nbStates; x++) {
				pxy[c][x] = Vdouble(_nbStates);
				for(unsigned int y = 0; y < _nbStates; y++) {
					pxy[c][x][y] = _model -> Pij_t(x, y, son -> getDistanceToFather() * _rateDistribution -> getCategory(c));
				}
			}
		}

		//Computes all dpxy/dt once for all:
		VVVdouble dpxy = VVVdouble(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			dpxy[c] = VVdouble(_nbStates);
			for(unsigned int x = 0; x < _nbStates; x++) {
				dpxy[c][x] = Vdouble(_nbStates);
				for(unsigned int y = 0; y < _nbStates; y++) {
					dpxy[c][x][y] = _model -> dPij_dt(x, y, son -> getDistanceToFather() * _rateDistribution -> getCategory(c));
				}
			}
		}

		for(unsigned int i = 0; i < nbSites; i++) {
			for(unsigned int c = 0; c < _nbClasses; c++) {
				for(unsigned int x = 0; x < _nbStates; x++) {
					double dl = 0;
					for(unsigned int y = 0; y < _nbStates; y++) {
						dl += (son == branch ?	dpxy[c][x][y] : pxy[c][x][y])
						   * _likelihoods[son][_patternLinks[father][son][i]][c][y];
					}
					_dLikelihoods[father][i][c][x] *= dl;
				}
			}
		}
	}
	computeDownSubtreeDLikelihood(father);
}

/******************************************************************************/

void HomogeneousTreeLikelihood::computeDownSubtreeDLikelihood(Node * node) {
	Node * father = node -> getFather();
	if(father == NULL) return;
		
	// Compute dLikelihoods array for the father node.
	// Fist initialize to 1:
	unsigned int nbSites = _likelihoods[father].size();
	for(unsigned int i = 0; i < nbSites; i++) {
		for(unsigned int c = 0; c < _nbClasses; c++) {
			for(unsigned int s = 0; s < _nbStates; s++) {
				_dLikelihoods[father][i][c][s] = 1.;	
			}
		}
	}

	unsigned int nbNodes = father -> getNumberOfSons();
	for(unsigned int l = 0; l < nbNodes; l++) {
		Node * son = father -> getSon(l);

		//Computes all pxy once for all:
		VVVdouble pxy = VVVdouble(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			pxy[c] = VVdouble(_nbStates);
			for(unsigned int x = 0; x < _nbStates; x++) {
				pxy[c][x] = Vdouble(_nbStates);
				for(unsigned int y = 0; y < _nbStates; y++) {
					pxy[c][x][y] = _model -> Pij_t(x, y, son -> getDistanceToFather() * _rateDistribution -> getCategory(c));
				}
			}
		}

		for(unsigned int i = 0; i < nbSites; i++) {
			for(unsigned int c = 0; c < _nbClasses; c++) {
				for(unsigned int x = 0; x < _nbStates; x++) {
					double dl = 0;
					for(unsigned int y = 0; y < _nbStates; y++) {
						dl += pxy[c][x][y] * (son == node ?
							_dLikelihoods[son][_patternLinks[father][son][i]][c][y] :
							 _likelihoods[son][_patternLinks[father][son][i]][c][y]);
					}
					_dLikelihoods[father][i][c][x] *= dl;
				}
			}
		}
	}
	computeDownSubtreeDLikelihood(father);
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
	ParameterList pl;
	for(unsigned int i = 0; i < _nbNodes; i++) {
		double d = _nodes[i] -> getDistanceToFather();
		if (d <= 0) {
			cout << "WARNING!!! Branch length " << i << " is <=0. Value is set to 0.000001." << endl;
			_nodes[i] -> setDistanceToFather(0.000001);
			d = 0.000001;
		}
		pl.addParameter(Parameter("BrLen" + TextTools::toString(i), d, & Parameter::R_PLUS_STAR));
	}
	return pl;
}

/******************************************************************************/

ParameterList HomogeneousTreeLikelihood::getSubstitutionModelParameters() const {
	return _model -> getParameters();
}

/******************************************************************************/

ParameterList HomogeneousTreeLikelihood::getRateDistributionParameters() const {
	return _rateDistribution -> getParameters();
}

/******************************************************************************/

void HomogeneousTreeLikelihood::initParameters()
{
	// Reset parameters:
	_parameters.reset();
	
	// Branch lengths:
	_parameters.addParameters(getBranchLengthsParameters());
	
	// Substitution model:
	_parameters.addParameters(_model -> getParameters());
	
	// Rate distribution:
	_parameters.addParameters(_rateDistribution -> getParameters());
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
	_likelihoods[node].resize(nbSites);
	_dLikelihoods[node].resize(nbSites);

	for(unsigned int i = 0; i < nbSites; i++) {
		_likelihoods[node][i].resize(_nbClasses);
		_dLikelihoods[node][i].resize(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			_likelihoods[node][i][c].resize(_nbStates);
			_dLikelihoods[node][i][c].resize(_nbStates);
			for(unsigned int s = 0; s < _nbStates; s++) {
				_likelihoods[node][i][c][s] = 1; //All likelihoods are initialized to 1.
				_dLikelihoods[node][i][c][s] = 0; //All dLikelihoods are initialized to 0.
			}
		}
	}

	//Now initialize likelihood values and pointers:
	
	if(node -> isLeaf()) {
		for(unsigned int i = 0; i < nbSites; i++) {
			for(unsigned int c = 0; c < _nbClasses; c++) {
				for(unsigned int s = 0; s < _nbStates; s++) {
					//Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
					//otherwise value set to 0:
					try {
						//cout << "i=" << i << "\tc=" << c << "\ts=" << s << endl;
						int state = sequences.getSequence(node -> getName()) -> getValue(i);
						_likelihoods[node][i][c][s] = _model -> getInitValue(s, state);
					} catch (SequenceNotFoundException snfe) {
						throw SequenceNotFoundException("HomogeneousTreeLikelihood::initTreelikelihoods. Leaf name in tree not found in site conainer: ", (node -> getName()));
					}
				}
			}
		}
	} else {
		//'node' is an internal node.
		unsigned int nbSonNodes = node -> getNumberOfSons();
		for(unsigned int l = 0; l < nbSonNodes; l++) {
			//For each son node,
			Node * son = (* node)[l];
			initTreeLikelihoods(son, sequences);

			//Init map:
			_patternLinks[node][son].resize(nbSites);

			for(unsigned int i = 0; i < nbSites; i++) {
				//_patternLinks[node][son][i] = & _likelihoods[son][i];
				_patternLinks[node][son][i] = i;
			}
		}
	}
	//displayLikelihood(node);
}

/******************************************************************************/

SiteContainer * HomogeneousTreeLikelihood::initTreeLikelihoodsWithPatterns( Node * node, const SiteContainer & sequences) throw (Exception)
{
	SiteContainer * tmp = PatternTools::getSequenceSubset(sequences, * node);
	SiteContainer * subSequences = PatternTools::shrinkSiteSet(* tmp);
	delete tmp;

	unsigned int nbSites = subSequences -> getNumberOfSites();

	//Initialize likelihood vector:
	_likelihoods[node].resize(nbSites);
	_dLikelihoods[node].resize(nbSites);

	for(unsigned int i = 0; i < nbSites; i++) {
		_likelihoods[node][i].resize(_nbClasses);
		_dLikelihoods[node][i].resize(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			_likelihoods[node][i][c].resize(_nbStates);
			_dLikelihoods[node][i][c].resize(_nbStates);
			for(unsigned int s = 0; s < _nbStates; s++) {
				_likelihoods[node][i][c][s] = 1; //All likelihoods are initialized to 1.
				_dLikelihoods[node][i][c][s] = 0; //All dLikelihoods are initialized to 0.
			}
		}
	}

	//Now initialize likelihood values and pointers:
	
	//For efficiency, we make a copy of the data for direct accesss to sequences:
	const AlignedSequenceContainer asc(* subSequences);
	
	if(node -> isLeaf()) {
		for(unsigned int i = 0; i < nbSites; i++) {
			for(unsigned int c = 0; c < _nbClasses; c++) {
				for(unsigned int s = 0; s < _nbStates; s++) {
					//Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
					//otherwise value set to 0:
					try {
						//cout << "i=" << i << "\tc=" << c << "\ts=" << s << endl;
						int state = asc.getSequence(node -> getName()) -> getValue(i);
						_likelihoods[node][i][c][s] = _model -> getInitValue(s, state);
					} catch (SequenceNotFoundException snfe) {
						throw SequenceNotFoundException("HomogeneousTreeLikelihood::initTreelikelihoodsWithPatterns. Leaf name in tree not found in site conainer: ", (node -> getName()));
					}
				}
			}
		}
	} else {
		//'node' is an internal node.
		//Now initialize pattern links:

		unsigned int nbSonNodes = node -> getNumberOfSons();
		for(unsigned int l = 0; l < nbSonNodes; l++) {
			//For each son node,
			Node * son = (* node)[l];

			//Init map:
			_patternLinks[node][son].resize(nbSites);

			//Initialize subtree 'l' and retrieves corresponding subSequences:
			SiteContainer * subSubSequencesShrunk = initTreeLikelihoodsWithPatterns(son, sequences);
			SiteContainer * subSubSequencesExpanded = PatternTools::getSequenceSubset(* subSequences, * son);

			for(unsigned int i = 0; i < nbSites; i++) {
				for(unsigned int ii = 0; ii < subSubSequencesShrunk -> getNumberOfSites(); ii++) {
					if(SiteTools::areSitesIdentical(* subSubSequencesShrunk -> getSite(ii), (* subSubSequencesExpanded -> getSite(i)))) {
						//_patternLinks[node][son][i] = & _likelihoods[son][ii];
						_patternLinks[node][son][i] = ii;
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
	for(unsigned int i = 0; i < nbSites; i++) {
		//For each site in the sequence,
		for(unsigned int c = 0; c < _nbClasses; c++) {
			//For each rate classe,
			for(unsigned int x = 0; x < _nbStates; x++) {
				//For each initial state,
				_likelihoods[node][i][c][x] = 1.;
			}
		}
	}

	for(unsigned int l = 0; l < nbNodes; l++) {
		//For each son node,
		Node * son = node -> getSon(l);
		computeSubtreeLikelihood(son); //Recursive method:

		//Computes all pxy once for all:
		VVVdouble pxy = VVVdouble(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			pxy[c] = VVdouble(_nbStates);
			for(unsigned int x = 0; x < _nbStates; x++) {
				pxy[c][x] = Vdouble(_nbStates);
				for(unsigned int y = 0; y < _nbStates; y++) {
					pxy[c][x][y] = _model -> Pij_t(x, y, son -> getDistanceToFather() * _rateDistribution -> getCategory(c));
				}
			}
		}

		for(unsigned int i = 0; i < nbSites; i++) {
			//For each site in the sequence,
			for(unsigned int c = 0; c < _nbClasses; c++) {
				//For each rate classe,
				for(unsigned int x = 0; x < _nbStates; x++) {
					//For each initial state,
					double likelihood = 0;
					for(unsigned int y = 0; y < _nbStates; y++) {
						//likelihood += pxy[c][x][y] * (* _patternLinks[node][son][i])[c][y];
						likelihood += pxy[c][x][y] * _likelihoods[son][_patternLinks[node][son][i]][c][y];
					}
					//cout << "Node " << node -> getName() << ", x = " << x << ", L = " << likelihood << endl;
					_likelihoods[node][i][c][x] *= likelihood;
				}
			}
		}
	}
	//displayLikelihood(node);
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
