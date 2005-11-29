//
// File: DistanceEstimation.cpp
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Wed jun 08 10:39 2005
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

#include "DistanceEstimation.h"
#include "DistanceMatrix.h"
#include "Tree.h"
#include "PatternTools.h"
#include "SitePatterns.h"

// From Utils:
#include <Utils/ApplicationTools.h>

// From SeqLib:
#include <Seq/SiteTools.h>
#include <Seq/Sequence.h>
#include <Seq/AlignedSequenceContainer.h>

// From the STL:
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

/******************************************************************************/

TwoTreeLikelihood::TwoTreeLikelihood(
	const string & seq1, const string & seq2,	
	const SiteContainer & data,
	SubstitutionModel * model,
	DiscreteDistribution * rDist,
	bool verbose
)	throw (Exception):
	AbstractDiscreteRatesAcrossSitesTreeLikelihood(rDist, verbose)
{
	_seqnames.resize(2);
	_seqnames[0] = seq1;
	_seqnames[1] = seq2;
	_data = PatternTools::getSequenceSubset(data, _seqnames);
	_model = model;
	if(_data -> getAlphabet() -> getAlphabetType()
			!= _model -> getAlphabet() -> getAlphabetType())
		throw AlphabetMismatchException("TwoTreeTreeLikelihood::TwoTreeTreeLikelihood. Data and model must have the same alphabet type.",
				_data -> getAlphabet(),
				_model -> getAlphabet());

	_nbSites   = _data -> getNumberOfSites();
	_nbClasses = _rateDistribution -> getNumberOfCategories();
	_nbStates  = _model -> getAlphabet() -> getSize();	
	if(verbose) ApplicationTools::message << "Double-Recursive Homogeneous Tree Likelihood" << endl;	
	
	_brLen = 0.000001;

	//Initialize root patterns:
	SitePatterns pattern(* _data);
	_shrunkData = pattern.getSites();
	_rootWeights = pattern.getWeights();
	_rootPatternLinks = pattern.getIndices();
	_nbDistinctSites = _shrunkData -> getNumberOfSites();
	if(verbose) ApplicationTools::displayResult("Number of distinct sites",	TextTools::toString(_nbDistinctSites));
	
	
	//Init _likelihoods:
	if(verbose) ApplicationTools::displayTask("Init likelihoods arrays recursively");
	// Clone data for more efficiency on sequences access:
	const SiteContainer * sequences = new AlignedSequenceContainer(* _shrunkData);
	initTreeLikelihoods(* sequences);
	delete sequences;

	if(verbose) ApplicationTools::displayTaskDone();
	
	// Now initializes all parameters:
	initParameters();
	fireParameterChanged(_parameters);
}

/******************************************************************************/

TwoTreeLikelihood::~TwoTreeLikelihood()
{
	delete _data; 
}

/******************************************************************************/

ParameterList TwoTreeLikelihood::getBranchLengthsParameters() const
{
	return _brLenParameters.getCommonParametersWith(_parameters);
}

/******************************************************************************/

ParameterList TwoTreeLikelihood::getSubstitutionModelParameters() const
{
	return _model -> getParameters().getCommonParametersWith(_parameters);
}

/******************************************************************************/

double TwoTreeLikelihood::getLikelihood() const
{
	double l = 1.;
	for(unsigned int i = 0; i < _nbDistinctSites; i++) {
		l *= std::pow(_rootLikelihoodsSR[i], (int)_rootWeights[i]);
	}
	return l;
}

/******************************************************************************/

double TwoTreeLikelihood::getLogLikelihood() const
{
	double ll = 0;
	for(unsigned int i = 0; i < _nbDistinctSites; i++) {
		ll += _rootWeights[i] * log(_rootLikelihoodsSR[i]);
	}
	return ll;
}

/******************************************************************************/

double TwoTreeLikelihood::getLikelihoodForASite(unsigned int site) const
{
	return _rootLikelihoodsSR[_rootPatternLinks[site]];
}

/******************************************************************************/

double TwoTreeLikelihood::getLogLikelihoodForASite(unsigned int site) const
{
	return log(_rootLikelihoodsSR[_rootPatternLinks[site]]);
}

/******************************************************************************/

double TwoTreeLikelihood::getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
	return _rootLikelihoodsS[_rootPatternLinks[site]][rateClass];
}

/******************************************************************************/

double TwoTreeLikelihood::getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
	return log(_rootLikelihoodsS[_rootPatternLinks[site]][rateClass]);
}

/******************************************************************************/	

double TwoTreeLikelihood::getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
	return _rootLikelihoods[_rootPatternLinks[site]][rateClass][state];
}

/******************************************************************************/

double TwoTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
	return log(_rootLikelihoods[_rootPatternLinks[site]][rateClass][state]);
}

/******************************************************************************/

void TwoTreeLikelihood::initParameters()
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

void TwoTreeLikelihood::ignoreParameter(const string & name)
throw (ParameterNotFoundException)
{
	_parameters.deleteParameter(name);
}

/******************************************************************************/

void TwoTreeLikelihood::applyParameters() throw (Exception)
{
	//Apply branch length:
	_brLen = _parameters.getParameter(string("BrLen")) -> getValue();
	//Apply substitution model parameters:
	_model -> matchParametersValues(_parameters);
	//Apply rate distribution parameters:
	_rateDistribution -> matchParametersValues(_parameters);
}

/******************************************************************************/

void TwoTreeLikelihood::initBranchLengthsParameters()
{
	if (_brLen <= 0) {
		cout << "WARNING!!! Branch length is <=0. Value is set to 0.000001." << endl;
		_brLen = 0.000001;
	}
	_brLenParameters.reset();
	_brLenParameters.addParameter(Parameter("BrLen", _brLen, & Parameter::R_PLUS_STAR));
}

/******************************************************************************/

void TwoTreeLikelihood::setParameters(const ParameterList & parameters)
	throw (ParameterNotFoundException, ConstraintException)
{
	setParametersValues(parameters);
}

/******************************************************************************/

void TwoTreeLikelihood::fireParameterChanged(const ParameterList & params)
{
	applyParameters();

	// For now we ignore the parameter that changed and we recompute all arrays...

	//Computes all pxy and pyx once for all:
	_pxy.resize(_nbClasses);
	for(unsigned int c = 0; c < _nbClasses; c++) {
		VVdouble * _pxy_c = & _pxy[c];
		_pxy_c -> resize(_nbStates);
		RowMatrix<double> Q = _model -> getPij_t(_brLen * _rateDistribution -> getCategory(c));
		for(unsigned int x = 0; x < _nbStates; x++) {
			Vdouble * _pxy_c_x = & (* _pxy_c)[x];
			_pxy_c_x -> resize(_nbStates);
			for(unsigned int y = 0; y < _nbStates; y++) {
				(* _pxy_c_x)[y] = Q(x, y);
			}
		}
		
		if(_computeDerivatives) {

			//Computes all dpxy/dt once for all:
			_dpxy.resize(_nbClasses);
			for(unsigned int c = 0; c < _nbClasses; c++) {
				VVdouble * _dpxy_c = & _dpxy[c];
				_dpxy_c -> resize(_nbStates);
				double rc = _rateDistribution -> getCategory(c);
				RowMatrix<double> dQ = _model -> getdPij_dt(_brLen * rc);  
				for(unsigned int x = 0; x < _nbStates; x++) {
					Vdouble * _dpxy_c_x = & (* _dpxy_c)[x];
					_dpxy_c_x -> resize(_nbStates);
					for(unsigned int y = 0; y < _nbStates; y++) {
						(* _dpxy_c_x)[y] =  rc * dQ(x, y); 
					}
				}
			}
			
			//Computes all d2pxy/dt2 once for all:
			_d2pxy.resize(_nbClasses);
			for(unsigned int c = 0; c < _nbClasses; c++) {
				VVdouble * _d2pxy_c = & _d2pxy[c];
				_d2pxy_c -> resize(_nbStates);
				double rc =  _rateDistribution -> getCategory(c);
				RowMatrix<double> d2Q = _model -> getd2Pij_dt2(_brLen * rc);
				for(unsigned int x = 0; x < _nbStates; x++) {
					Vdouble * _d2pxy_c_x = & (* _d2pxy_c)[x];
					_d2pxy_c_x -> resize(_nbStates);
					for(unsigned int y = 0; y < _nbStates; y++) {
						(* _d2pxy_c_x)[y] =  rc * rc * d2Q(x, y);
					}
				}
			}
		}
	}

	computeTreeLikelihood();
	if(_computeDerivatives) {
		computeTreeDLikelihood();	
		computeTreeD2Likelihood();
	}
}

/******************************************************************************/

double TwoTreeLikelihood::getValue() const
throw (Exception)
{
	return - getLogLikelihood();
}

/******************************************************************************/

void TwoTreeLikelihood::initTreeLikelihoods(const SequenceContainer & sequences) throw (Exception)
{
	const Sequence * seq1 = sequences.getSequence(_seqnames[0]);
	const Sequence * seq2 = sequences.getSequence(_seqnames[1]);
	_leafLikelihoods1.resize(_nbDistinctSites);
	_leafLikelihoods2.resize(_nbDistinctSites);
	for(unsigned int i = 0; i < _nbDistinctSites; i++) {
		Vdouble * _leafLikelihoods1_i = & _leafLikelihoods1[i];
		Vdouble * _leafLikelihoods2_i = & _leafLikelihoods2[i];
		_leafLikelihoods1_i -> resize(_nbStates);
		_leafLikelihoods2_i -> resize(_nbStates);
		for(unsigned int s = 0; s < _nbStates; s++) {
			//Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
			//otherwise value set to 0:
			try {
				int state1 = seq1 -> getValue(i);
				int state2 = seq2 -> getValue(i);
				(* _leafLikelihoods1_i)[s] = _model -> getInitValue(s, state1);
				(* _leafLikelihoods2_i)[s] = _model -> getInitValue(s, state2);
			} catch (SequenceNotFoundException & snfe) {
				throw SequenceNotFoundException("TwoTreeLikelihood::initTreelikelihoods. Leaf name in tree not found in site conainer: ", snfe.getSequenceId());
			}
		}
	}

	//Initialize likelihood vector:
	_rootLikelihoods.resize(_nbDistinctSites);
	_rootLikelihoodsS.resize(_nbDistinctSites);
	_rootLikelihoodsSR.resize(_nbDistinctSites);
	for(unsigned int i = 0; i < _nbDistinctSites; i++) {
		VVdouble * _rootLikelihoods_i = & _rootLikelihoods[i];
		Vdouble * _rootLikelihoodsS_i = & _rootLikelihoodsS[i];
		_rootLikelihoods_i -> resize(_nbClasses);
		_rootLikelihoodsS_i -> resize(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			Vdouble * _rootLikelihoods_i_c = & (* _rootLikelihoods_i)[c];
			_rootLikelihoods_i_c -> resize(_nbStates);
			for(unsigned int s = 0; s < _nbStates; s++) {
				(* _rootLikelihoods_i_c)[s] = 1.; //All likelihoods are initialized to 1.
			}
		}
	}

	// Initialize d and d2 likelihoods:
	_dLikelihoods.resize(_nbDistinctSites);
	_d2Likelihoods.resize(_nbDistinctSites);
}

/******************************************************************************/

void TwoTreeLikelihood::computeTreeLikelihood()
{
	for(unsigned int i = 0; i < _nbDistinctSites; i++) {
		VVdouble * _rootLikelihoods_i = & _rootLikelihoods[i];
		Vdouble * _leafLikelihoods1_i = & _leafLikelihoods1[i];
		Vdouble * _leafLikelihoods2_i = & _leafLikelihoods2[i];
		for(unsigned int c = 0; c < _nbClasses; c++) {
			Vdouble * _rootLikelihoods_i_c = & (* _rootLikelihoods_i)[c];
			VVdouble * _pxy_c = & _pxy[c];
			for(unsigned int x = 0; x < _nbStates; x++) {
				Vdouble * _pxy_c_x = & (* _pxy_c)[x];
				double l = 0;
				double l1 = (* _leafLikelihoods1_i)[x];
				for(unsigned int y = 0; y < _nbStates; y++) {
					double l2 = (* _leafLikelihoods2_i)[y];
					l += l1 * l2 * (* _pxy_c_x)[y];
				}
				(* _rootLikelihoods_i_c)[x] = l;
			}
		}
	}
	
	Vdouble f = _model -> getFrequencies();
	Vdouble p = _rateDistribution -> getProbabilities();
	for(unsigned int i = 0; i < _nbDistinctSites; i++) {
		//For each site in the sequence,
		VVdouble * _rootLikelihoods_i = & _rootLikelihoods[i];
		Vdouble * _rootLikelihoodsS_i = & _rootLikelihoodsS[i];
		_rootLikelihoodsSR[i] = 0;
		for(unsigned int c = 0; c < _nbClasses; c++) {
			(* _rootLikelihoodsS_i)[c] = 0;
			//For each rate classe,
			Vdouble * _rootLikelihoods_i_c = & (* _rootLikelihoods_i)[c];
			for(unsigned int x = 0; x < _nbStates; x++) {
				//For each initial state,
				(* _rootLikelihoodsS_i)[c] += f[x] * (* _rootLikelihoods_i_c)[x];
			}
			_rootLikelihoodsSR[i] += p[c] * (* _rootLikelihoodsS_i)[c];
		}
	}
}

/******************************************************************************/

void TwoTreeLikelihood::computeTreeDLikelihood()
{
	for(unsigned int i = 0; i < _nbDistinctSites; i++) {
		Vdouble * _leafLikelihoods1_i = & _leafLikelihoods1[i];
		Vdouble * _leafLikelihoods2_i = & _leafLikelihoods2[i];
		double dli = 0;
		for(unsigned int c = 0; c < _nbClasses; c++) {
			VVdouble * _dpxy_c = & _dpxy[c];
			double dlic = 0;
			for(unsigned int x = 0; x < _nbStates; x++) {
				Vdouble * _dpxy_c_x = & (* _dpxy_c)[x];
				double l1 = (* _leafLikelihoods1_i)[x];
				double dlicx = 0;
				for(unsigned int y = 0; y < _nbStates; y++) {
					double l2 = (* _leafLikelihoods2_i)[y];
					dlicx += l1 * l2 * (* _dpxy_c_x)[y];
				}
				dlic += dlicx * _model -> freq(x);
			}
			dli += dlic * _rateDistribution -> getProbability(c);
		}
		_dLikelihoods[i] = dli / _rootLikelihoodsSR[i];
	}
}

/******************************************************************************/

void TwoTreeLikelihood::computeTreeD2Likelihood()
{
	for(unsigned int i = 0; i < _nbDistinctSites; i++) {
		Vdouble * _leafLikelihoods1_i = & _leafLikelihoods1[i];
		Vdouble * _leafLikelihoods2_i = & _leafLikelihoods2[i];
		double d2li = 0;
		for(unsigned int c = 0; c < _nbClasses; c++) {
			VVdouble * _d2pxy_c = & _d2pxy[c];
			double d2lic = 0;
			for(unsigned int x = 0; x < _nbStates; x++) {
				Vdouble * _d2pxy_c_x = & (* _d2pxy_c)[x];
				double l1 = (* _leafLikelihoods1_i)[x];
				double d2licx = 0;
				for(unsigned int y = 0; y < _nbStates; y++) {
					double l2 = (* _leafLikelihoods2_i)[y];
					d2licx += l1 * l2 * (* _d2pxy_c_x)[y];
				}
				d2lic += d2licx * _model -> freq(x);
			}
			d2li += d2lic * _rateDistribution -> getProbability(c);
		}
		_d2Likelihoods[i] = d2li / _rootLikelihoodsSR[i];
	}
}

/******************************************************************************/

double TwoTreeLikelihood::getFirstOrderDerivative(const string & variable) const
throw (Exception)
{ 
	Parameter * p = _parameters.getParameter(variable);
	if(p == NULL) throw ParameterNotFoundException("TwoTreeLikelihood::getFirstOrderDerivative", variable);
	if(getRateDistributionParameters().getParameter(variable) != NULL) {
		cout << "DEBUB: WARNING!!! Derivatives respective to rate distribution parameter are not implemented." << endl;
		return log(-1.);
	}
	if(getSubstitutionModelParameters().getParameter(variable) != NULL) {
		cout << "DEBUB: WARNING!!! Derivatives respective to substitution model parameters are not implemented." << endl;
		return log(-1.);
	}
	
	//
	// Computation for branch lengths:
	//
	
	// Get the node with the branch whose length must be derivated:
	double d = 0;
	for(unsigned int i = 0; i < _nbDistinctSites; i++) d += _rootWeights[i] * _dLikelihoods[i];
	return -d;
}

/******************************************************************************/

double TwoTreeLikelihood::getSecondOrderDerivative(const string & variable) const
throw (Exception)
{
	Parameter * p = _parameters.getParameter(variable);
	if(p == NULL) throw ParameterNotFoundException("TwoTreeLikelihood::getSecondOrderDerivative", variable);
	if(getRateDistributionParameters().getParameter(variable) != NULL) {
		cout << "DEBUB: WARNING!!! Derivatives respective to rate distribution parameter are not implemented." << endl;
		return log(-1.);
	}
	if(getSubstitutionModelParameters().getParameter(variable) != NULL) {
		cout << "DEBUB: WARNING!!! Derivatives respective to substitution model parameters are not implemented." << endl;
		return log(-1.);
	}
	
	//
	// Computation for branch lengths:
	//
	
	// Get the node with the branch whose length must be derivated:
	double d2 = 0;
	for(unsigned int i = 0; i < _nbDistinctSites; i++) d2 += _rootWeights[i] * (_d2Likelihoods[i] - pow(_dLikelihoods[i], 2));
	return -d2;
}

/******************************************************************************/

DistanceEstimation::DistanceEstimation(SubstitutionModel * model, DiscreteDistribution * rateDist, const SiteContainer * sites, unsigned int verbose, bool computeMat) :
	_model(model),
	_rateDist(rateDist),
	_sites(sites),
	_verbose(verbose)
{
	//_defaultOptimizer = new SimpleMultiDimensions(NULL);
	_defaultOptimizer = new NewtonBrentMetaOptimizer(NULL);
	_defaultOptimizer -> setMessageHandler(NULL);
	_defaultOptimizer -> setProfiler(NULL);
	_optimizer = _defaultOptimizer;
	if(computeMat) computeMatrix();
}

DistanceEstimation::DistanceEstimation(unsigned int verbose) :
	_model(NULL),
	_rateDist(NULL),
	_sites(NULL),
	_verbose(verbose)
{
	//_defaultOptimizer = new SimpleMultiDimensions(NULL);
	_defaultOptimizer = new NewtonBrentMetaOptimizer(NULL);
	_defaultOptimizer -> setMessageHandler(NULL);
	_defaultOptimizer -> setProfiler(NULL);
	_optimizer = _defaultOptimizer;
}

/******************************************************************************/

DistanceMatrix * DistanceEstimation::getMatrix() const { return new DistanceMatrix(* _dist); }

/******************************************************************************/

void DistanceEstimation::computeMatrix() throw (NullPointerException)
{	
	unsigned int n = _sites -> getNumberOfSequences();
	vector<string> names = _sites -> getSequencesNames();
	_dist = new DistanceMatrix(names);
	_optimizer -> getStopCondition() -> setTolerance(0.0001);
	_optimizer -> setVerbose(max(_verbose - 2, 0));
	for(unsigned int i = 0; i < n; i++) {
		(* _dist)(i, i) = 0;
		if(_verbose > 0) { cout << "*"; cout.flush(); }
		for(unsigned int j = i + 1; j < n; j++) {
			if(_verbose > 1) { cout << "."; cout.flush(); }
			//Node * n0 = new Node(0, names[i]);
			//Node * n1 = new Node(1, names[j]);
			//n1 -> setDistanceToFather(0.1);
			//n0 -> addSon(*n1);
			//Tree<Node> * tree = new Tree<Node>(*n0);

			// Likelihood function:
			//DRHomogeneousTreeLikelihood * lik = 
			//	new DRHomogeneousTreeLikelihood(*tree, *_sites, _model, _rateDist, _verbose > 3);
			TwoTreeLikelihood * lik = 
				new TwoTreeLikelihood(names[i], names[j], *_sites, _model, _rateDist, _verbose > 3);
			lik -> setComputeDerivatives(true);
			const Sequence * seqi = _sites -> getSequence(names[i]);
			const Sequence * seqj = _sites -> getSequence(names[j]);
			unsigned int d = SymbolListTools::getNumberOfDistinctPositions(* seqi, * seqj);
			unsigned int g = SymbolListTools::getNumberOfPositionsWithoutGap(* seqi, * seqj);
			lik -> setParameterValue("BrLen", g == 0 ? 0 : (double)d/(double)g);
			// Optimization:
			_optimizer -> setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_AUTO);
			_optimizer -> setFunction(lik);
			ParameterList params = lik -> getBranchLengthsParameters();
			params.addParameters(_parameters);
			_optimizer -> init(params);
			_optimizer -> optimize();
			// Store results:
			//(* _dist)(i, j) = (* _dist)(j, i) = lik -> getParameter("BrLen0");
			(* _dist)(i, j) = (* _dist)(j, i) = lik -> getParameterValue("BrLen");
			delete lik;
			//delete tree;
		}
	}
}

/******************************************************************************/

