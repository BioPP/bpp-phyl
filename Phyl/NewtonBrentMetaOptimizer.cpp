//
// File: NewtonBrentMetaOptimizer.cpp
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: ue Nov 17 17:22 2004
//

/*
Copyright ou © ou Copr. CNRS, (17 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant à fournir des classes
pour l'analyse de données phylogénétiques.

Ce logiciel est régi par la licence CeCILL soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accepté les
termes.
*/

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

Julien.Dutheil@univ-montp2.fr

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

/**************************************************************************/

#include "NewtonBrentMetaOptimizer.h"

/**************************************************************************/

double NewtonBrentMetaOptimizer::BRANCH_LENGTHS_TOL = 10;
double NewtonBrentMetaOptimizer::RATE_DISTRIBUTION_TOL = 0.01;
double NewtonBrentMetaOptimizer::SUBSTITUTION_MODEL_TOL = 0.01;

/**************************************************************************/

NewtonBrentMetaOptimizer::NewtonBrentMetaOptimizer(DiscreteRatesAcrossSitesTreeLikelihood * tl):
	AbstractOptimizer(tl)
{
	_defaultStopCondition = new FunctionStopCondition(this);
	_stopCondition = _defaultStopCondition;
	_rateDistributionOptimizer = NULL;
	_substitutionModelOptimizer = NULL;
	_branchLengthsOptimizer = NULL;
	_rough=true;
}

/**************************************************************************/

NewtonBrentMetaOptimizer::~NewtonBrentMetaOptimizer()
{
	// Delete all optimizers:
	delete _rateDistributionOptimizer;
	delete _substitutionModelOptimizer;
  delete _branchLengthsOptimizer;
	delete _defaultStopCondition;
}

/**************************************************************************/

void NewtonBrentMetaOptimizer::init(const ParameterList & parameters)
	throw (Exception)
{
	DiscreteRatesAcrossSitesTreeLikelihood * tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(_function);
	_stopCondition -> init();
	_parameters = parameters;
	unsigned int nbParams = _parameters.size();

	// Some cleaning first.
	// This is useful only if the MetaOptimizer have been initialized once before this time.
	delete _rateDistributionOptimizer;
	delete _substitutionModelOptimizer;
	delete _branchLengthsOptimizer;
	
	// Check parameters class:
	ParameterList rateDistParams = tl -> getRateDistributionParameters();
	_nbRateDistParams = 0;
	_rateDistributionParameters.reset();
	for(unsigned int i = 0; i < rateDistParams.size(); i++) {
		Parameter * rateDistParam = rateDistParams[i];
		if(parameters.getParameter(rateDistParam -> getName()) != NULL) {
			_rateDistributionParameters.addParameter(* rateDistParam);
			_nbRateDistParams++;
		}
	}

	ParameterList subsModParams = tl -> getSubstitutionModelParameters();
	_nbSubsModParams = 0;
	_substitutionModelParameters.reset();
	for(unsigned int i = 0; i < subsModParams.size(); i++) {
		Parameter * subsModParam = subsModParams[i];
		if(parameters.getParameter(subsModParam -> getName()) != NULL) {
			_substitutionModelParameters.addParameter(* subsModParam);
			_nbSubsModParams++;
		}
	}

	ParameterList branchLengthParams = tl -> getBranchLengthsParameters();
	_nbBranchLengths = 0;
	_branchLengthsParameters.reset();
	for(unsigned int i = 0; i < branchLengthParams.size(); i++) {
		Parameter * branchLengthParam = branchLengthParams[i];
		if(parameters.getParameter(branchLengthParam -> getName()) != NULL) {
			_branchLengthsParameters.addParameter(* branchLengthParam);
			_nbBranchLengths++;
		}
	}

	
	// Initialize optimizers:
	if(_nbRateDistParams > 0) {
		_rateDistributionOptimizer = new SimpleMultiDimensions(tl);
		_rateDistributionOptimizer -> setProfiler(_profiler);
		_rateDistributionOptimizer -> setMessageHandler(_messageHandler);
		_rateDistributionOptimizer -> setConstraintPolicy(_constraintPolicy);
		_rateDistributionOptimizer -> setVerbose(_verbose > 1 ? 1 : 0);
		//_rateDistributionOptimizer -> setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_IGNORE);
	}
	
	if(_nbSubsModParams > 0) {
		_substitutionModelOptimizer = new SimpleMultiDimensions(tl);
		_substitutionModelOptimizer -> setProfiler(_profiler);
		_substitutionModelOptimizer -> setMessageHandler(_messageHandler);
		_substitutionModelOptimizer -> setConstraintPolicy(_constraintPolicy);
		_substitutionModelOptimizer -> setVerbose(_verbose > 1 ? 1 : 0);
	}

	if(_nbBranchLengths > 0) {
		_branchLengthsOptimizer = new PseudoNewtonOptimizer(tl);
		_branchLengthsOptimizer -> setProfiler(_profiler);
		_branchLengthsOptimizer -> setMessageHandler(_messageHandler);
		_branchLengthsOptimizer -> setConstraintPolicy(_constraintPolicy);
		_branchLengthsOptimizer -> setVerbose(_verbose > 1 ? 1 : 0);
	}
	
	// Dump to profile:
	for(unsigned int i = 0; i < nbParams; i++) {
		profile(_parameters[i] -> getName() + "\t"); 
	}
	profileln("Function");

	printPoint(_parameters, _function -> f(_parameters));

	// Initialize stop condition:
 _stopCondition -> isToleranceReached();

}

/**************************************************************************/

double NewtonBrentMetaOptimizer::optimize()
	throw (Exception)
{
	DiscreteRatesAcrossSitesTreeLikelihood * tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(_function);
	_nbEval = 0;
	
	if(_rough) {
		// First adjust roughly the branch lengths:
		if(_nbBranchLengths > 0) {
			if(_verbose > 0) {
				cout << endl << "Branch lengths (rough):" << endl;
				cout.flush();
			}
			tl -> setComputeDerivatives(true);
			_branchLengthsOptimizer -> getStopCondition() -> setTolerance(BRANCH_LENGTHS_TOL);
			_branchLengthsOptimizer -> init(_branchLengthsParameters);
			_branchLengthsOptimizer -> optimize();
			_nbEval += _branchLengthsOptimizer -> getNumberOfEvaluations();
			tl -> setComputeDerivatives(false);
		}
		
		// Then adjust rate distribution parameters:
		if(_nbRateDistParams > 0) {
			if(_verbose > 0) {
				cout << endl << "Rate distribution (rough):" << endl;
				cout.flush();
			}
			_rateDistributionOptimizer -> getStopCondition() -> setTolerance(RATE_DISTRIBUTION_TOL);
			_rateDistributionOptimizer -> init(_rateDistributionParameters);
			_rateDistributionOptimizer -> optimize();
			_nbEval += _rateDistributionOptimizer -> getNumberOfEvaluations();
		}
		
		// Now adjust substitution parameters:
		if(_nbSubsModParams > 0) {
			if(_verbose > 0) {
				cout << endl << "Substitution model (rough):" << endl;
				cout.flush();
			}
			_substitutionModelOptimizer -> getStopCondition() -> setTolerance(SUBSTITUTION_MODEL_TOL);
			_substitutionModelOptimizer -> init(_substitutionModelParameters);
			_substitutionModelOptimizer -> optimize();
			_nbEval += _substitutionModelOptimizer -> getNumberOfEvaluations();	
		}
		
	}
	
	// Finally adjust all in one until precision:
	double tol = _stopCondition -> getTolerance();
	if(_nbBranchLengths  > 0) _branchLengthsOptimizer     -> getStopCondition() -> setTolerance(tol);
	if(_nbRateDistParams > 0) _rateDistributionOptimizer  -> getStopCondition() -> setTolerance(tol);
	if(_nbSubsModParams  > 0) _substitutionModelOptimizer -> getStopCondition() -> setTolerance(tol);

	// Actualize parameters:
	_parameters.matchParametersValues(tl -> getParameters());
	
	_tolIsReached = false;
	while(_nbEval < _nbEvalMax && !_tolIsReached) {
		
		if(_nbBranchLengths > 0) {
			if(_verbose > 0) {
				cout << endl << "Branch lengths:" << endl;
				cout.flush();
			}
			tl -> setComputeDerivatives(true);
			_branchLengthsParameters.matchParametersValues(_parameters);
			_branchLengthsOptimizer  -> init(_branchLengthsParameters);
			_branchLengthsOptimizer  -> optimize();
			_nbEval += _branchLengthsOptimizer -> getNumberOfEvaluations();
			tl -> setComputeDerivatives(false);
		}

		if(_nbRateDistParams > 0) {
			if(_verbose > 0) {
				cout << endl << "Rate distribution:" << endl;
				cout.flush();
			}
			_rateDistributionParameters.matchParametersValues(_parameters);
			_rateDistributionOptimizer -> init(_rateDistributionParameters);
			_rateDistributionOptimizer -> optimize();
			_nbEval += _rateDistributionOptimizer -> getNumberOfEvaluations();
		}

		if(_nbSubsModParams > 0) {
			if(_verbose > 0) {
				cout << endl << "Substitution model:" << endl;
				cout.flush();
			}
			_substitutionModelParameters.matchParametersValues(_parameters);
			_substitutionModelOptimizer -> init(_substitutionModelParameters);
			_substitutionModelOptimizer -> optimize();
			_nbEval += _substitutionModelOptimizer -> getNumberOfEvaluations();
		}
		
		// Actualize parameters:
		_parameters.matchParametersValues(tl -> getParameters());

		_tolIsReached = _stopCondition -> isToleranceReached();
	}

	return tl -> getValue();
}

/**************************************************************************/

//double NewtonBrentMetaOptimizer::getFunctionValue() const
//{
//	return _tl -> getValue();
//}

/**************************************************************************/

