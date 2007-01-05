//
// File: OptimizationTools.cpp
// Created by: Julien Dutheil
// Created on: Sun Dec 14 09:43:32 2003
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

#include "OptimizationTools.h"
#include "PseudoNewtonOptimizer.h"
#include "NewtonBrentMetaOptimizer.h"
#include "NNISearchable.h"
#include "NNITopologySearch.h"

// From Utils:
#include <Utils/ApplicationTools.h>

// From NumCalc:
#include <NumCalc/ParameterList.h>
#include <NumCalc/PowellMultiDimensions.h>
#include <NumCalc/DownhillSimplexMethod.h>
#include <NumCalc/BrentOneDimension.h>
#include <NumCalc/OptimizationStopCondition.h>
#include <NumCalc/FivePointsNumericalDerivative.h>
#include <NumCalc/ThreePointsNumericalDerivative.h>

/******************************************************************************/

OptimizationTools::OptimizationTools() {}

OptimizationTools::~OptimizationTools() {}
	
/******************************************************************************/

OptimizationTools::ScaleFunction::ScaleFunction(TreeLikelihood * tl): _tl(tl) {
	// We work only on the branch lengths:
	_brLen = tl -> getBranchLengthsParameters();
  _lambda.addParameter(Parameter("scale factor", 1)); 
}
	
OptimizationTools::ScaleFunction::~ScaleFunction() {}

void OptimizationTools::ScaleFunction::setParameters(const ParameterList & lambda)
throw (ParameterNotFoundException, ConstraintException)
{
	if(lambda.size() != 1) throw Exception("OptimizationTools::ScaleFunction::f(). This is a one parameter function!");
	_lambda.setParametersValues(lambda);
}

double OptimizationTools::ScaleFunction::getValue() const
throw (ParameterException)
{
	// Scale the tree:
	ParameterList brLen = _brLen;
	for(unsigned int i = 0; i < brLen.size(); i++) {
		brLen[i] -> setValue(brLen[i] -> getValue() * _lambda[0] -> getValue());
	}
	return _tl -> f(brLen);
}

/******************************************************************************/

int OptimizationTools::optimizeTreeScale(
	TreeLikelihood * tl,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler,
	unsigned int verbose
	)	throw (Exception)
{
	ScaleFunction sf(tl);
	BrentOneDimension bod(&sf);
	bod.setMessageHandler(messageHandler);
	bod.setProfiler(profiler);
	bod.setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_IGNORE);
	bod.setInitialInterval(0.5, 1.5);
	ParameterList singleParameter;
	singleParameter.addParameter(Parameter("scale factor", 1.));
	bod.init(singleParameter);
	ParametersStopCondition PS(&bod, tolerance);
	bod.setStopCondition(PS);
	bod.setMaximumNumberOfEvaluations(tlEvalMax);
	bod.optimize();
	return bod.getNumberOfEvaluations();
}

/******************************************************************************/

int OptimizationTools::optimizeNumericalParameters(
	DiscreteRatesAcrossSitesTreeLikelihood * tl,
  unsigned int nstep,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler,
	unsigned int verbose
	)	throw (Exception)
{
	// Build optimizer:
	NewtonBrentMetaOptimizer * optimizer = new NewtonBrentMetaOptimizer(tl, nstep);
	optimizer->setVerbose(verbose);
	optimizer->setProfiler(profiler);
	optimizer->setMessageHandler(messageHandler);
	optimizer->setMaximumNumberOfEvaluations(tlEvalMax);
	optimizer->getStopCondition()->setTolerance(tolerance);
	
	// Optimize TreeLikelihood function:
	try {
		ParameterList pl = tl->getParameters();
		optimizer->setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_AUTO);
		optimizer->init(pl);
		optimizer->optimize();
	} catch(Exception e) {
		cout << e.what() << endl;
		exit(-1);
	}
	// We're done.
	int n = optimizer->getNumberOfEvaluations(); 
	// Delete optimizer:
	delete optimizer;
	// Send number of evaluations done:
	return n;
}
	
/******************************************************************************/

int OptimizationTools::optimizeNumericalParameters2(
	DiscreteRatesAcrossSitesTreeLikelihood * tl,
  string method,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler,
	unsigned int verbose
	)	throw (Exception)
{
  AbstractNumericalDerivative * fun = NULL;
  if(method == "3points")
  {
    fun = new ThreePointsNumericalDerivative(tl);
  }
  else if(method == "5points")
  {
    fun = new FivePointsNumericalDerivative(tl);
  }
  else throw Exception("Unknow numerical derivative method: " + method + ".");

  int n = 0;
  
  //Numerical derivatives:
  vector<string> variables;
  ParameterList subs = tl->getSubstitutionModelParameters();
  ParameterList dist = tl->getRateDistributionParameters();
  for(unsigned int i = 0; i < subs.size(); i++)
  {
    variables.push_back(subs[i]->getName());
  }
  for(unsigned int i = 0; i < dist.size(); i++)
  {
    variables.push_back(dist[i]->getName());
  }
  fun->setParametersToDerivate(variables);
  
  // Build optimizer:
	PseudoNewtonOptimizer * optimizer = new PseudoNewtonOptimizer(fun);
	optimizer->setVerbose(verbose);
	optimizer->setProfiler(profiler);
	optimizer->setMessageHandler(messageHandler);
	optimizer->setMaximumNumberOfEvaluations(tlEvalMax);
	optimizer->getStopCondition()->setTolerance(tolerance);
	
	// Optimize TreeLikelihood function:
	try {
		ParameterList pl = tl->getParameters();
		optimizer->setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_AUTO);
		optimizer->init(pl);
		optimizer->optimize();
	} catch(Exception e) {
    cout << "Error in optimization process." << endl;
		cout << e.what() << endl;
		exit(-1);
	}
  //Set function at the best point:
  tl->f(optimizer->getParameters());
	// We're done.
	n += optimizer->getNumberOfEvaluations(); 
	// Delete optimizer:
	delete optimizer;
	// Send number of evaluations done:
	return n;
}

/******************************************************************************/

int OptimizationTools::optimizeBranchLengthsParameters(
	DiscreteRatesAcrossSitesTreeLikelihood * tl,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler,
	unsigned int verbose
	)	throw (Exception)
{
	// Build optimizer:
	PseudoNewtonOptimizer * optimizer = new PseudoNewtonOptimizer(tl);
	optimizer -> setVerbose(verbose);
	optimizer -> setProfiler(profiler);
	optimizer -> setMessageHandler(messageHandler);
	optimizer -> setMaximumNumberOfEvaluations(tlEvalMax);
	optimizer -> getStopCondition() -> setTolerance(tolerance);
	
	// Optimize TreeLikelihood function:
	try
  {
		ParameterList pl = tl->getBranchLengthsParameters();
		optimizer->setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_AUTO);
		optimizer->init(pl);
		optimizer->optimize();
	}
  catch(Exception e)
  {
		cout << e.what() << endl;
		exit(-1);
	}
	// We're done.
	int n = optimizer->getNumberOfEvaluations(); 
	// Delete optimizer:
	delete optimizer;
	// Send number of evaluations done:
	return n;
}
	
/******************************************************************************/	

void NNITopologyListener::topologyChangeTested(const TopologyChangeEvent & event)
{
}

/******************************************************************************/	

void NNITopologyListener::topologyChangeSuccessful(const TopologyChangeEvent & event)
{
  _optimizeCounter++;
  if(_optimizeCounter == _optimizeNumerical)
  {
    DiscreteRatesAcrossSitesTreeLikelihood * likelihood = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(_topoSearch->getSearchableObject());
    OptimizationTools::optimizeNumericalParameters(likelihood, 1, _tolerance, 1000000, _messenger, _profiler, _verbose);
    _optimizeCounter = 0;
  }
}

/******************************************************************************/	

DiscreteRatesAcrossSitesTreeLikelihood * OptimizationTools::optimizeTreeNNI(
    DiscreteRatesAcrossSitesTreeLikelihood * tl,
		double tolBefore,
		double tolDuring,
		int tlEvalMax,
    unsigned int numStep,
		ostream * messageHandler,
		ostream * profiler,
		unsigned int verbose)
  throw (Exception)
{
  //Roughly optimize parameter
  OptimizationTools::optimizeNumericalParameters(tl, 1, tolBefore, 1000000, messageHandler, profiler, verbose);
  //Begin topo search:
  NNISearchable *topo = dynamic_cast<NNISearchable *>(tl);
  NNITopologySearch topoSearch(*topo, NNITopologySearch::PHYML, verbose);
  NNITopologyListener *topoListener = new NNITopologyListener(&topoSearch, tolDuring, messageHandler, profiler, verbose);
  topoListener->setNumericalOptimizationCounter(numStep);
  topoSearch.addTopologyListener(*topoListener);
  topoSearch.search();
  delete topoListener;
  return dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(topoSearch.getSearchableObject());
}

/******************************************************************************/	

