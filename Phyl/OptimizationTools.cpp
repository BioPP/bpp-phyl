//
// File: OptimizationTools.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Sun Dec 14 09:43:32 2003
//

#include "OptimizationTools.h"

// From NumCalc:
#include <NumCalc/ParameterList.h>
#include <NumCalc/PowellMultiDimensions.h>
#include <NumCalc/DownhillSimplexMethod.h>
#include <NumCalc/BrentOneDimension.h>
#include "GaltierNewtonOptimizer.h"
#include "NewtonBrentMetaOptimizer.h"
#include <NumCalc/OptimizationStopCondition.h>

/******************************************************************************/

OptimizationTools::OptimizationTools() {}

OptimizationTools::~OptimizationTools() {}
	
/******************************************************************************/

int OptimizationTools::optimizeWithDownhillSimplexMethod(
	TreeLikelihood * tl,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler)
	throw (Exception)
{
	// Build optimizer:
	DownhillSimplexMethod * optimizer = new DownhillSimplexMethod(tl);
	optimizer -> setProfiler(profiler);
	optimizer -> setMessageHandler(messageHandler);
	optimizer -> setMaximumNumberOfEvaluations(tlEvalMax);
	ParametersStopCondition * PS = NULL;	
	// Optimize TreeLikelihood function:
	try {
		ParameterList pl = tl -> getParameters();
		optimizer -> setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_AUTO);
		optimizer -> init(pl);
		PS = new ParametersStopCondition(optimizer, tolerance, 30);
		optimizer -> setStopCondition(PS);
		optimizer -> optimize();
	} catch(Exception e) {
		cout << e.what() << endl;
		exit(-1);
	}
	// We're done.
	int n = optimizer -> getNumberOfEvaluations(); 
	// Delete optimizer:
	delete optimizer;
	// Delete StopCondition:
	delete PS;
	// Send number of evaluations done:
	return n;
}

/******************************************************************************/	
	
int OptimizationTools::optimizeWithPowellMethod(
	TreeLikelihood * tl,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler)
	throw (Exception)
{
	// Build optimizer:
	PowellMultiDimensions * optimizer = new PowellMultiDimensions(tl);
	optimizer -> setProfiler(profiler);
	optimizer -> setMessageHandler(messageHandler);
	optimizer -> setMaximumNumberOfEvaluations(tlEvalMax);
	ParametersStopCondition * PS = NULL;
	
	// Optimize TreeLikelihood function:
	try {
		ParameterList pl = tl -> getParameters();
		optimizer -> setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_IGNORE);
		optimizer -> init(pl);
		PS = new ParametersStopCondition(optimizer, tolerance);
		optimizer -> setStopCondition(PS);
		optimizer -> optimize();
	} catch(Exception e) {
		cout << e.what() << endl;
		exit(-1);
	}
	// We're done.
	int n = optimizer -> getNumberOfEvaluations(); 
	// Delete optimizer:
	delete optimizer;
	// Delete StopCondition:
	delete PS;
	// Send number of evaluations done:
	return n;
}
	
/******************************************************************************/

int OptimizationTools::optimizeWithDownhillSimplexAndPowellMethod(
	TreeLikelihood * tl,
	double simplexTolerance,
	double powellTolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler)
	throw (Exception)
{
	int n = optimizeWithDownhillSimplexMethod(
		tl,
		simplexTolerance,
		tlEvalMax,
		messageHandler,
		profiler);
	int n2 = optimizeWithPowellMethod(
		tl,
		powellTolerance,
		tlEvalMax - n,
		messageHandler,
		profiler);
	// We're done.
	return n + n2;
}
	
/******************************************************************************/

OptimizationTools::ScaleFunction::ScaleFunction(TreeLikelihood * tl): _tl(tl) {
	// We work only on the branch lengths:
	_brLen = tl -> getBranchLengthsParameters();
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

int OptimizationTools::optimizeTreeScale(
	TreeLikelihood * tl,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler)
	throw (Exception)
{
	ScaleFunction sf(tl);
	BrentOneDimension bod(&sf);
	bod.setMessageHandler(messageHandler);
	bod.setProfiler(profiler);
	bod.setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_IGNORE);
	bod.setInitialInterval(0.5, 1.5);
	ParameterList singleParameter;
	singleParameter.addParameter(Parameter("scale factor", 0));
	bod.init(singleParameter);
	ParametersStopCondition * PS = new ParametersStopCondition(&bod, tolerance);
	bod.setStopCondition(PS);
	bod.setTolerance(tolerance);
	bod.setMaximumNumberOfEvaluations(tlEvalMax);
	bod.optimize();
	delete PS;
	return bod.getNumberOfEvaluations();
}

/******************************************************************************/

int OptimizationTools::optimizeWithDownhillSimplexMethodAlphaSeparately(
	TreeLikelihood * tl,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler,
	ostream * profilerAlpha)
	throw (Exception)
{
	// Build optimizers:

	DownhillSimplexMethod * optimizer1 = new DownhillSimplexMethod(tl);
	optimizer1 -> setProfiler(profiler);
	optimizer1 -> setMessageHandler(messageHandler);
	optimizer1 -> setMaximumNumberOfEvaluations(tlEvalMax);

	BrentOneDimension * optimizer2 = new BrentOneDimension(tl);
	optimizer2 -> setProfiler(profilerAlpha);
	optimizer2 -> setMessageHandler(messageHandler);
	optimizer2 -> setMaximumNumberOfEvaluations(tlEvalMax);
	optimizer2 -> setTolerance(tolerance);
	optimizer2 -> setInitialInterval(0.1, 0.6);
	
	ParametersStopCondition * PS1 = NULL, * PS2 = NULL;
	
	int nbEval = 0;
	// Optimize TreeLikelihood function:
	try {
		ParameterList pl1 = tl -> getParameters();
		Parameter * p = pl1.getParameter("alpha");
		ParameterList pl2;
		pl2.addParameter(*p);
		pl1.deleteParameter("alpha");
		// Now pl1 contains all parameters without 'alpha'
		// and pl2 contains only parameter 'alpha'.
		
		optimizer1 -> setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_IGNORE);

		optimizer2 -> setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_IGNORE);

		// First we adjust roughtly the branch lengths:
		cout << "Branch lengths: first try..." << endl;
		optimizer1 -> init(pl1);
		PS1 = new ParametersStopCondition(optimizer1, tolerance, 5);
		optimizer1 -> setStopCondition(PS1);
		PS1 -> setTolerance(0.05);
		optimizer1 -> optimize();
		nbEval += optimizer1 -> getNumberOfEvaluations();
		
		// Then we adjust roughtly alpha:
		cout << endl << "Alpha: first try..." << endl;
		optimizer2 -> init(pl2);
		PS2 = new ParametersStopCondition(optimizer2, tolerance, 5);
		optimizer2 -> setStopCondition(PS2);
		PS2 -> setTolerance(0.001);
		optimizer2 -> optimize();
		nbEval += optimizer2 -> getNumberOfEvaluations();
		
		if(nbEval == tlEvalMax) {
			cout << "WARNING!!! 'tlEvalMax' reached before end of full ";
			cout << "optimization in 'optimizeWithDownhillSimplexMethodAlphaSeparately'.";
			cout << endl;
			return nbEval;
		}
		
		// Now we adjust more accurately the branch lengths:
		cout << endl << "Branch lengths: second try..." << endl;
		PS1 -> setTolerance(tolerance);
		PS1 -> resetCounter();
		optimizer1 -> optimize();
		nbEval += optimizer1 -> getNumberOfEvaluations();
		
		if(nbEval == tlEvalMax) {
			cout << "WARNING!!! 'tlEvalMax' reached before end of full ";
			cout << "optimization in 'optimizeWithDownhillSimplexMethodAlphaSeparately'.";
			cout << endl;
			return nbEval;
		}
		
		// Finally we adjust all parameters in the same time:
		bool tolIsReached = false;
		PS2 -> setTolerance(tolerance);
		PS2 -> resetCounter();
		PS1 -> resetCounter();
		PS1 -> setBurnin(pl1.size() + 1);
		cout << endl << "All in one..." << endl;
		while(nbEval < tlEvalMax && !tolIsReached) {
			PS2 -> resetCounter();
			optimizer2 -> optimize();
			nbEval += optimizer2 -> getNumberOfEvaluations();
			optimizer1 -> step();
			tolIsReached = optimizer1 -> isToleranceReached();
			nbEval ++;
		}
	} catch(Exception e) {
		cout << e.what() << endl;
		exit(-1);
	}
	// We're done.

	// Delete optimizer:
	delete optimizer1;
	delete optimizer2;
	// Delete StopCondition:
	delete PS1;
	delete PS2;
	// Send number of evaluations done:
	return nbEval;
}

/******************************************************************************/

int OptimizationTools::optimizeWithPowellMethodAlphaSeparately(
	TreeLikelihood * tl,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler,
	ostream * profilerAlpha)
	throw (Exception)
{
	// Build optimizers:

	PowellMultiDimensions * optimizer1 = new PowellMultiDimensions(tl);
	optimizer1 -> setProfiler(profiler);
	optimizer1 -> setMessageHandler(messageHandler);
	optimizer1 -> setMaximumNumberOfEvaluations(tlEvalMax);

	BrentOneDimension * optimizer2 = new BrentOneDimension(tl);
	optimizer2 -> setProfiler(profilerAlpha);
	optimizer2 -> setMessageHandler(messageHandler);
	optimizer2 -> setMaximumNumberOfEvaluations(tlEvalMax);
	optimizer2 -> setTolerance(tolerance);
	optimizer2 -> setInitialInterval(0.1, 0.6);
	
	ParametersStopCondition * PS1 = NULL, * PS2 = NULL;
	
	int nbEval = 0;
	// Optimize TreeLikelihood function:
	try {
		ParameterList pl1 = tl -> getParameters();
		Parameter * p = pl1.getParameter("alpha");
		ParameterList pl2;
		pl2.addParameter(*p);
		pl1.deleteParameter("alpha");
		// Now pl1 contains all parameters without 'alpha'
		// and pl2 contains only parameter 'alpha'.
		
		optimizer1 -> setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_IGNORE);

		optimizer2 -> setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_IGNORE);

		// First we adjust roughtly the branch lengths:
		cout << "Branch lengths: first try..." << endl;
		optimizer1 -> init(pl1);
		PS1 = new ParametersStopCondition(optimizer1, tolerance, 5);
		optimizer1 -> setStopCondition(PS1);
		PS1 -> setTolerance(0.05);
		optimizer1 -> optimize();
		nbEval += optimizer1 -> getNumberOfEvaluations();
		
		// Then we adjust roughtly alpha:
		cout << endl << "Alpha: first try..." << endl;
		optimizer2 -> init(pl2);
		PS2 = new ParametersStopCondition(optimizer2, tolerance, 5);
		optimizer2 -> setStopCondition(PS2);
		PS2 -> setTolerance(0.001);
		optimizer2 -> optimize();
		nbEval += optimizer2 -> getNumberOfEvaluations();
		
		if(nbEval == tlEvalMax) {
			cout << "WARNING!!! 'tlEvalMax' reached before end of full ";
			cout << "optimization in 'optimizeWithPowellMethodAlphaSeparately'.";
			cout << endl;
			return nbEval;
		}
		
		// Now we adjust more accurately the branch lengths:
		cout << endl << "Branch lengths: second try..." << endl;
		PS1 -> setTolerance(tolerance);
		PS1 -> resetCounter();
		optimizer1 -> optimize();
		nbEval += optimizer1 -> getNumberOfEvaluations();
		
		if(nbEval == tlEvalMax) {
			cout << "WARNING!!! 'tlEvalMax' reached before end of full ";
			cout << "optimization in 'optimizeWithPowellMethodAlphaSeparately'.";
			cout << endl;
			return nbEval;
		}
		
		// Finally we adjust all parameters in the same time:
		bool tolIsReached = false;
		PS2 -> setTolerance(tolerance);
		PS2 -> resetCounter();
		PS1 -> resetCounter();
		PS1 -> setBurnin(pl1.size() + 1);
		cout << endl << "All in one..." << endl;
		while(nbEval < tlEvalMax && !tolIsReached) {
			PS2 -> resetCounter();
			optimizer2 -> optimize();
			nbEval += optimizer2 -> getNumberOfEvaluations();
			optimizer1 -> step();
			tolIsReached = optimizer1 -> isToleranceReached();
			nbEval ++;
		}
	} catch(Exception e) {
		cout << e.what() << endl;
		exit(-1);
	}
	// We're done.

	// Delete optimizer:
	delete optimizer1;
	delete optimizer2;
	// Delete StopCondition:
	delete PS1;
	delete PS2;
	// Send number of evaluations done:
	return nbEval;
}

/******************************************************************************/
	
int OptimizationTools::optimizeWithNewtonMethod(
	TreeLikelihood * tl,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler)
	throw (Exception)
{
	// Build optimizer:
	GaltierNewtonOptimizer * optimizer = new GaltierNewtonOptimizer(tl);
	optimizer -> setProfiler(profiler);
	optimizer -> setMessageHandler(messageHandler);
	optimizer -> setMaximumNumberOfEvaluations(tlEvalMax);
	optimizer -> getStopCondition() -> setTolerance(tolerance);
	
	// Optimize TreeLikelihood function:
	try {
		ParameterList pl = tl -> getParameters();
		optimizer -> setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_AUTO);
		optimizer -> init(pl);
		optimizer -> optimize();
	} catch(Exception e) {
		cout << e.what() << endl;
		exit(-1);
	}
	// We're done.
	int n = optimizer -> getNumberOfEvaluations(); 
	// Delete optimizer:
	delete optimizer;
	// Send number of evaluations done:
	return n;
}
	
/******************************************************************************/	

int OptimizationTools::optimizeWithNewtonBrentMethod(
	AbstractHomogeneousTreeLikelihood * tl,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler)
	throw (Exception)
{
	// Build optimizer:
	NewtonBrentMetaOptimizer * optimizer = new NewtonBrentMetaOptimizer(tl);
	optimizer -> setProfiler(profiler);
	optimizer -> setMessageHandler(messageHandler);
	optimizer -> setMaximumNumberOfEvaluations(tlEvalMax);
	optimizer -> getStopCondition() -> setTolerance(tolerance);
	
	// Optimize TreeLikelihood function:
	try {
		ParameterList pl = tl -> getParameters();
		optimizer -> setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_AUTO);
		optimizer -> init(pl);
		optimizer -> optimize();
	} catch(Exception e) {
		cout << e.what() << endl;
		exit(-1);
	}
	// We're done.
	int n = optimizer -> getNumberOfEvaluations(); 
	// Delete optimizer:
	delete optimizer;
	// Send number of evaluations done:
	return n;
}
	
/******************************************************************************/	

