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
#include "ClockMetaOptimizer.h"
#include "NNISearchable.h"
#include "NNITopologySearch.h"

// From Utils:
#include <Utils/ApplicationTools.h>

// From NumCalc:
#include <NumCalc/ParameterList.h>
#include <NumCalc/PowellMultiDimensions.h>
#include <NumCalc/SimpleNewtonMultiDimensions.h>
#include <NumCalc/DownhillSimplexMethod.h>
#include <NumCalc/BrentOneDimension.h>
#include <NumCalc/OptimizationStopCondition.h>
#include <NumCalc/FivePointsNumericalDerivative.h>
#include <NumCalc/ThreePointsNumericalDerivative.h>

/******************************************************************************/

OptimizationTools::OptimizationTools() {}

OptimizationTools::~OptimizationTools() {}
	
/******************************************************************************/

OptimizationTools::ScaleFunction::ScaleFunction(TreeLikelihood * tl): _tl(tl)
{
	// We work only on the branch lengths:
	_brLen = tl->getBranchLengthsParameters();
  _lambda.addParameter(Parameter("scale factor", 1, &Parameter::R_PLUS_STAR)); 
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
	for(unsigned int i = 0; i < brLen.size(); i++)
  {
		brLen[i]->setValue(brLen[i]->getValue() * _lambda[0]->getValue());
	}
	return _tl->f(brLen);
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
	bod.setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
	bod.setInitialInterval(0.99, 1.01);
	ParameterList singleParameter;
	singleParameter.addParameter(Parameter("scale factor", 1.));
	ParametersStopCondition PS(&bod, tolerance);
	bod.setStopCondition(PS);
	bod.setMaximumNumberOfEvaluations(tlEvalMax);
	bod.init(singleParameter);
	bod.optimize();
	return bod.getNumberOfEvaluations();
}

/******************************************************************************/

int OptimizationTools::optimizeNumericalParameters(
	DiscreteRatesAcrossSitesTreeLikelihood * tl,
  OptimizationListener * listener,
  unsigned int nstep,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler,
	unsigned int verbose
	)	throw (Exception)
{
	// Build optimizer:
	NewtonBrentMetaOptimizer * optimizer = new NewtonBrentMetaOptimizer(tl, NewtonBrentMetaOptimizer::TYPE_PSEUDONEWTON, NewtonBrentMetaOptimizer::IT_TYPE_FULL, nstep);
	optimizer->setVerbose(verbose);
	optimizer->setProfiler(profiler);
	optimizer->setMessageHandler(messageHandler);
	optimizer->setMaximumNumberOfEvaluations(tlEvalMax);
	optimizer->getStopCondition()->setTolerance(tolerance);
	
	// Optimize TreeLikelihood function:
	try
  {
    optimizer->setDerivableParameters(tl->getDerivableParameters().getParameterNames());
    optimizer->setNonDerivableParameters(tl->getNonDerivableParameters().getParameterNames());
    ParameterList pl = tl->getParameters();
		optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
		optimizer->init(pl);
    if(listener) optimizer->addOptimizationListener(listener);
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

int OptimizationTools::optimizeNumericalParameters2(
	DiscreteRatesAcrossSitesTreeLikelihood * tl,
  OptimizationListener * listener,
  string method,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler,
	unsigned int verbose
	)	throw (Exception)
{
  int n = 0;
  
  // Build optimizer:
	PseudoNewtonOptimizer * optimizer = new PseudoNewtonOptimizer(tl);
	optimizer->setVerbose(verbose);
	optimizer->setProfiler(profiler);
	optimizer->setMessageHandler(messageHandler);
	optimizer->setMaximumNumberOfEvaluations(tlEvalMax);
	optimizer->getStopCondition()->setTolerance(tolerance);
	
	// Optimize TreeLikelihood function:
	try
  {
		ParameterList pl = tl->getParameters();
		optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
		optimizer->init(pl);
    if(listener) optimizer->addOptimizationListener(listener);
		optimizer->optimize();
	}
  catch(Exception e)
  {
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
	optimizer->setVerbose(verbose);
	optimizer->setProfiler(profiler);
	optimizer->setMessageHandler(messageHandler);
	optimizer->setMaximumNumberOfEvaluations(tlEvalMax);
	optimizer->getStopCondition()->setTolerance(tolerance);
	
	// Optimize TreeLikelihood function:
	try
  {
		ParameterList pl = tl->getBranchLengthsParameters();
		optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
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

int OptimizationTools::optimizeNumericalParametersWithGlobalClock(
	ClockTreeLikelihood * cl,
  unsigned int nstep,
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
    fun = new ThreePointsNumericalDerivative(cl);
  }
  else if(method == "5points")
  {
    fun = new FivePointsNumericalDerivative(cl);
  }
  else throw Exception("Unknow numerical derivative method: " + method + ".");

  //Numerical derivatives:
  vector<string> variables = cl->getTotalHeightParameter().getParameterNames();
  fun->setParametersToDerivate(variables);

  ClockMetaOptimizer * optimizer = new ClockMetaOptimizer(cl, nstep);
	optimizer->setVerbose(verbose);
	optimizer->setProfiler(profiler);
	optimizer->setMessageHandler(messageHandler);
	optimizer->setMaximumNumberOfEvaluations(tlEvalMax);
	optimizer->getStopCondition()->setTolerance(tolerance*100);
	// Optimize TreeLikelihood function:
	try
  {
    ParameterList pl = cl->getParameters();
		optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
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

int OptimizationTools::optimizeNumericalParametersWithGlobalClock2(
	ClockTreeLikelihood * cl,
  unsigned int nstep,
  string method,
	double tolerance,
	int tlEvalMax,
	ostream * messageHandler,
	ostream * profiler,
	unsigned int verbose
	)	throw (Exception)
{
  cl->enableFirstOrderDerivatives(true);
  cl->enableSecondOrderDerivatives(true);
  AbstractNumericalDerivative * fun = NULL;
  if(method == "3points")
  {
    fun = new ThreePointsNumericalDerivative(cl);
  }
  else if(method == "5points")
  {
    fun = new FivePointsNumericalDerivative(cl);
  }
  else throw Exception("Unknow numerical derivative method: " + method + ".");

  //Numerical derivatives:
  vector<string> variables = cl->getTotalHeightParameter().getParameterNames();
  fun->setParametersToDerivate(variables);

	// Build optimizer:
  ClockOptimizationListener * listener = new ClockOptimizationListener(cl);
	NewtonBrentMetaOptimizer * optimizer = new NewtonBrentMetaOptimizer(fun, NewtonBrentMetaOptimizer::TYPE_PSEUDONEWTON, NewtonBrentMetaOptimizer::IT_TYPE_FULL, nstep);
	optimizer->updateParameters(true);
  optimizer->setVerbose(verbose);
	optimizer->setProfiler(profiler);
	optimizer->setMessageHandler(messageHandler);
	optimizer->setMaximumNumberOfEvaluations(tlEvalMax);
	optimizer->getStopCondition()->setTolerance(tolerance);
	// Optimize TreeLikelihood function:
	try
  {
    ParameterList dpl = cl->getBranchLengthsParameters();
    dpl.addParameters(cl->getTotalHeightParameter());
    optimizer->setDerivableParameters(dpl.getParameterNames());
    
    ParameterList ndpl = cl->getSubstitutionModelParameters();
    ndpl.addParameters(cl->getRateDistributionParameters());
    optimizer->setNonDerivableParameters(ndpl.getParameterNames());
		
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
		ParameterList pl = cl->getParameters();
		optimizer->init(pl);
	  optimizer->getNewtonOptimizer()->addOptimizationListener(listener);
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
  delete listener;
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
    OptimizationTools::optimizeNumericalParameters(likelihood, NULL, 1, _tolerance, 1000000, _messenger, _profiler, _verbose);
    _optimizeCounter = 0;
  }
}

/******************************************************************************/	

NNIHomogeneousTreeLikelihood * OptimizationTools::optimizeTreeNNI(
    NNIHomogeneousTreeLikelihood * tl,
    bool optimizeNumFirst,
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
  if(optimizeNumFirst)
    OptimizationTools::optimizeNumericalParameters(tl, NULL, 1, tolBefore, 1000000, messageHandler, profiler, verbose);
  //Begin topo search:
  NNITopologySearch topoSearch(*tl, NNITopologySearch::PHYML, verbose > 2 ? verbose - 2 : 0);
  NNITopologyListener *topoListener = new NNITopologyListener(&topoSearch, tolDuring, messageHandler, profiler, verbose);
  topoListener->setNumericalOptimizationCounter(numStep);
  topoSearch.addTopologyListener(*topoListener);
  topoSearch.search();
  delete topoListener;
  return dynamic_cast<NNIHomogeneousTreeLikelihood *>(topoSearch.getSearchableObject());
}

/******************************************************************************/	

DRTreeParsimonyScore * OptimizationTools::optimizeTreeNNI(
        DRTreeParsimonyScore * tp,
        unsigned int verbose)	
{
  NNISearchable *topo = dynamic_cast<NNISearchable *>(tp);
  NNITopologySearch topoSearch(*topo, NNITopologySearch::PHYML, verbose);
  topoSearch.search();
  return dynamic_cast<DRTreeParsimonyScore *>(topoSearch.getSearchableObject());
};

/******************************************************************************/	

string OptimizationTools::DISTANCEMETHOD_INIT       = "init";
string OptimizationTools::DISTANCEMETHOD_PAIRWISE   = "pairwise";
string OptimizationTools::DISTANCEMETHOD_ITERATIONS = "iterations";

/******************************************************************************/	

TreeTemplate<Node> * OptimizationTools::buildDistanceTree(
    DistanceEstimation & estimationMethod,
    AgglomerativeDistanceMethod & reconstructionMethod,
    const ParameterList & parametersToIgnore,
    bool optimizeBrLen,
    bool rooted,
    const string & param,
    double tolerance,
    unsigned int tlEvalMax,
    ostream * profiler,
    ostream * messenger,
    unsigned int verbose) throw (Exception)
{
  estimationMethod.resetAdditionalParameters();
  estimationMethod.setVerbose(verbose);
  if(param == DISTANCEMETHOD_PAIRWISE)
  {
    ParameterList tmp = estimationMethod.getModel()->getParameters();
    tmp.addParameters(estimationMethod.getRateDistribution()->getParameters());
    tmp.deleteParameters(parametersToIgnore.getParameterNames());
    estimationMethod.setAdditionalParameters(tmp);
  }
  TreeTemplate<Node> * tree = NULL;
  bool test = true;
  double previousLL = -log(0.);
  double currentLL = -log(0.);
  while(test)
  {
    //Compute matrice:
    if(verbose > 0)
      ApplicationTools::displayTask("Estimating distance matrix", true); 
    estimationMethod.computeMatrix();
    DistanceMatrix * matrix = estimationMethod.getMatrix();
    if(verbose > 0)
      ApplicationTools::displayTaskDone();

    //Compute tree:
    if(verbose > 0)
      ApplicationTools::displayTask("Building tree");
    reconstructionMethod.setDistanceMatrix(*matrix);
    reconstructionMethod.computeTree(rooted);
    if(tree) delete tree;
    delete matrix;
    tree = dynamic_cast<TreeTemplate<Node> *>(reconstructionMethod.getTree());
    if(verbose > 0)
      ApplicationTools::displayTaskDone();
    if(param != DISTANCEMETHOD_ITERATIONS) break; //Ends here.
    
    //Now, re-estimate parameters:
    DRHomogeneousTreeLikelihood tl(*tree, *estimationMethod.getData(), estimationMethod.getModel(), estimationMethod.getRateDistribution(), true, verbose > 2);
    tl.initialize();
    if(!optimizeBrLen)
    {
      vector<string> vs = tl.getBranchLengthsParameters().getParameterNames();
      for(unsigned int i = 0; i < vs.size(); i++)
      {
        tl.ignoreParameter(vs[i]);
      }
    }
    for(unsigned int i = 0; i < parametersToIgnore.size(); i++)
      tl.ignoreParameter(parametersToIgnore[i]->getName());
    optimizeNumericalParameters(&tl, NULL, 0, tolerance, tlEvalMax, messenger, profiler, verbose - 1);
    previousLL = currentLL;
    currentLL = tl.getLogLikelihood();
    test = (std::abs(currentLL - previousLL) > tolerance);
  }
  return tree;
}

/******************************************************************************/	

