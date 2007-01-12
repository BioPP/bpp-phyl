//
// File: NewtonBrentMetaOptimizer.cpp
// Created by: Julien Dutheil
// Created on: ue Nov 17 17:22 2004
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

NewtonBrentMetaOptimizer::NewtonBrentMetaOptimizer(DiscreteRatesAcrossSitesTreeLikelihood * tl, unsigned int n):
  AbstractOptimizer(tl)
{
  _defaultStopCondition = new FunctionStopCondition(this);
  _stopCondition = _defaultStopCondition->clone();
  _newtonOptimizer = NULL;
  _brentOptimizer = NULL;
  _n = n;
  _precisionStep = log10(_stopCondition->getTolerance()) / _n;
  _stepCount = 0;
}

/**************************************************************************/

NewtonBrentMetaOptimizer::NewtonBrentMetaOptimizer(
    const NewtonBrentMetaOptimizer & opt):
  AbstractOptimizer(opt)
{
  _newtonParameters   = opt._newtonParameters;
  _brentParameters    = opt._brentParameters;
  _newtonOptimizer    = opt._newtonOptimizer->clone();
  _brentOptimizer     = opt._brentOptimizer->clone();
  _nbNewtonParameters = opt._nbNewtonParameters;
  _nbBrentParameters  = opt._nbBrentParameters;
  _n                  = opt._n;
  _precisionStep      = opt._precisionStep;
  _stepCount          = opt._stepCount;
}

/**************************************************************************/

NewtonBrentMetaOptimizer & NewtonBrentMetaOptimizer::operator=(
    const NewtonBrentMetaOptimizer & opt)
{
  AbstractOptimizer::operator=(opt);
  _newtonParameters   = opt._newtonParameters;
  _brentParameters    = opt._brentParameters;
  _newtonOptimizer    = opt._newtonOptimizer->clone();
  _brentOptimizer     = opt._brentOptimizer->clone();
  _nbNewtonParameters = opt._nbNewtonParameters;
  _nbBrentParameters  = opt._nbBrentParameters;
  _n                  = opt._n;
  _precisionStep      = opt._precisionStep;
  _stepCount          = opt._stepCount;
  return *this;
}

/**************************************************************************/

NewtonBrentMetaOptimizer::~NewtonBrentMetaOptimizer()
{
  // Delete all optimizers:
  delete _newtonOptimizer;
  delete _brentOptimizer;
}

/**************************************************************************/

void NewtonBrentMetaOptimizer::init(const ParameterList & parameters)
  throw (Exception)
{
  DiscreteRatesAcrossSitesTreeLikelihood * tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(_function);
  _stopCondition->init();
  _parameters = parameters;
  unsigned int nbParams = _parameters.size();
  
  // Some cleaning first.
  // This is useful only if the MetaOptimizer have been initialized once before this time.
  delete _newtonOptimizer;
  delete _brentOptimizer;
  
  // Check parameters class:
  ParameterList branchLengthParams = tl->getBranchLengthsParameters();
  _nbNewtonParameters = 0;
  _newtonParameters.reset();
  for(unsigned int i = 0; i < branchLengthParams.size(); i++)
  {
    Parameter * branchLengthParam = branchLengthParams[i];
    if(parameters.getParameter(branchLengthParam->getName()) != NULL)
    {
      _newtonParameters.addParameter(* branchLengthParam);
      _nbNewtonParameters++;
    }
  }

  _nbBrentParameters = 0;
  _brentParameters.reset();
  ParameterList rateDistParams = tl->getRateDistributionParameters();
  for(unsigned int i = 0; i < rateDistParams.size(); i++)
  {
    Parameter * rateDistParam = rateDistParams[i];
    if(parameters.getParameter(rateDistParam->getName()) != NULL)
    {
      _brentParameters.addParameter(* rateDistParam);
      _nbBrentParameters++;
    }
  }
  ParameterList subsModParams = tl->getSubstitutionModelParameters();
  for(unsigned int i = 0; i < subsModParams.size(); i++)
  {
    Parameter * subsModParam = subsModParams[i];
    if(parameters.getParameter(subsModParam->getName()) != NULL)
    {
      _brentParameters.addParameter(* subsModParam);
      _nbBrentParameters++;
    }
  }
  
  // Initialize optimizers:
  if(_nbNewtonParameters > 0)
  {
    _newtonOptimizer = new PseudoNewtonOptimizer(tl);
    _newtonOptimizer->setProfiler(_profiler);
    _newtonOptimizer->setMessageHandler(_messageHandler);
    _newtonOptimizer->setConstraintPolicy(_constraintPolicy);
    _newtonOptimizer->setVerbose(_verbose > 0 ? _verbose - 1 : 0);
  }
  if(_nbBrentParameters > 0)
  {
    _brentOptimizer = new SimpleMultiDimensions(tl);
    _brentOptimizer->setProfiler(_profiler);
    _brentOptimizer->setMessageHandler(_messageHandler);
    _brentOptimizer->setConstraintPolicy(_constraintPolicy);
    _brentOptimizer->setVerbose(_verbose > 0 ? _verbose - 1 : 0);
  }
  
  // Dump to profile:
  for(unsigned int i = 0; i < nbParams; i++)
  {
    profile(_parameters[i]->getName() + "\t"); 
  }
  profileln("Function");

  printPoint(_parameters, _function->f(_parameters));
  
  // Actualize parameters:
  _parameters.matchParametersValues(tl->getParameters());
  
  _tolIsReached = false;
  
  // Reset counter:
  _stepCount = 0;
  // Recompute step if precision has changed:
  _precisionStep = log10(_stopCondition->getTolerance()) / _n;
}

/**************************************************************************/

double NewtonBrentMetaOptimizer::step() throw (Exception)
{
  DiscreteRatesAcrossSitesTreeLikelihood * tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(_function);
  _stepCount++;
  
  if(_verbose == 0) { cout << "*"; cout.flush(); }

  double tol = _stopCondition->getTolerance();
  if(_stepCount < _n)
  {
    tol = pow(10, _stepCount*_precisionStep);
  }
  
  if(_nbNewtonParameters > 0)
  {
    if(_verbose - 1 > 0)
    {
      cout << endl << "Branch lengths: ";
      cout.flush();
    }
    tl->setComputeDerivatives(true);
    _newtonParameters.matchParametersValues(_parameters);
    _newtonOptimizer->getStopCondition()->setTolerance(tol);
    _newtonOptimizer->init(_newtonParameters);
    _newtonOptimizer->optimize();
    _nbEval += _newtonOptimizer->getNumberOfEvaluations();
     tl->setComputeDerivatives(false);
     if(_verbose - 1 > 0) cout << endl;
  }

  if(_nbBrentParameters > 0)
  {
    _brentParameters.matchParametersValues(_parameters);
    _brentOptimizer->getStopCondition()->setTolerance(tol);
    _brentOptimizer->init(_brentParameters);
    _brentOptimizer->step();
    _nbEval += _brentOptimizer->getNumberOfEvaluations();
     if(_verbose - 1 > 0) cout << endl;
  }
    
  // Actualize parameters:
  _parameters.matchParametersValues(tl->getParameters());

  _tolIsReached = _stopCondition->isToleranceReached();
  return tl->getValue();
}

/**************************************************************************/

double NewtonBrentMetaOptimizer::optimize()
  throw (Exception)
{
  DiscreteRatesAcrossSitesTreeLikelihood * tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(_function);
  _nbEval = 0;
  
  while(_nbEval < _nbEvalMax && !_tolIsReached)
  {  
    step();
  }

  return tl->getValue();
}

/**************************************************************************/

