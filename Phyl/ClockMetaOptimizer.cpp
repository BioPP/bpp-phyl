//
// File: ClockMetaOptimizer.cpp
// Created by: Julien Dutheil
// Created on: Wed May 02 16:33 2007
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

#include "ClockMetaOptimizer.h"

/**************************************************************************/

ClockMetaOptimizer::ClockMetaOptimizer(ClockTreeLikelihood * function, unsigned int n):
  AbstractOptimizer(function), _newtonListener(function), _numDer(function)
{
  _defaultStopCondition = new FunctionStopCondition(this);
  _stopCondition = _defaultStopCondition->clone();
  _newton1Optimizer = NULL;
  _newton2Optimizer = NULL;
  _brentOptimizer = NULL;
  _nbNewton2Parameters = 0;
  _n = n;
  _precisionStep = log10(_stopCondition->getTolerance()) / _n;
  _stepCount = 0;
  _newton1Parameters = function->getDerivableParameters();
  _newton1Parameters.addParameters(function->getTotalHeightParameter());
  _nbNewton1Parameters = _newton1Parameters.size();
  _brentParameters = function->getSubstitutionModelParameters();
  _brentParameters.addParameters(function->getRateDistributionParameters());

  _nbBrentParameters = _brentParameters.size();
  _stepChar = "";
}

/**************************************************************************/

ClockMetaOptimizer::ClockMetaOptimizer(
    const ClockMetaOptimizer & opt):
  AbstractOptimizer(opt), _newtonListener(opt._newtonListener), _numDer(opt._numDer)
{
  _newton1Parameters   = opt._newton1Parameters;
  _newton2Parameters   = opt._newton2Parameters;
  _brentParameters     = opt._brentParameters;
  _newton1Optimizer    = opt._newton1Optimizer->clone();
  _newton2Optimizer    = opt._newton2Optimizer->clone();
  _brentOptimizer      = opt._brentOptimizer->clone();
  _nbNewton1Parameters = opt._nbNewton1Parameters;
  _nbNewton2Parameters = opt._nbNewton2Parameters;
  _nbBrentParameters   = opt._nbBrentParameters;
  _n                   = opt._n;
  _precisionStep       = opt._precisionStep;
  _stepCount           = opt._stepCount;
}

/**************************************************************************/

ClockMetaOptimizer & ClockMetaOptimizer::operator=(
    const ClockMetaOptimizer & opt)
{
  AbstractOptimizer::operator=(opt);
  _newtonListener      = opt._newtonListener;
  _newton1Parameters   = opt._newton1Parameters;
  _newton2Parameters   = opt._newton2Parameters;
  _brentParameters     = opt._brentParameters;
  _newton1Optimizer    = opt._newton1Optimizer->clone();
  _newton2Optimizer    = opt._newton2Optimizer->clone();
  _brentOptimizer      = opt._brentOptimizer->clone();
  _nbNewton1Parameters = opt._nbNewton1Parameters;
  _nbNewton2Parameters = opt._nbNewton2Parameters;
  _nbBrentParameters   = opt._nbBrentParameters;
  _n                   = opt._n;
  _precisionStep       = opt._precisionStep;
  _stepCount           = opt._stepCount;
  _numDer              = opt._numDer;
  return *this;
}

/**************************************************************************/

ClockMetaOptimizer::~ClockMetaOptimizer()
{
  // Delete all optimizers:
  if(_newton1Optimizer) delete _newton1Optimizer;
  if(_newton2Optimizer) delete _newton2Optimizer;
  if(_brentOptimizer)   delete _brentOptimizer;
}

/**************************************************************************/

void ClockMetaOptimizer::doInit(const ParameterList & parameters)
  throw (Exception)
{
  // Some cleaning first.
  // This is useful only if the MetaOptimizer have been initialized once before this time.
  if(_newton1Optimizer) delete _newton1Optimizer;
  if(_newton2Optimizer) delete _newton2Optimizer;
  if(_brentOptimizer)   delete _brentOptimizer;
  
  // Initialize optimizers:
  if(_nbNewton1Parameters == 0 && _nbBrentParameters == 0)
    throw Exception("ClockMetaOptimizer::init(). No derivable and non-derivable parameter set.");
  if(_nbNewton1Parameters > 0)
  {
    _newton1Optimizer = new PseudoNewtonOptimizer(&_numDer);
    //Numerical derivatives:
    if(_function->getParameters().getParameter("TotalHeight") != NULL)
    {
      vector<string> variables = dynamic_cast<ClockTreeLikelihood *>(_function)->getTotalHeightParameter().getParameterNames();
      _numDer.setParametersToDerivate(variables);
    }
    _newton1Optimizer->setProfiler(_profiler);
    _newton1Optimizer->setMessageHandler(_messageHandler);
    _newton1Optimizer->setConstraintPolicy(_constraintPolicy);
    _newton1Optimizer->setVerbose(_verbose > 0 ? _verbose - 1 : 0);
    _newton1Optimizer->addOptimizationListener(&_newtonListener);
    
    _newton2Optimizer = new SimpleNewtonMultiDimensions(dynamic_cast<DerivableSecondOrder *>(_function));
    dynamic_cast<SimpleNewtonMultiDimensions *>(_newton2Optimizer)->updateParameters(true);
    _newton2Optimizer->setProfiler(_profiler);
    _newton2Optimizer->setMessageHandler(_messageHandler);
    _newton2Optimizer->setConstraintPolicy(_constraintPolicy);
    _newton2Optimizer->setVerbose(_verbose > 0 ? _verbose - 1 : 0);
  }
  if(_nbBrentParameters > 0)
  {
    _brentOptimizer = new SimpleMultiDimensions(_function);
    _brentOptimizer->setProfiler(_profiler);
    _brentOptimizer->setMessageHandler(_messageHandler);
    _brentOptimizer->setConstraintPolicy(_constraintPolicy);
    _brentOptimizer->setVerbose(_verbose > 0 ? _verbose - 1 : 0);
  }
  
  // Actualize parameters:
  _parameters.matchParametersValues(_function->getParameters());
  
  // Reset counter:
  _stepCount = 0;
  // Recompute step if precision has changed:
  _precisionStep = log10(_stopCondition->getTolerance()) / _n;
}

/**************************************************************************/

double ClockMetaOptimizer::doStep() throw (Exception)
{
   _stepCount++;
  
  double tol = _stopCondition->getTolerance();
  if(_stepCount < _n)
  {
    tol = pow(10, _stepCount * _precisionStep);
  }
  
  if(_nbNewton1Parameters > 0)
  {
    if(_verbose > 1)
    {
      cout << endl << "Derivable parameters, method 1:" << endl;
    }
    dynamic_cast<DerivableSecondOrder *>(_function)->enableSecondOrderDerivatives(true);
    _newton1Parameters.matchParametersValues(_parameters);
    _newtonListener.resetParameters();
    _newton1Optimizer->setVerbose(max((int)_verbose - 1, 0));
    _newton1Optimizer->getStopCondition()->setTolerance(tol);
    _newton1Optimizer->init(_newton1Parameters);
    _newton1Optimizer->optimize();
    _nbEval += _newton1Optimizer->getNumberOfEvaluations();
    
    //Check for conflicting parameters:
    ParameterList pl = _newtonListener.getParameters();
    for(unsigned int i = pl.size(); i > 0; i--)
      if(_newton1Parameters.getParameter(pl[i-1]->getName()) == NULL)
        pl.deleteParameter(i-1);
    for(unsigned int i = 0; i < pl.size(); i++)
      _newton1Parameters.deleteParameter(pl[i]->getName());
    _newton2Parameters.addParameters(pl);
    _nbNewton1Parameters -= pl.size();
    _nbNewton2Parameters += pl.size();
    
    dynamic_cast<DerivableSecondOrder *>(_function)->enableSecondOrderDerivatives(false);
    _parameters.matchParametersValues(_function->getParameters());
    if(_verbose > 1) cout << endl;
  }

  if(_nbNewton2Parameters > 0)
  {
    if(_verbose > 1)
    {
      cout << endl << "Derivable parameters, method 2:" << endl;
    }
    dynamic_cast<ClockTreeLikelihood *>(_function)->enableSecondOrderDerivatives(true);
    dynamic_cast<ClockTreeLikelihood *>(_function)->resetHeightsConstraints();
    _newton2Parameters.matchParametersValues(_parameters);
    _newton2Optimizer->setVerbose(max((int)_verbose - 1, 0));
    _newton2Optimizer->getStopCondition()->setTolerance(tol);
    _newton2Optimizer->init(_newton2Parameters);
    _newton2Optimizer->step();
    _nbEval += _newton2Optimizer->getNumberOfEvaluations();

    //for(unsigned int i = _newton2Parameters.size(); i > 0; i--)
    //{
    //  if(std::abs(_function->getParameterValue(_newton2Parameters[i-1]->getName()) - _newton2Parameters[i-1]->getValue()) < 0.000000001)
    //  {
    //    _newton1Parameters.addParameter(*_newton2Parameters[i-1]);
    //    _newton2Parameters.deleteParameter(i-1);
    //  }
    //}
    
    dynamic_cast<DerivableSecondOrder *>(_function)->enableSecondOrderDerivatives(false);
    _parameters.matchParametersValues(_function->getParameters());
    if(_verbose > 1) cout << endl;
  }

  if(_nbBrentParameters > 0)
  {
    if(_verbose > 1)
    {
      cout << endl << "Non-derivable parameters:" << endl;
    }
    _brentParameters.matchParametersValues(_parameters);
    _brentOptimizer->setVerbose(max((int)_verbose - 1, 0));
    _brentOptimizer->getStopCondition()->setTolerance(tol);
    _brentOptimizer->init(_brentParameters);
    _brentOptimizer->step();
    _nbEval += _brentOptimizer->getNumberOfEvaluations();
    _parameters.matchParametersValues(_function->getParameters());
     if(_verbose > 1) cout << endl;
  }
    
  return _function->getValue();
}

/**************************************************************************/

