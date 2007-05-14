//
// File: NewtonBrentMetaOptimizer.cpp
// Created by: Julien Dutheil
// Created on: Tue Nov 17 17:22 2004
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

string NewtonBrentMetaOptimizer::TYPE_SIMPLENEWTON= "simple";
string NewtonBrentMetaOptimizer::TYPE_PSEUDONEWTON= "pseudo";
string NewtonBrentMetaOptimizer::IT_TYPE_STEP= "step";
string NewtonBrentMetaOptimizer::IT_TYPE_FULL= "full";

/**************************************************************************/

NewtonBrentMetaOptimizer::NewtonBrentMetaOptimizer(DerivableSecondOrder * function, const string & type, const string & itType, unsigned int n):
  AbstractOptimizer(function)
{
  _defaultStopCondition = new FunctionStopCondition(this);
  _stopCondition = _defaultStopCondition->clone();
  _newtonOptimizer = NULL;
  _brentOptimizer = NULL;
  _nbNewtonParameters = 0;
  _nbBrentParameters = 0;
  _n = n;
  _precisionStep = log10(_stopCondition->getTolerance()) / _n;
  _stepCount = 0;
  _type = type;
  _itType = itType;
  _stepChar = "";
}

/**************************************************************************/

NewtonBrentMetaOptimizer::NewtonBrentMetaOptimizer(
    const NewtonBrentMetaOptimizer & opt):
  AbstractOptimizer(opt)
{
  _newtonParameters   = opt._newtonParameters;
  _brentParameters    = opt._brentParameters;
  _newtonOptimizer    = opt._newtonOptimizer ? opt._newtonOptimizer->clone() : NULL;
  _brentOptimizer     = opt._brentOptimizer? opt._brentOptimizer->clone() : NULL;
  _nbNewtonParameters = opt._nbNewtonParameters;
  _nbBrentParameters  = opt._nbBrentParameters;
  _n                  = opt._n;
  _precisionStep      = opt._precisionStep;
  _stepCount          = opt._stepCount;
  _type               = opt._type;
  _itType             = opt._itType;
}

/**************************************************************************/

NewtonBrentMetaOptimizer & NewtonBrentMetaOptimizer::operator=(
    const NewtonBrentMetaOptimizer & opt)
{
  AbstractOptimizer::operator=(opt);
  _newtonParameters   = opt._newtonParameters;
  _brentParameters    = opt._brentParameters;
  _newtonOptimizer    = opt._newtonOptimizer ? opt._newtonOptimizer->clone() : NULL;
  _brentOptimizer     = opt._brentOptimizer? opt._brentOptimizer->clone() : NULL;
  _nbNewtonParameters = opt._nbNewtonParameters;
  _nbBrentParameters  = opt._nbBrentParameters;
  _n                  = opt._n;
  _precisionStep      = opt._precisionStep;
  _stepCount          = opt._stepCount;
  _type               = opt._type;
  _itType             = opt._itType;
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

void NewtonBrentMetaOptimizer::doInit(const ParameterList & parameters)
  throw (Exception)
{
  // Some cleaning first.
  // This is useful only if the MetaOptimizer have been initialized once before this time.
  if(_newtonOptimizer) delete _newtonOptimizer;
  if(_brentOptimizer)  delete _brentOptimizer;

  _newtonParameters.reset();
  for(unsigned int i = 0; i < _newtonParametersNames.size(); i++)
  {
    const Parameter *p = parameters.getParameter(_newtonParametersNames[i]);
    if(p) _newtonParameters.addParameter(*p);
  }
  _nbNewtonParameters = _newtonParameters.size();

  _brentParameters.reset();
  for(unsigned int i = 0; i < _brentParametersNames.size(); i++)
  {
    const Parameter *p = parameters.getParameter(_brentParametersNames[i]);
    if(p) _brentParameters.addParameter(*p);
  }
  _nbBrentParameters = _brentParameters.size();
  
  // Initialize optimizers:
  if(_nbNewtonParameters == 0 && _nbBrentParameters == 0)
    throw Exception("NewtonBrentMetaOptimizer::init(). No derivable and non-derivable parameter set.");
  if(_nbNewtonParameters > 0)
  {
    if(_type == TYPE_PSEUDONEWTON) 
      _newtonOptimizer = new PseudoNewtonOptimizer(dynamic_cast<DerivableSecondOrder *>(_function));
    else if(_type == TYPE_SIMPLENEWTON)
      _newtonOptimizer = new SimpleNewtonMultiDimensions(dynamic_cast<DerivableSecondOrder *>(_function));
    else throw Exception("NewtonBrentMetaOptimizer::init. Unknown type of newton optimizer.");

    dynamic_cast<AbstractOptimizer *>(_newtonOptimizer)->updateParameters(_updateParameters);
    _newtonOptimizer->setProfiler(_profiler);
    _newtonOptimizer->setMessageHandler(_messageHandler);
    _newtonOptimizer->setConstraintPolicy(_constraintPolicy);
    _newtonOptimizer->setVerbose(_verbose > 0 ? _verbose - 1 : 0);
  }
  if(_nbBrentParameters > 0)
  {
    _brentOptimizer = new SimpleMultiDimensions(_function);
    dynamic_cast<AbstractOptimizer *>(_brentOptimizer)->updateParameters(_updateParameters);
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
  _function->setParameters(_parameters);
}

/**************************************************************************/

double NewtonBrentMetaOptimizer::doStep() throw (Exception)
{
   _stepCount++;
  
  double tol = _stopCondition->getTolerance();
  if(_stepCount < _n)
  {
    tol = pow(10, _stepCount * _precisionStep);
  }
  
  if(_nbNewtonParameters > 0)
  {
    if(_verbose > 1)
    {
      cout << endl << "Derivable parameters:";
      cout.flush();
    }
    dynamic_cast<DerivableSecondOrder *>(_function)->enableSecondOrderDerivatives(true);
    _newtonParameters.matchParametersValues(_parameters);
    _newtonOptimizer->getStopCondition()->setTolerance(tol);
    _newtonOptimizer->init(_newtonParameters);
    if(_itType == IT_TYPE_STEP)
      _newtonOptimizer->step();
    else if(_itType == IT_TYPE_FULL)
      _newtonOptimizer->optimize();
    else throw Exception("NewtonBrentMetaOptimizer::step. Unknown iteration type specified.");
    _nbEval += _newtonOptimizer->getNumberOfEvaluations();
    dynamic_cast<DerivableSecondOrder *>(_function)->enableSecondOrderDerivatives(false);
    if(_verbose > 1) cout << endl;
  }

  if(_nbBrentParameters > 0)
  {
    if(_verbose > 1)
    {
      cout << endl << "Non-derivable parameters:" << endl;
    }
    _brentParameters.matchParametersValues(_parameters);
    _brentOptimizer->getStopCondition()->setTolerance(tol);
    _brentOptimizer->init(_brentParameters);
    _brentOptimizer->step();
    _nbEval += _brentOptimizer->getNumberOfEvaluations();
     if(_verbose > 1) cout << endl;
  }
    
  // Actualize parameters:
  _parameters.matchParametersValues(_function->getParameters());
  return _function->getValue();
}

/**************************************************************************/

