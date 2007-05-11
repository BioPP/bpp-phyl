//
// File: ClockMetaOptimizer.h
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

#ifndef _CLOCKMETAOPTIMIZER_H_
#define _CLOCKMETAOPTIMIZER_H_

#include "AbstractHomogeneousTreeLikelihood.h"
#include "PseudoNewtonOptimizer.h"
#include "ClockTreeLikelihood.h"

// From NumCalc:
#include <NumCalc/SimpleMultiDimensions.h>
#include <NumCalc/SimpleNewtonMultiDimensions.h>
#include <NumCalc/AbstractOptimizer.h>
#include <NumCalc/ThreePointsNumericalDerivative.h>

// From the STL:
#include <vector>
using namespace std;

/**
 */
class ClockMetaOptimizer:
  public AbstractOptimizer
{
  protected:
    class Listener:
      public OptimizationListener
      {
        protected:
          ClockTreeLikelihood * _tl;
          ParameterList _parameters;

        public:
          Listener(ClockTreeLikelihood * tl): _tl(tl) {}

        public:
          void optimizationInitializationPerformed(const OptimizationEvent & event)
          {
            _tl->updateHeightsConstraints();
          }
          void optimizationStepPerformed(const OptimizationEvent & event)
          {
            ParameterList tmp = _tl->getConflictingParameters();
            for(unsigned int i = 0; i < tmp.size(); i++)
            {
              if(_parameters.getParameter(tmp[i]->getName()) == NULL)
              {
                _parameters.addParameter(*tmp[i]);
              }
            }
            _tl->updateHeightsConstraints();
          }
          bool listenerModifiesParameters() const { return false; }
          ParameterList getParameters() const { return _parameters; }
          void resetParameters() { _parameters.reset(); }
      };

	protected:
		ParameterList _newton1Parameters;
		ParameterList _newton2Parameters;
		ParameterList _brentParameters;
		Optimizer * _newton1Optimizer;
		Optimizer * _newton2Optimizer;
		Optimizer * _brentOptimizer;
		unsigned int _nbNewton1Parameters;
		unsigned int _nbNewton2Parameters;
		unsigned int _nbBrentParameters;
    unsigned int _n;
    double _precisionStep;
    unsigned int _stepCount;
    string _type;
    string _itType;
    Listener _newtonListener;
    ThreePointsNumericalDerivative _numDer;
		
	public:
    /**
     * @brief Build a new NewtonBrentMetaOptimizer object.
     *
     * @param tl     A likelihood function.
     * @param n      The number of progressive steps to use in optimization (only used with IT_TYPE_FULL method).
     */
		ClockMetaOptimizer(ClockTreeLikelihood * function, unsigned int n = 1);

		virtual ~ClockMetaOptimizer();

    ClockMetaOptimizer(const ClockMetaOptimizer & opt);

    ClockMetaOptimizer & operator=(const ClockMetaOptimizer & opt);

    ClockMetaOptimizer * clone() const { return new ClockMetaOptimizer(*this); }

	public:
		
		void doInit(const ParameterList & parameters) throw (Exception);
    
    double doStep() throw (Exception);

    Optimizer * getNewton1Optimizer() { return _newton1Optimizer; }
    const Optimizer * getNewton1Optimizer() const { return _newton1Optimizer; }
    Optimizer * getNewton2Optimizer() { return _newton2Optimizer; }
    const Optimizer * getNewton2Optimizer() const { return _newton2Optimizer; }
    Optimizer * getBrentOptimizer() { return _brentOptimizer; }
    const Optimizer * getBrentOptimizer() const { return _brentOptimizer; }

};

#endif //_CLOCKMETAOPTIMIZER_H_

