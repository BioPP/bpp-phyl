//
// File: OptimizationTools.h
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

#ifndef _OPTIMIZATIONTOOLS_H_
#define _OPTIMIZATIONTOOLS_H_

#include "TreeLikelihood.h"
#include "HomogeneousTreeLikelihood.h"
#include "AbstractHomogeneousTreeLikelihood.h"

// From the STL:
#include <iostream>
using namespace std;

/**
 * @brief Optimization methods for phylogenetic inference.
 *
 * This class provides optimization methods for phylogenetics.
 * Parameters of the optimization methods are set to work with TreeLikelihood
 * object. Some non trivial parameters are left to the user choice (tolerance, maximum
 * number of function evaluation, output streams).
 */
class OptimizationTools
{
	public:
		OptimizationTools();
		virtual ~OptimizationTools();
	
	public:
		
		
		/**
		 * @brief Optimize numerical parameters (branch length, substitution model & rate distribution) of a TreeLikelihood function.
		 *
		 * Uses Newton's method for branch length and Brent's one dimensional method for other parameters
		 * (NewtonBrentMetaOptimizer).
		 *
		 * A condition over function values is used as a stop condition for the algorithm.
		 *
		 * @see NewtonBrentMetaOptimizer
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @param verbose        The verbose level.
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static int optimizeNumericalParameters(
			DiscreteRatesAcrossSitesTreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			unsigned int verbose = 1)
			throw (Exception);
	
		
	private:
		
		class ScaleFunction: public virtual Function {
				
			protected:
				TreeLikelihood * _tl;
				mutable ParameterList _brLen, _lambda;
				
			public:
				ScaleFunction(TreeLikelihood * tl);
				~ScaleFunction();
				
			public:
				void setParameters(const ParameterList & lambda) throw (ParameterNotFoundException, ConstraintException);
				double getValue() const throw (ParameterException);
				ParameterList getParameters() const throw (Exception) { return _lambda; }
				double getParameterValue(const string & name) const throw (ParameterNotFoundException) { return _lambda.getParameter(name) -> getValue(); };
				void setAllParametersValues(const ParameterList & params) 
					throw (ParameterNotFoundException, ConstraintException) {}
				void setParameterValue(const string & name, double value) 
					throw (ParameterNotFoundException, ConstraintException) {}
				void setParametersValues(const ParameterList & params)
					throw (ParameterNotFoundException, ConstraintException) {}
				void matchParametersValues(const ParameterList & params)
					throw (ConstraintException) {};
		};
	
	public:

		/**
		 * @brief Optimize the scale of a TreeLikelihood.
		 *
		 * This method only works on branch lengths parameters.
		 * It multiply all branch length by a factor 'x' which is optimized
		 * using Brent's algorithm in one dimension.
		 * This method may be usefull for scaling a tree whose branch lengths
		 * come from the Neighbor-Joining algorithm for instance.
		 *
		 * Practically, and contrarily to what one may expect, this does not
		 * speed up the optimization!
		 *
		 * A condition over parameters is used as a stop condition for the algorithm.
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @param verbose        The verbose level.
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static int optimizeTreeScale(
			TreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			unsigned int verbose = 1
			)	throw (Exception);
	
};


#endif	//_OPTIMIZATIONTOOLS_H_

