//
// File: OptimizationTools.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Sun Dec 14 09:43:32 2003
//

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
 * NB: For now, only numerical parameter optimization is performed.
 *
 * This class provides optimization methods for phylogenetics.
 * They are combinations of classical optimization methods that can be found
 * in the NumCalc library.
 *
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
		 * @brief Optimize a TreeLikelihood object with the Downhill Simplex method.
		 *
		 * This method optimize all parameters with the Downhill Simple method.
		 * The constraint mode is set to 'AUTO', meaning that all parameters are
		 * replaced by parameters of class AutoParameter.
		 * Indeed, this method does not work on likelihood objects with values 
		 * out of constraint: some null value for branch length for instance lead
		 * to NaN values for likelihood,  which 'attract' the algorithm.
		 *
		 * A condition over parameters is used as a stop condition for the algorithm.
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static int optimizeWithDownhillSimplexMethod(
			TreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout)
			throw (Exception);

		/**
		 * @brief Optimize a TreeLikelihood object with Powell's multidimensions method.
		 *
		 * This method optimize all parameters with Powell's multidimensions method.
		 * The constraint mode is set to 'IGNORE', meaning that all constraints on
		 * parameters are removed.
		 * This is because Powell's algorithm may lead to transitionnal value out of
		 * parameter constraints. AutoParameter objects hence may lead to values that
		 * does not match the local - i.e. one-dimensional - optimum.
		 *
		 * A condition over parameters is used as a stop condition for the algorithm.
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static int optimizeWithPowellMethod(
			TreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout)
			throw (Exception);
		
		/**
		 * @brief Optimize a TreeLikelihood object with Newton's method.
		 *
		 * This method optimize all parameters with Galtier's modified Newton's method.
		 * The constraint mode is set to 'AUTO', meaning that all constraints on
		 * parameters are taken into acount (see the AutoParameter class).
		 *
		 * A condition over function values is used as a stop condition for the algorithm.
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static int optimizeWithNewtonMethod(
			TreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout)
			throw (Exception);
	
		/**
		 * @brief Optimize a TreeLikelihood object with Newton's method for branch length
		 * and Brent's one dimensional method for other parameters.
		 *
		 * A condition over function values is used as a stop condition for the algorithm.
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static int optimizeWithNewtonBrentMethod(
			AbstractHomogeneousTreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout)
			throw (Exception);
	
		static int optimizeWithDownhillSimplexAndPowellMethod(
			TreeLikelihood * tl,
			double ratio,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout)
			throw (Exception);
		
	private:
		
		class ScaleFunction: public Function {
				
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
				double getParameter(const string & name) const throw (ParameterNotFoundException) { return _lambda.getParameter(name) -> getValue(); };
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
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static int optimizeTreeScale(
			TreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout)
			throw (Exception);
	
		static int optimizeWithDownhillSimplexMethodAlphaSeparately(
			TreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			ostream * profilerAlpha  = &cout)
			throw (Exception);	

		static int optimizeWithPowellMethodAlphaSeparately(
			TreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			ostream * profilerAlpha  = &cout)
			throw (Exception);	
};


#endif	//_OPTIMIZATIONTOOLS_H_
