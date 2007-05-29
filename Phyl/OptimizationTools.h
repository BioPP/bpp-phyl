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
#include "DRHomogeneousTreeLikelihood.h"
#include "NNIHomogeneousTreeLikelihood.h"
#include "ClockTreeLikelihood.h"
#include "NNITopologySearch.h"
#include "DRTreeParsimonyScore.h"
#include "TreeTemplate.h"
#include "DistanceEstimation.h"
#include "AgglomerativeDistanceMethod.h"

// From NumCalc:
#include <NumCalc/SimpleNewtonMultiDimensions.h>

// From the STL:
#include <iostream>
using namespace std;

/**
 * @brief Listener used internally by the optimizeTreeNNI method.
 */
class NNITopologyListener:
  public TopologyListener
{
  protected:
    NNITopologySearch * _topoSearch;
    double _tolerance;
    ostream *_messenger;
    ostream *_profiler;
    unsigned int _verbose;
    unsigned int _optimizeCounter;
    unsigned int _optimizeNumerical;

  public:
    NNITopologyListener(NNITopologySearch * ts, double tolerance, ostream *messenger, ostream *profiler, unsigned int verbose):
      _topoSearch(ts), _tolerance(tolerance),
      _messenger(messenger), _profiler(profiler),
      _verbose(verbose),
      _optimizeCounter(0), _optimizeNumerical(1) {}
    virtual ~NNITopologyListener() {}

  public:
    void topologyChangeTested(const TopologyChangeEvent & event);
    void topologyChangeSuccessful(const TopologyChangeEvent & event);
    void setNumericalOptimizationCounter(unsigned int c) { _optimizeNumerical = c; }

};






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

    static string OPTIMIZATION_GRADIENT;
    static string OPTIMIZATION_NEWTON;
    static string OPTIMIZATION_2POINTS;
    static string OPTIMIZATION_3POINTS;
    static string OPTIMIZATION_5POINTS;
		
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
     * @param listener       A pointer toward an optimization listener, if needed.
     * @param nstep          The number of progressive steps to perform (see NewtonBrentMetaOptimizer). 1 means full precision from start.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @param verbose        The verbose level.
     * @param optMethod      Optimization type for derivable parameters (first or second order derivatives).
     * @see OPTIMIZATION_NEWTON, OPTIMIZATION_GRADIENT
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static unsigned int optimizeNumericalParameters(
			DiscreteRatesAcrossSitesTreeLikelihood * tl,
      OptimizationListener * listener = NULL,
      unsigned int nstep = 1,
			double tolerance = 0.000001,
			unsigned int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			unsigned int verbose = 1,
      const string & method = OPTIMIZATION_NEWTON)
			throw (Exception);
	
		/**
		 * @brief Optimize numerical parameters (branch length, substitution model & rate distribution) of a TreeLikelihood function.
		 *
		 * Uses Newton's method for all parameters, branch length derivatives are computed analytically, derivatives for other parameters numerically.
		 *
		 * @see PseudoNewtonOptimizer
		 * @see ThreePointsNumericalDerivative
		 * @see FivePointsNumericalDerivative
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
     * @param listener       A pointer toward an optimization listener, if needed.
     * @param derMethod      Numerical derivative computation method. Must be one of "3points" or "5points", otherwise an exception is thrown.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @param verbose        The verbose level.
     * @param optMethod      Optimization type for derivable parameters (first or second order derivatives).
     * @see OPTIMIZATION_3POINTS, OPTIMIZATION_3POINTS, OPTIMIZATION_5POINTS, OPTIMIZATION_NEWTON, OPTIMIZATION_GRADIENT
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static unsigned int optimizeNumericalParameters2(
			DiscreteRatesAcrossSitesTreeLikelihood * tl,
      OptimizationListener * listener = NULL,
      const string & derMethod = OPTIMIZATION_3POINTS,
			double tolerance = 0.000001,
			unsigned int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			unsigned int verbose = 1,
      const string & optMethod = OPTIMIZATION_NEWTON)
			throw (Exception);

    /**
		 * @brief Optimize branch lengths parameters of a TreeLikelihood function.
		 *
		 * Uses Newton's method.
		 *
		 * A condition over function values is used as a stop condition for the algorithm.
		 *
		 * @see NewtonBrentMetaOptimizer
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
     * @param listener       A pointer toward an optimization listener, if needed.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @param verbose        The verbose level.
     * @param optMethod      Optimization type for derivable parameters (first or second order derivatives).
     * @see OPTIMIZATION_NEWTON, OPTIMIZATION_GRADIENT
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static unsigned int optimizeBranchLengthsParameters(
			DiscreteRatesAcrossSitesTreeLikelihood * tl,
      OptimizationListener * listener = NULL,
			double tolerance = 0.000001,
			unsigned int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			unsigned int verbose = 1,
      const string & optMethod = OPTIMIZATION_NEWTON)
			throw (Exception);
		
		/**
		 * @brief Optimize numerical parameters assuming a global clock (branch heights, substitution model & rate distribution) of a ClockTreeLikelihood function.
		 *
		 * Uses Newton or conjugate gradient method for branch length and Brent's one dimensional method for other parameters
		 * (NewtonBrentMetaOptimizer).
     * Derivatives are computed analytically.
		 *
		 * A condition over function values is used as a stop condition for the algorithm.
		 *
		 * @see NewtonBrentMetaOptimizer
		 * @see TwoPointsNumericalDerivative
		 * @see ThreePointsNumericalDerivative
		 * @see FivePointsNumericalDerivative
		 *
		 * @param cl             A pointer toward the ClockTreeLikelihood object to optimize.
     * @param listener       A pointer toward an optimization listener, if needed.
     * @param nstep          The number of progressive steps to perform (see NewtonBrentMetaOptimizer). 1 means full precision from start.
     * @param derMethod      Numerical derivative computation method. Must be one of "2points", "3points" or "5points", otherwise an exception is thrown.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @param verbose        The verbose level.
     * @param optMethod      Optimization type for derivable parameters (first or second order derivatives).
     * @see OPTIMIZATION_2POINTS, OPTIMIZATION_3POINTS, OPTIMIZATION_5POINTS, OPTIMIZATION_NEWTON, OPTIMIZATION_GRADIENT
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static unsigned int optimizeNumericalParametersWithGlobalClock(
			ClockTreeLikelihood * cl,
      OptimizationListener * listener = NULL,
      unsigned int nstep = 1,
      const string & method = OPTIMIZATION_3POINTS,
			double tolerance = 0.000001,
			unsigned int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			unsigned int verbose = 1,
      const string & optMethod = OPTIMIZATION_GRADIENT)
			throw (Exception);

    /**
		 * @brief Optimize numerical parameters assuming a global clock (branch heights, substitution model & rate distribution) of a ClockTreeLikelihood function.
		 *
		 * Uses Newton or conjugate gradient method for all parameters, branch length derivatives are computed analytically, derivatives for other parameters numerically.
		 *
		 * @see PseudoNewtonOptimizer
		 * @see TwoPointsNumericalDerivative
		 * @see ThreePointsNumericalDerivative
		 * @see FivePointsNumericalDerivative
		 *
		 * @param cl             A pointer toward the ClockTreeLikelihood object to optimize.
     * @param derMethod      Numerical derivative computation method. Must be one of "2points", "3points" or "5points", otherwise an exception is thrown.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @param verbose        The verbose level.
     * @param optMethod      Optimization type for derivable parameters (first or second order derivatives).
     * @see OPTIMIZATION_2POINTS, OPTIMIZATION_3POINTS, OPTIMIZATION_5POINTS, OPTIMIZATION_NEWTON, OPTIMIZATION_GRADIENT
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static unsigned int optimizeNumericalParametersWithGlobalClock2(
			ClockTreeLikelihood * cl,
      OptimizationListener * listener = NULL,
      const string & method = OPTIMIZATION_3POINTS,
			double tolerance = 0.000001,
			unsigned int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			unsigned int verbose = 1,
      const string & optMethod = OPTIMIZATION_GRADIENT
      )
			throw (Exception);


	private:
		
		class ScaleFunction:
      public Function
    {
				
			protected:
				TreeLikelihood * _tl;
				mutable ParameterList _brLen, _lambda;
				
			public:
				ScaleFunction(TreeLikelihood * tl);
				virtual ~ScaleFunction();

        ScaleFunction * clone() const { return new ScaleFunction(*this); }
				
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
		 * @throw Exception any exception thrown by the optimizer.
		 */
		static unsigned int optimizeTreeScale(
			  TreeLikelihood * tl,
			  double tolerance = 0.000001,
			  int tlEvalMax = 1000000,
			  ostream * messageHandler = &cout,
			  ostream * profiler       = &cout,
			  unsigned int verbose = 1)
      throw (Exception);

    /**
     * @brief Optimize all parameters from a TreeLikelihood object, including tree topology using Nearest Neighbor Interchanges.
     *
     * This function takes as input a TreeLikelihood object implementing the NNISearchable interface.
     *
     * Details:
     * A NNITopologySearch object is instanciated and is associated an additional TopologyListener.
     * This listener is used to re-estimate numerical parameters after one or several topology change.
     * By default, the PHYML option is used for the NNITopologySearch object, and numerical parameters are re-estimated
     * every 4 NNI runs (as in the phyml software).
     *
     * The optimizeNumericalParameters method is used for estimating numerical parameters.
     * The tolerance passed to this function is specified as input parameters.
     * They are generally very high to avoid local optima.
     *
		 * @param tl               A pointer toward the TreeLikelihood object to optimize.
     * @param optimizeNumFirst Tell if we must optimize numerical parameters before searching topology.
		 * @param tolBefore        The tolerance to use when estimating numerical parameters before topology search (if optimizeNumFirst is set to 'true').
		 * @param tolDuring        The tolerance to use when estimating numerical parameters during the topology search.
		 * @param tlEvalMax        The maximum number of function evaluations.
     * @param numStep          Number of NNI rounds before re-estimating numerical parameters.
		 * @param messageHandler   The massage handler.
		 * @param profiler         The profiler.
		 * @param verbose          The verbose level.
     * @return A pointer toward the final likelihood object.
     * This pointer may be the same as passed in argument (tl), but in some cases the algorithm
     * clone this object. We may change this bahavior in the future...
     * You hence should write something like
     * @code
     * tl = OptimizationTools::optimizeTreeNNI(tl, ...);
     * @endcode
		 * @throw Exception any exception thrown by the optimizer.
     */
    static NNIHomogeneousTreeLikelihood * optimizeTreeNNI(
        NNIHomogeneousTreeLikelihood * tl,
        bool optimizeNumFirst = true,
			  double tolBefore = 100,
			  double tolDuring = 100,
			  int tlEvalMax = 1000000,
        unsigned int numStep = 1,
			  ostream * messageHandler = &cout,
			  ostream * profiler       = &cout,
			  unsigned int verbose = 1)
      throw (Exception);

    static DRTreeParsimonyScore * optimizeTreeNNI(
        DRTreeParsimonyScore * tp,
        unsigned int verbose = 1);

    /**
     * @brief Build a tree using a distance method.
     *
     * This method estimate a distance matrix using a DistanceEstimation object, and then builds the phylogenetic tree using a AgglomerativeDistanceMethod object.
     * The main issue here is to estimate non-branch lengths parameters, as substitution model and rate distribution parameters.
     * Three options are provideed here:
     * - DISTANCEMETHOD_INIT (default) keep parameters to there initial value,
     * - DISTANCEMETHOD_PAIRWISE estimated parameters in a pairwise manner, which is standard but not good at all...
     * - DISTANCEMETHOD_ITERATIONS uses Ninio et al's iterative algorithm, which uses Maximum Likelihood to estimate these parameters, and then update the distance matrix.
     * Ninio M, Privman E, Pupko T, Friedman N.	
     * Phylogeny reconstruction: increasing the accuracy of pairwise distance estimation using Bayesian inference of evolutionary rates.
     * Bioinformatics. 2007 Jan 15;23(2):e136-41.
     *
     * @param estimationMethod The distance estimation object to use.
     * @param reconstructionMethod The tree reconstruction object to use.
     * @param parametersToIgnore A list of parameters to ignore while optimizing parameters.
     * @param optimizeBrLen Tell if branch lengths should be optimized together with other parameters. This may lead to more accurate parameter estimation, but is slower.
     * @param rooted Tell if rooted trees must be produced.
     * @param param String describing the type of optimization to use.
     * @param tolerance Threshold on likelihood for stopping the iterative procedure. Used only with param=DISTANCEMETHOD_ITERATIONS.
     * @param tlEvalMax Maximum number of likelihood computations to perform when optimizing parameters. Used only with param=DISTANCEMETHOD_ITERATIONS.
     * @param profiler Output stream used by optimizer. Used only with param=DISTANCEMETHOD_ITERATIONS.
     * @param messenger Output stream used by optimizer. Used only with param=DISTANCEMETHOD_ITERATIONS.
     * @param verbose Verbose level.
     */
    static TreeTemplate<Node> * buildDistanceTree(
        DistanceEstimation & estimationMethod,
        AgglomerativeDistanceMethod & reconstructionMethod,
        const ParameterList & parametersToIgnore,
        bool optimizeBrLen = false,
        bool rooted = false,
        const string & param = DISTANCEMETHOD_INIT,
        double tolerance = 0.000001,
        unsigned int tlEvalMax = 1000000,
        ostream * profiler = NULL,
        ostream * messenger = NULL,
        unsigned int verbose = 0) throw (Exception);

  public:
    static string DISTANCEMETHOD_INIT;
    static string DISTANCEMETHOD_PAIRWISE;
    static string DISTANCEMETHOD_ITERATIONS;
};


#endif	//_OPTIMIZATIONTOOLS_H_

