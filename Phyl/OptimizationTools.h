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

#include "ClockTreeLikelihood.h"
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

namespace bpp
{

/**
 * @brief Listener used internally by the optimizeTreeNNI method.
 */
class NNITopologyListener:
  public TopologyListener
{
  protected:
    NNITopologySearch * _topoSearch;
    ParameterList parameters_;
    double _tolerance;
    ostream *_messenger;
    ostream *_profiler;
    unsigned int _verbose;
    unsigned int _optimizeCounter;
    unsigned int _optimizeNumerical;
    string _optMethod;
    unsigned int _nStep;

  public:
    /**
     * @brief Build a new NNITopologyListener object.
     *
     * This listener listens to a NNITopologySearch object, and optimizes numerical parameters every *n* topological movements.
     * Optimization is performed using the optimizeNumericalParameters method (see there documentation for more details).
     *
     * @param ts         The NNITopologySearch object attached to this listener.
     * @param parameters The list of parameters to optimize. Use tl->getIndependentParameters() in order to estimate all parameters.
     * @param tolerance  Tolerance to use during optimizaton.
     * @param messenger  Where to output messages.
     * @param profiler   Where to output optimization steps.
     * @param verbose    Verbose level during optimization.
     * @param optMethod  Optimization method to use.
     * @param nStep      The number of optimization steps to perform.
     */
    NNITopologyListener(NNITopologySearch * ts, const ParameterList& parameters, double tolerance, ostream *messenger, ostream *profiler, unsigned int verbose, const string & optMethod, unsigned int nStep):
      _topoSearch(ts), parameters_(parameters), _tolerance(tolerance),
      _messenger(messenger), _profiler(profiler),
      _verbose(verbose),
      _optimizeCounter(0), _optimizeNumerical(1),
      _optMethod(optMethod), _nStep(nStep) {}

    virtual ~NNITopologyListener() {}

  public:
    void topologyChangeTested(const TopologyChangeEvent & event) {}
    void topologyChangeSuccessful(const TopologyChangeEvent & event);
    void setNumericalOptimizationCounter(unsigned int c) { _optimizeNumerical = c; }

};

/**
 * @brief Listener used internally by the optimizeTreeNNI2 method.
 */
class NNITopologyListener2:
  public TopologyListener
{
  protected:
    NNITopologySearch * _topoSearch;
    ParameterList parameters_;
    double _tolerance;
    ostream *_messenger;
    ostream *_profiler;
    unsigned int _verbose;
    unsigned int _optimizeCounter;
    unsigned int _optimizeNumerical;
    string _optMethod;

  public:
    /**
     * @brief Build a new NNITopologyListener2 object.
     *
     * This listener listens to a NNITopologySearch object, and optimizes numerical parameters every *n* topological movements.
     * Optimization is performed using the optimizeNumericalParameters2 method (see there documentation for more details).
     *
     * @param ts         The NNITopologySearch object attached to this listener.
     * @param parameters The list of parameters to optimize. Use ts->getIndependentParameters() in order to estimate all parameters.
     * @param tolerance  Tolerance to use during optimizaton.
     * @param messenger  Where to output messages.
     * @param profiler   Where to output optimization steps.
     * @param verbose    Verbose level during optimization.
     * @param optMethod  Optimization method to use.
     */
    NNITopologyListener2(NNITopologySearch * ts, const ParameterList& parameters, double tolerance, ostream *messenger, ostream *profiler, unsigned int verbose, const string & optMethod):
      _topoSearch(ts), parameters_(parameters), _tolerance(tolerance),
      _messenger(messenger), _profiler(profiler),
      _verbose(verbose),
      _optimizeCounter(0), _optimizeNumerical(1),
      _optMethod(optMethod) {}

    virtual ~NNITopologyListener2() {}

  public:
    void topologyChangeTested(const TopologyChangeEvent & event) {}
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
     * @param parameters     The list of parameters to optimize. Use tl->getIndependentParameters() in order to estimate all parameters.
     * @param listener       A pointer toward an optimization listener, if needed.
     * @param nstep          The number of progressive steps to perform (see NewtonBrentMetaOptimizer). 1 means full precision from start.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @param verbose        The verbose level.
     * @param method         Optimization type for derivable parameters (first or second order derivatives).
     * @see OPTIMIZATION_NEWTON, OPTIMIZATION_GRADIENT
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static unsigned int optimizeNumericalParameters(
			DiscreteRatesAcrossSitesTreeLikelihood * tl,
      const ParameterList& parameters,
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
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
     * @param parameters     The list of parameters to optimize. Use tl->getIndependentParameters() in order to estimate all parameters.
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
		static unsigned int optimizeNumericalParameters2(
			DiscreteRatesAcrossSitesTreeLikelihood * tl,
      const ParameterList& parameters,
      OptimizationListener * listener = NULL,
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
     * @param parameters     The list of parameters to optimize. The intersection of branch length parameters and the input set will be used. Use tl->getBranchLengthsParameters() in order to estimate all branch length parameters.
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
      const ParameterList& parameters,
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
		 *
		 * @param cl             A pointer toward the ClockTreeLikelihood object to optimize.
     * @param parameters     The list of parameters to optimize. Use cl->getIndependentParameters() in order to estimate all parameters.
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
		static unsigned int optimizeNumericalParametersWithGlobalClock(
			DiscreteRatesAcrossSitesClockTreeLikelihood * cl,
      const ParameterList& parameters,
      OptimizationListener * listener = NULL,
      unsigned int nstep = 1,
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
		 *
		 * @param cl             A pointer toward the ClockTreeLikelihood object to optimize.
     * @param parameters     The list of parameters to optimize. Use cl->getIndependentParameters() in order to estimate all parameters.
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
		static unsigned int optimizeNumericalParametersWithGlobalClock2(
			DiscreteRatesAcrossSitesClockTreeLikelihood * cl,
      const ParameterList& parameters,
      OptimizationListener * listener = NULL,
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
      public Function,
      public ParametrizableAdapter
    {
				
			protected:
				TreeLikelihood * _tl;
				mutable ParameterList _brLen, _lambda;
				
			public:
				ScaleFunction(TreeLikelihood * tl);
				virtual ~ScaleFunction();

#ifndef NO_VIRTUAL_COV
        ScaleFunction*
#else
        Clonable*
#endif
        clone() const { return new ScaleFunction(*this); }
				
			public:
				void setParameters(const ParameterList & lambda) throw (ParameterNotFoundException, ConstraintException);
				double getValue() const throw (ParameterException);
				const ParameterList & getParameters() const throw (Exception) { return _lambda; }
        const Parameter & getParameter(const string & name) const throw (ParameterNotFoundException)
        {
          if(name == "lambda") return *_lambda[0];
          else throw ParameterNotFoundException("ScaleFunction::getParameter.", name);
        }
				double getParameterValue(const string & name) const throw (ParameterNotFoundException)
        {
          return _lambda.getParameter(name).getValue();
        }
		    unsigned int getNumberOfParameters() const { return 1; }
		    unsigned int getNumberOfIndependentParameters() const { return 1; }
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
     * @param parameters       The list of parameters to optimize. Use tl->getIndependentParameters() in order to estimate all parameters.
     * @param optimizeNumFirst Tell if we must optimize numerical parameters before searching topology.
		 * @param tolBefore        The tolerance to use when estimating numerical parameters before topology search (if optimizeNumFirst is set to 'true').
		 * @param tolDuring        The tolerance to use when estimating numerical parameters during the topology search.
		 * @param tlEvalMax        The maximum number of function evaluations.
     * @param numStep          Number of NNI rounds before re-estimating numerical parameters.
		 * @param messageHandler   The massage handler.
		 * @param profiler         The profiler.
		 * @param verbose          The verbose level.
     * @param optMethod        Option passed to optimizeNumericalParameters.
     * @param nStep            Option passed to optimizeNumericalParameters.
     * @param nniMethod        NNI algorithm to use.
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
        const ParameterList& parameters,
        bool optimizeNumFirst = true,
			  double tolBefore = 100,
			  double tolDuring = 100,
			  int tlEvalMax = 1000000,
        unsigned int numStep = 1,
			  ostream * messageHandler = &cout,
			  ostream * profiler       = &cout,
			  unsigned int verbose     = 1,
        const string & optMethod = OptimizationTools::OPTIMIZATION_NEWTON,
        unsigned int nStep       = 1,
        const string & nniMethod = NNITopologySearch::PHYML)
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
     * The optimizeNumericalParameters2 method is used for estimating numerical parameters.
     * The tolerance passed to this function is specified as input parameters.
     * They are generally very high to avoid local optima.
     *
		 * @param tl               A pointer toward the TreeLikelihood object to optimize.
     * @param parameters       The list of parameters to optimize. Use tl->getIndependentParameters() in order to estimate all parameters.
     * @param optimizeNumFirst Tell if we must optimize numerical parameters before searching topology.
		 * @param tolBefore        The tolerance to use when estimating numerical parameters before topology search (if optimizeNumFirst is set to 'true').
		 * @param tolDuring        The tolerance to use when estimating numerical parameters during the topology search.
		 * @param tlEvalMax        The maximum number of function evaluations.
     * @param numStep          Number of NNI rounds before re-estimating numerical parameters.
		 * @param messageHandler   The massage handler.
		 * @param profiler         The profiler.
		 * @param verbose          The verbose level.
     * @param optMethod        Option passed to optimizeNumericalParameters2.
     * @param nniMethod        NNI algorithm to use.
     * @return A pointer toward the final likelihood object.
     * This pointer may be the same as passed in argument (tl), but in some cases the algorithm
     * clone this object. We may change this bahavior in the future...
     * You hence should write something like
     * @code
     * tl = OptimizationTools::optimizeTreeNNI2(tl, ...);
     * @endcode
		 * @throw Exception any exception thrown by the optimizer.
     */
    static NNIHomogeneousTreeLikelihood * optimizeTreeNNI2(
        NNIHomogeneousTreeLikelihood * tl,
        const ParameterList& parameters,
        bool optimizeNumFirst = true,
			  double tolBefore = 100,
			  double tolDuring = 100,
			  int tlEvalMax = 1000000,
        unsigned int numStep = 1,
			  ostream * messageHandler = &cout,
			  ostream * profiler       = &cout,
			  unsigned int verbose     = 1,
        const string & optMethod = OptimizationTools::OPTIMIZATION_NEWTON,
        const string & nniMethod = NNITopologySearch::PHYML)
      throw (Exception);

    /**
     * @brief Optimize tree topology from a DRTreeParsimonyScore using Nearest Neighbor Interchanges.
     *
		 * @param tp               A pointer toward the DRTreeParsimonyScore object to optimize.
		 * @param verbose          The verbose level.
     * @return A pointer toward the final parsimony score object.
     * This pointer may be the same as passed in argument (tl), but in some cases the algorithm
     * clone this object. We may change this bahavior in the future...
     * You hence should write something like
     * @code
     * tp = OptimizationTools::optimizeTreeNNI(tp, ...);
     * @endcode
     */
    
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


} //end of namespace bpp.

#endif	//_OPTIMIZATIONTOOLS_H_

