// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_OPTIMIZATIONTOOLS_H
#define BPP_PHYL_OPTIMIZATIONTOOLS_H

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/OutputStream.h>
#include <Bpp/Numeric/Function/SimpleNewtonMultiDimensions.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>

#include "Distance/DistanceEstimation.h"
#include "Distance/AbstractAgglomerativeDistanceMethod.h"


namespace bpp
{
/**
 * @brief A listener which capture NaN function values and throw an exception in case this happens.
 */
class NaNListener : public OptimizationListener
{
private:
  OptimizerInterface* optimizer_;
  FunctionInterface* function_;

public:
  NaNListener(OptimizerInterface* optimizer, FunctionInterface* function) : optimizer_(optimizer), function_(function) {}

  NaNListener(const NaNListener& lr) :
    optimizer_(lr.optimizer_),
    function_(lr.function_)
  {}

  NaNListener& operator=(const NaNListener& lr)
  {
    optimizer_ = lr.optimizer_;
    function_  = lr.function_;
    return *this;
  }

public:
  void optimizationInitializationPerformed(const OptimizationEvent& event) {}
  void optimizationStepPerformed(const OptimizationEvent& event)
  {
    if (std::isnan(optimizer_->getFunction()->getValue()))
    {
      cerr << "Oups... something abnormal happened!" << endl;
      function_->getParameters().printParameters(cerr);
      throw Exception("Optimization failed because likelihood function returned NaN.");
    }
  }
  bool listenerModifiesParameters () const { return false; }
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
  static std::string OPTIMIZATION_GRADIENT;
  static std::string OPTIMIZATION_NEWTON;
  static std::string OPTIMIZATION_BRENT;
  static std::string OPTIMIZATION_BFGS;

  /**
   * Class used to register Optimization options
   * \c parameters     The list of parameters to optimize. Use tl->getIndependentParameters() in order to estimate all parameters.
   * \c listener       A pointer toward an optimization listener, if needed.
   * \c nstep          The number of progressive steps to perform (see NewtonBrentMetaOptimizer). 1 means full precision from start.
   * \c tolerance      The tolerance to use in the algorithm.
   * \c nbEvalMax      The maximum number of function evaluations.
   * \c messenger The massage handler.
   * \c profiler       The profiler.
   * \c reparametrization Tell if parameters should be transformed in order to remove constraints.
   *                          This can improve optimization, but is a bit slower.
   * \c useClock       True if follows the tree clock
   * \c verbose        The verbose level.
   * \c optMethodDeriv Optimization type for derivable parameters (first or second order derivatives).
   * @see OPTIMIZATION_NEWTON, OPTIMIZATION_GRADIENT
   * \c optMethodModel Optimization type for model parameters (Brent or BFGS).
   * @see OPTIMIZATION_BRENT, OPTIMIZATION_BFGS
   * @throw Exception any exception thrown by the Optimizer.
   **/
  
  class OptimizationOptions
  {
  public:
    ParameterList parameters;
    std::shared_ptr<OptimizationListener> listener;
    std::string backupFile;
    unsigned int nstep;
    double tolerance; 
    unsigned int nbEvalMax;
    std::shared_ptr<OutputStream> messenger;
    std::shared_ptr<OutputStream> profiler;
    bool reparametrization;
    bool useClock;
    unsigned int verbose;
    std::string optMethodDeriv;
    std::string optMethodModel;


    OptimizationOptions() :
      parameters(),
      listener(nullptr),
      backupFile("none"),
      nstep(1),
      tolerance(0.000001),
      nbEvalMax(1000000),
      messenger(ApplicationTools::message),
      profiler(ApplicationTools::message),
      reparametrization(false),
      useClock(0),
      verbose(1),
      optMethodDeriv(OPTIMIZATION_NEWTON),
      optMethodModel(OPTIMIZATION_BRENT){};

    OptimizationOptions(
      std::shared_ptr<PhyloLikelihoodInterface> lik,
      const std::map<std::string, std::string>& params,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verb = true,
      int warn = 1);

  };
  
  /**
   * @brief Optimize numerical parameters (branch length, substitution model & rate distribution) of a TreeLikelihood function.
   *
   * Uses Newton's method for branch length and Brent or BFGS one dimensional method for other parameters.
   *
   * A condition over function values is used as a stop condition for the algorithm.
   *
   * @see BrentOneDimension, BFGSMultiDimensions
   *
   * @param lik             A pointer toward the PhyloLikelihood object to optimize.
   * @param optopt   Optimization Options
   */
  
  static unsigned int optimizeNumericalParameters(
      std::shared_ptr<PhyloLikelihoodInterface> lik,
      const OptimizationOptions& optopt);
  
  /**
   * @brief Optimize numerical parameters (branch length, substitution model & rate distribution) of a TreeLikelihood function.
   *
   * Uses Newton's method for all parameters, branch length derivatives are computed analytically, derivatives for other parameters numerically.
   *
   * @see PseudoNewtonOptimizer
   *
   * @param lik            A pointer toward the PhyloLikelihood object to optimize.
   * @param optopt   Optimization Options
   * @throw Exception any exception thrown by the Optimizer.
   */

  static unsigned int optimizeNumericalParameters2(
      std::shared_ptr<PhyloLikelihoodInterface> lik,
      const OptimizationOptions& optopt);

  static unsigned int optimizeNumericalParameters2(
      std::shared_ptr<SingleProcessPhyloLikelihood> lik,
      const OptimizationOptions& optopt);

  /**
   * @brief Estimate a distance matrix using maximum likelihood.
   *
   * This method estimate a distance matrix using a DistanceEstimation object.
   * The main issue here is to estimate non-branch lengths parameters, as substitution model and rate distribution parameters.
   * Twoe options are provideed here:
   * - DISTANCEMETHOD_INIT (default) keep parameters to there initial value,
   * - DISTANCEMETHOD_PAIRWISE estimated parameters in a pairwise manner, which is standard but not that satisfying...
   *
   * @param estimationMethod The distance estimation object to use.
   * @param parametersToIgnore A list of parameters to ignore while optimizing parameters.
   * @param param String describing the type of optimization to use.
   * @param verbose Verbose level.
   *
   * @see buildDistanceTree for a procedure to jointly estimate the distance matrix and underlying tree.
   */
  static std::unique_ptr<DistanceMatrix> estimateDistanceMatrix(
      DistanceEstimation& estimationMethod,
      const ParameterList& parametersToIgnore,
      const std::string& param = DISTANCEMETHOD_INIT,
      unsigned int verbose = 0);

  /**
   * @brief Build a tree using a distance method.
   *
   * This method estimate a distance matrix using a DistanceEstimation object, and then builds the phylogenetic tree using a AgglomerativeDistanceMethod object.
   * The main issue here is to estimate non-branch lengths parameters, as substitution model and rate distribution parameters.
   * Three options are provideed here:
   * - DISTANCEMETHOD_INIT (default) keep parameters to there initial value,
   * - DISTANCEMETHOD_PAIRWISE estimated parameters in a pairwise manner, which is standard but not that satisfying...
   * - DISTANCEMETHOD_ITERATIONS uses Ninio et al's iterative algorithm, which uses Maximum Likelihood to estimate these parameters, and then update the distance matrix.
   * Ninio M, Privman E, Pupko T, Friedman N.
   * Phylogeny reconstruction: increasing the accuracy of pairwise distance estimation using Bayesian inference of evolutionary rates.
   * Bioinformatics. 2007 Jan 15;23(2):e136-41.
   *
   * @param estimationMethod The distance estimation object to use.
   * @param reconstructionMethod The tree reconstruction object to use.
   * @param param  Method of for distance calculation
   * @param optopt   Optimization Options
   */

  static std::unique_ptr<TreeTemplate<Node>> buildDistanceTree(
      DistanceEstimation& estimationMethod,
      AgglomerativeDistanceMethodInterface& reconstructionMethod,
      const std::string& param,
      OptimizationOptions& optopt);

public:
  static std::string DISTANCEMETHOD_INIT;
  static std::string DISTANCEMETHOD_PAIRWISE;
  static std::string DISTANCEMETHOD_ITERATIONS;
};
} // end of namespace bpp.
#endif // BPP_PHYL_OPTIMIZATIONTOOLS_H
