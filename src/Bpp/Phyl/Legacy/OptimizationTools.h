//
// File: OptimizationTools.h
// Authors:
//   Julien Dutheil
// Created: 2003-12-14 09:43:32
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef BPP_PHYL_LEGACY_OPTIMIZATIONTOOLS_H
#define BPP_PHYL_LEGACY_OPTIMIZATIONTOOLS_H

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/OutputStream.h>
#include <Bpp/Numeric/Function/SimpleNewtonMultiDimensions.h>
#include <Bpp/Phyl/OptimizationTools.h>

#include "../Parsimony/DRTreeParsimonyScore.h"
#include "../Tree/TreeTemplate.h"
#include "Likelihood/NNIHomogeneousTreeLikelihood.h"
#include "Tree/NNITopologySearch.h"

namespace bpp
{

/**
 * @brief Listener used internally by the optimizeTreeNNI method.
 */
class NNITopologyListener :
  public virtual TopologyListener
{
private:
  std::shared_ptr<NNITopologySearch> topoSearch_;
  ParameterList parameters_;
  double tolerance_;
  std::shared_ptr<OutputStream> messenger_;
  std::shared_ptr<OutputStream> profiler_;
  unsigned int verbose_;
  unsigned int optimizeCounter_;
  unsigned int optimizeNumerical_;
  std::string optMethod_;
  unsigned int nStep_;
  bool reparametrization_;

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
   * @param reparametrization Tell if parameters should be transformed in order to remove constraints.
   *                          This can improve optimization, but is a bit slower.
   */
  NNITopologyListener(
      std::shared_ptr<NNITopologySearch> ts,
      const ParameterList& parameters,
      double tolerance,
      std::shared_ptr<OutputStream> messenger,
      std::shared_ptr<OutputStream> profiler,
      unsigned int verbose,
      const std::string& optMethod,
      unsigned int nStep,
      bool reparametrization) :
    topoSearch_(ts),
    parameters_(parameters),
    tolerance_(tolerance),
    messenger_(messenger),
    profiler_(profiler),
    verbose_(verbose),
    optimizeCounter_(0),
    optimizeNumerical_(1),
    optMethod_(optMethod),
    nStep_(nStep),
    reparametrization_(reparametrization) {}

  NNITopologyListener(const NNITopologyListener& tl) :
    topoSearch_(tl.topoSearch_),
    parameters_(tl.parameters_),
    tolerance_(tl.tolerance_),
    messenger_(tl.messenger_),
    profiler_(tl.profiler_),
    verbose_(tl.verbose_),
    optimizeCounter_(tl.optimizeCounter_),
    optimizeNumerical_(tl.optimizeNumerical_),
    optMethod_(tl.optMethod_),
    nStep_(tl.nStep_),
    reparametrization_(tl.reparametrization_)
  {}

  NNITopologyListener& operator=(const NNITopologyListener& tl)
  {
    topoSearch_        = tl.topoSearch_;
    parameters_        = tl.parameters_;
    tolerance_         = tl.tolerance_;
    messenger_         = tl.messenger_;
    profiler_          = tl.profiler_;
    verbose_           = tl.verbose_;
    optimizeCounter_   = tl.optimizeCounter_;
    optimizeNumerical_ = tl.optimizeNumerical_;
    optMethod_         = tl.optMethod_;
    nStep_             = tl.nStep_;
    reparametrization_ = tl.reparametrization_;
    return *this;
  }

  NNITopologyListener* clone() const { return new NNITopologyListener(*this); }

  virtual ~NNITopologyListener() {}

public:
  void topologyChangeTested(const TopologyChangeEvent& event) {}
  void topologyChangeSuccessful(const TopologyChangeEvent& event);
  void setNumericalOptimizationCounter(unsigned int c) { optimizeNumerical_ = c; }
};

/**
 * @brief Listener used internally by the optimizeTreeNNI2 method.
 */
class NNITopologyListener2 :
  public TopologyListener
{
private:
  std::shared_ptr<NNITopologySearch> topoSearch_;
  ParameterList parameters_;
  double tolerance_;
  std::shared_ptr<OutputStream> messenger_;
  std::shared_ptr<OutputStream> profiler_;
  unsigned int verbose_;
  unsigned int optimizeCounter_;
  unsigned int optimizeNumerical_;
  std::string optMethod_;
  bool reparametrization_;

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
   * @param reparametrization Tell if parameters should be transformed in order to remove constraints.
   *                          This can improve optimization, but is a bit slower.
   */
  NNITopologyListener2(
      std::shared_ptr<NNITopologySearch> ts,
      const ParameterList& parameters,
      double tolerance,
      std::shared_ptr<OutputStream> messenger,
      std::shared_ptr<OutputStream> profiler,
      unsigned int verbose,
      const std::string& optMethod,
      bool reparametrization) :
    topoSearch_(ts),
    parameters_(parameters),
    tolerance_(tolerance),
    messenger_(messenger),
    profiler_(profiler),
    verbose_(verbose),
    optimizeCounter_(0),
    optimizeNumerical_(1),
    optMethod_(optMethod),
    reparametrization_(reparametrization) {}

  NNITopologyListener2(const NNITopologyListener2& tl) :
    topoSearch_(tl.topoSearch_),
    parameters_(tl.parameters_),
    tolerance_(tl.tolerance_),
    messenger_(tl.messenger_),
    profiler_(tl.profiler_),
    verbose_(tl.verbose_),
    optimizeCounter_(tl.optimizeCounter_),
    optimizeNumerical_(tl.optimizeNumerical_),
    optMethod_(tl.optMethod_),
    reparametrization_(tl.reparametrization_)
  {}

  NNITopologyListener2& operator=(const NNITopologyListener2& tl)
  {
    topoSearch_        = tl.topoSearch_;
    parameters_        = tl.parameters_;
    tolerance_         = tl.tolerance_;
    messenger_         = tl.messenger_;
    profiler_          = tl.profiler_;
    verbose_           = tl.verbose_;
    optimizeCounter_   = tl.optimizeCounter_;
    optimizeNumerical_ = tl.optimizeNumerical_;
    optMethod_         = tl.optMethod_;
    reparametrization_ = tl.reparametrization_;
    return *this;
  }

  NNITopologyListener2* clone() const { return new NNITopologyListener2(*this); }

  virtual ~NNITopologyListener2() {}

public:
  void topologyChangeTested(const TopologyChangeEvent& event) {}
  void topologyChangeSuccessful(const TopologyChangeEvent& event);
  void setNumericalOptimizationCounter(unsigned int c) { optimizeNumerical_ = c; }
};


/**
 * @brief Optimization methods for phylogenetic inference.
 *
 * This class provides optimization methods for phylogenetics.
 * Parameters of the optimization methods are set to work with TreeLikelihood
 * object. Some non trivial parameters are left to the user choice (tolerance, maximum
 * number of function evaluation, output streams).
 */
class LegacyOptimizationTools :
  public OptimizationTools
{
public:
  LegacyOptimizationTools();
  virtual ~LegacyOptimizationTools();

public:
  /**
   * @brief Optimize numerical parameters (branch length, substitution model & rate distribution) of a TreeLikelihood function.
   *
   * Uses Newton's method for branch length and Brent or BFGS one dimensional method for other parameters.
   *
   * A condition over function values is used as a stop condition for the algorithm.
   *
   * @see BrentOneDimension, BFGSMultiDimensions
   *
   * @param tl             A pointer toward the TreeLikelihood object to optimize.
   * @param parameters     The list of parameters to optimize. Use tl->getIndependentParameters() in order to estimate all parameters.
   * @param listener       A pointer toward an optimization listener, if needed.
   * @param nstep          The number of progressive steps to perform (see NewtonBrentMetaOptimizer). 1 means full precision from start.
   * @param tolerance      The tolerance to use in the algorithm.
   * @param tlEvalMax      The maximum number of function evaluations.
   * @param messageHandler The massage handler.
   * @param profiler       The profiler.
   * @param reparametrization Tell if parameters should be transformed in order to remove constraints.
   *                          This can improve optimization, but is a bit slower.
   * @param verbose        The verbose level.
   * @param optMethodDeriv Optimization type for derivable parameters (first or second order derivatives).
   * @see OPTIMIZATION_NEWTON, OPTIMIZATION_GRADIENT
   * @param optMethodModel Optimization type for model parameters (Brent or BFGS).
   * @see OPTIMIZATION_BRENT, OPTIMIZATION_BFGS
   * @throw Exception any exception thrown by the Optimizer.
   */
  static unsigned int optimizeNumericalParameters(
    std::shared_ptr<DiscreteRatesAcrossSitesTreeLikelihoodInterface> tl,
    const ParameterList& parameters,
    std::shared_ptr<OptimizationListener> listener = nullptr,
    unsigned int nstep                             = 1,
    double tolerance                               = 0.000001,
    unsigned int tlEvalMax                         = 1000000,
    std::shared_ptr<OutputStream> messageHandler   = ApplicationTools::message,
    std::shared_ptr<OutputStream> profiler         = ApplicationTools::message,
    bool reparametrization                         = false,
    unsigned int verbose                           = 1,
    const std::string& optMethodDeriv              = OPTIMIZATION_NEWTON,
    const std::string& optMethodModel              = OPTIMIZATION_BRENT);

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
   * @param reparametrization Tell if parameters should be transformed in order to remove constraints.
   *                          This can improve optimization, but is a bit slower.
   * @param useClock       Tell if branch lengths have to be optimized under a global molecular clock constraint.
   * @param verbose        The verbose level.
   * @param optMethodDeriv Optimization type for derivable parameters (first or second order derivatives).
   * @see OPTIMIZATION_NEWTON, OPTIMIZATION_GRADIENT
   * @throw Exception any exception thrown by the Optimizer.
   */
  static unsigned int optimizeNumericalParameters2(
    std::shared_ptr<DiscreteRatesAcrossSitesTreeLikelihoodInterface> tl,
    const ParameterList& parameters,
    std::shared_ptr<OptimizationListener> listener = nullptr,
    double tolerance                               = 0.000001,
    unsigned int tlEvalMax                         = 1000000,
    std::shared_ptr<OutputStream> messageHandler   = ApplicationTools::message,
    std::shared_ptr<OutputStream> profiler         = ApplicationTools::message,
    bool reparametrization                         = false,
    bool useClock                                  = false,
    unsigned int verbose                           = 1,
    const std::string& optMethodDeriv              = OPTIMIZATION_NEWTON);

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
   * @param optMethodDeriv Optimization type for derivable parameters (first or second order derivatives).
   * @see OPTIMIZATION_NEWTON, OPTIMIZATION_GRADIENT
   * @throw Exception any exception thrown by the Optimizer.
   */
  static unsigned int optimizeBranchLengthsParameters(
    std::shared_ptr<DiscreteRatesAcrossSitesTreeLikelihoodInterface> tl,
    const ParameterList& parameters,
    std::shared_ptr<OptimizationListener> listener = nullptr,
    double tolerance                               = 0.000001,
    unsigned int tlEvalMax                         = 1000000,
    std::shared_ptr<OutputStream> messageHandler   = ApplicationTools::message,
    std::shared_ptr<OutputStream> profiler         = ApplicationTools::message,
    unsigned int verbose                           = 1,
    const std::string& optMethodDeriv              = OPTIMIZATION_NEWTON);




private:
  class ScaleFunction :
    public virtual FunctionInterface,
    public ParametrizableAdapter
  {
private:
    std::shared_ptr<TreeLikelihoodInterface> tl_;
    mutable ParameterList brLen_, lambda_;

public:
    ScaleFunction(std::shared_ptr<TreeLikelihoodInterface> tl);

    ScaleFunction(const ScaleFunction& sf) :
      tl_(sf.tl_),
      brLen_(sf.brLen_),
      lambda_(sf.lambda_)
    {}

    ScaleFunction& operator=(const ScaleFunction& sf)
    {
      tl_     = sf.tl_;
      brLen_  = sf.brLen_;
      lambda_ = sf.lambda_;
      return *this;
    }

    virtual ~ScaleFunction();

    ScaleFunction* clone() const { return new ScaleFunction(*this); }

public:
    void setParameters(const ParameterList& lambda);
    double getValue() const;
    const ParameterList& getParameters() const { return lambda_; }
    const Parameter& getParameter(const std::string& name) const
    {
      if (name == "lambda") return lambda_[0];
      else throw ParameterNotFoundException("ScaleFunction::getParameter.", name);
    }
    double getParameterValue(const std::string& name) const
    {
      return lambda_.getParameter(name).getValue();
    }
    size_t getNumberOfParameters() const { return 1; }
    size_t getNumberOfIndependentParameters() const { return 1; }

protected:
    ParameterList& getParameters_() { return lambda_; }
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
    std::shared_ptr<TreeLikelihoodInterface> tl,
    double tolerance                             = 0.000001,
    unsigned int tlEvalMax                       = 1000000,
    std::shared_ptr<OutputStream> messageHandler = ApplicationTools::message,
    std::shared_ptr<OutputStream> profiler       = ApplicationTools::message,
    unsigned int verbose                         = 1);

  
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
   * @param tl                A pointer toward the TreeLikelihood object to optimize.
   * @param parameters        The list of parameters to optimize. Use tl->getIndependentParameters() in order to estimate all parameters.
   * @param optimizeNumFirst  Tell if we must optimize numerical parameters before searching topology.
   * @param tolBefore         The tolerance to use when estimating numerical parameters before topology search (if optimizeNumFirst is set to 'true').
   * @param tolDuring         The tolerance to use when estimating numerical parameters during the topology search.
   * @param tlEvalMax         The maximum number of function evaluations.
   * @param numStep           Number of NNI rounds before re-estimating numerical parameters.
   * @param messageHandler    The massage handler.
   * @param profiler          The profiler.
   * @param reparametrization Tell if parameters should be transformed in order to remove constraints.
   *                          This can improve optimization, but is a bit slower.
   * @param verbose           The verbose level.
   * @param optMethod         Option passed to optimizeNumericalParameters.
   * @param nStep             Option passed to optimizeNumericalParameters.
   * @param nniMethod         NNI algorithm to use.
   * @return A pointer toward the final likelihood object.
   * This pointer may be the same as passed in argument (tl), but in some cases the algorithm
   * clone this object. We may change this bahavior in the future...
   * You hence should write something like
   * @code
   * tl = OptimizationTools::optimizeTreeNNI(tl, ...);
   * @endcode
   * @throw Exception any exception thrown by the optimizer.
   */
  static std::shared_ptr<NNIHomogeneousTreeLikelihood> optimizeTreeNNI(
    std::shared_ptr<NNIHomogeneousTreeLikelihood> tl,
    const ParameterList& parameters,
    bool optimizeNumFirst                        = true,
    double tolBefore                             = 100,
    double tolDuring                             = 100,
    unsigned int tlEvalMax                       = 1000000,
    unsigned int numStep                         = 1,
    std::shared_ptr<OutputStream> messageHandler = ApplicationTools::message,
    std::shared_ptr<OutputStream> profiler       = ApplicationTools::message,
    bool reparametrization                       = false,
    unsigned int verbose                         = 1,
    const std::string& optMethod                 = OptimizationTools::OPTIMIZATION_NEWTON,
    unsigned int nStep                           = 1,
    const std::string& nniMethod                 = NNITopologySearch::PHYML);

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
   * @param tl                A pointer toward the TreeLikelihood object to optimize.
   * @param parameters        The list of parameters to optimize. Use tl->getIndependentParameters() in order to estimate all parameters.
   * @param optimizeNumFirst  Tell if we must optimize numerical parameters before searching topology.
   * @param tolBefore         The tolerance to use when estimating numerical parameters before topology search (if optimizeNumFirst is set to 'true').
   * @param tolDuring         The tolerance to use when estimating numerical parameters during the topology search.
   * @param tlEvalMax         The maximum number of function evaluations.
   * @param numStep           Number of NNI rounds before re-estimating numerical parameters.
   * @param messageHandler    The massage handler.
   * @param profiler          The profiler.
   * @param reparametrization Tell if parameters should be transformed in order to remove constraints.
   *                          This can improve optimization, but is a bit slower.
   * @param verbose           The verbose level.
   * @param optMethod         Option passed to optimizeNumericalParameters2.
   * @param nniMethod         NNI algorithm to use.
   * @return A pointer toward the final likelihood object.
   * This pointer may be the same as passed in argument (tl), but in some cases the algorithm
   * clone this object. We may change this bahavior in the future...
   * You hence should write something like
   * @code
   * tl = OptimizationTools::optimizeTreeNNI2(tl, ...);
   * @endcode
   * @throw Exception any exception thrown by the optimizer.
   */
  static std::shared_ptr<NNIHomogeneousTreeLikelihood> optimizeTreeNNI2(
    std::shared_ptr<NNIHomogeneousTreeLikelihood> tl,
    const ParameterList& parameters,
    bool optimizeNumFirst                        = true,
    double tolBefore                             = 100,
    double tolDuring                             = 100,
    unsigned int tlEvalMax                       = 1000000,
    unsigned int numStep                         = 1,
    std::shared_ptr<OutputStream> messageHandler = ApplicationTools::message,
    std::shared_ptr<OutputStream> profiler       = ApplicationTools::message,
    bool reparametrization                       = false,
    unsigned int verbose                         = 1,
    const std::string& optMethod                 = OptimizationTools::OPTIMIZATION_NEWTON,
    const std::string& nniMethod                 = NNITopologySearch::PHYML);

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
  static std::shared_ptr<DRTreeParsimonyScore> optimizeTreeNNI(
    std::shared_ptr<DRTreeParsimonyScore> tp,
    unsigned int verbose = 1);

};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_OPTIMIZATIONTOOLS_H
