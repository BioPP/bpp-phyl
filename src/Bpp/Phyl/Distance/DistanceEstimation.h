// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_DISTANCE_DISTANCEESTIMATION_H
#define BPP_PHYL_DISTANCE_DISTANCEESTIMATION_H

#include <Bpp/Clonable.h>
#include <Bpp/Numeric/Function/MetaOptimizer.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Function/SimpleMultiDimensions.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

#include "../Likelihood/SubstitutionProcess.h"
#include "../PseudoNewtonOptimizer.h"
#include "../Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h"

// From bpp-seq:
#include <Bpp/Seq/Container/AlignmentData.h>

// From the STL:
#include <memory>

namespace bpp
{
class DistanceEstimation :
  public virtual Clonable
{
private:
  std::shared_ptr<SubstitutionProcessInterface> process_;
  size_t numProc_; // needed if in Collection Process
  std::shared_ptr<const AlignmentDataInterface> sites_;
  std::shared_ptr<DistanceMatrix> dist_;
  std::shared_ptr<OptimizerInterface> optimizer_;
  std::shared_ptr<MetaOptimizer> defaultOptimizer_;
  size_t verbose_;
  ParameterList parameters_;

public:
  /**
   * @brief Create a new DistanceEstimation object according to a given substitution model and a rate distribution.
   *
   * @param process  The substitution process to use.
   * @param rateDist The discrete rate distribution to use.
   * @param verbose  The verbose level:
   *  - 0=Off,
   *  - 1=one * by row computation
   *  - 2=one * by row computation and one . by column computation
   *  - 3=2 + optimization verbose enabled
   *  - 4=3 + likelihood object verbose enabled
   */
  DistanceEstimation(
      std::shared_ptr<SubstitutionProcessInterface> process,
      size_t verbose = 1) :
    process_(process),
    numProc_(0),
    sites_(0),
    dist_(0),
    optimizer_(0),
    defaultOptimizer_(0),
    verbose_(verbose),
    parameters_()
  {
    init_();
  }

  /**
   * @brief Create a new DistanceEstimation object and compute distances
   * according to a given substitution model and a rate distribution.
   *
   * This instance will own the model and distribution, and will take car of their recopy and destruction.
   *
   * @param process  The substitution process to use.
   * @param sites    The sequence data.
   * @param verbose  The verbose level:
   *  - 0=Off,
   *  - 1=one * by row computation
   *  - 2=one * by row computation and one . by column computation
   *  - 3=2 + optimization verbose enabled
   *  - 4=3 + likelihood object verbose enabled
   *  @param computeMat if true the computeMatrix() method is called.
   */
  DistanceEstimation(
      std::shared_ptr<SubstitutionProcessInterface> process,
      std::shared_ptr<const AlignmentDataInterface> sites,
      size_t verbose = 1,
      bool computeMat = true) :
    process_(process),
    numProc_(0),
    sites_(sites),
    dist_(0),
    optimizer_(0),
    defaultOptimizer_(0),
    verbose_(verbose),
    parameters_()
  {
    init_();
    if (computeMat) computeMatrix();
  }

  /**
   * @brief Copy constructor.
   *
   * Only the distance matrix is hard-copied, if there is one.
   *
   * @param distanceEstimation The object to copy.
   */
  DistanceEstimation(const DistanceEstimation& distanceEstimation) :
    process_(distanceEstimation.process_),
    numProc_(distanceEstimation.numProc_),
    sites_(distanceEstimation.sites_),
    dist_(0),
    optimizer_(distanceEstimation.optimizer_),
    defaultOptimizer_(distanceEstimation.defaultOptimizer_),
    verbose_(distanceEstimation.verbose_),
    parameters_(distanceEstimation.parameters_)
  {
    if (distanceEstimation.dist_ != 0)
      dist_ = std::make_shared<DistanceMatrix>(*distanceEstimation.dist_);
    else
      dist_ = 0;
  }

  /**
   * @brief Assignment operator.
   *
   * Only the distance matrix is hard-copied, if there is one.
   *
   * @param distanceEstimation The object to copy.
   * @return A reference toward this object.
   */
  DistanceEstimation& operator=(const DistanceEstimation& distanceEstimation)
  {
    process_    = distanceEstimation.process_;
    numProc_    = distanceEstimation.numProc_;
    sites_      = distanceEstimation.sites_;
    if (distanceEstimation.dist_)
      dist_     = std::make_shared<DistanceMatrix>(*distanceEstimation.dist_);
    else
      dist_     = 0;
    optimizer_  = distanceEstimation.optimizer_;
    // _defaultOptimizer has already been initialized since the default constructor has been called.
    verbose_    = distanceEstimation.verbose_;
    parameters_ = distanceEstimation.parameters_;
    return *this;
  }

  virtual ~DistanceEstimation() {}

  DistanceEstimation* clone() const override { return new DistanceEstimation(*this); }

private:
  void init_();

public:
  /**
   * @brief Perform the distance computation.
   *
   * Result can be called by the getMatrix() method.
   *
   * @throw NullPointerException if at least one of the model,
   * rate distribution or data are not initialized.
   */
  void computeMatrix();

  /**
   * @brief Get the distance matrix.
   *
   * @return A pointer toward the computed distance matrix.
   */
  std::unique_ptr<DistanceMatrix> getMatrix() const
  {
    return dist_ == nullptr ? nullptr : std::make_unique<DistanceMatrix>(*dist_);
  }

  bool hasProcess() const { return process_.get(); }

  const SubstitutionProcessInterface& process() const
  {
    if (hasProcess())
      return *process_;
    else
      throw Exception("DistanceEstimation::getSubstitutionModel(). No process associated to this instance.");
  }

  std::shared_ptr<const SubstitutionProcessInterface> getProcess() const
  {
    return process_;
  }

  void setProcess(std::shared_ptr<SubstitutionProcessInterface> process = nullptr) { process_ = process; }

  void setData(std::shared_ptr<const AlignmentDataInterface> sites = nullptr) { sites_ = sites; }

  std::shared_ptr<const AlignmentDataInterface> getData() const { return sites_; }

  const AlignmentDataInterface& data() const { return *sites_; }

  void setOptimizer(std::shared_ptr<OptimizerInterface> optimizer)
  {
    optimizer_ = optimizer;
  }

  std::shared_ptr<const OptimizerInterface> getOptimizer() const { return optimizer_; }

  std::shared_ptr<OptimizerInterface> getOptimizer() { return optimizer_; }

  const OptimizerInterface& optimizer() const { return *optimizer_; }

  OptimizerInterface& optimizer() { return *optimizer_; }

  void resetOptimizer() { optimizer_ = dynamic_pointer_cast<OptimizerInterface>(defaultOptimizer_); }

  bool matchParametersValues(const ParameterList& parameters)
  {
    if (hasProcess())
      return process_->matchParametersValues(parameters);
    return false;
  }
  
  /**
   * @brief Specify a list of parameters to be estimated.
   *
   * Parameters will be estimated separately for each distance.
   *
   * @param parameters A list of parameters to estimate.
   */
  void setAdditionalParameters(const ParameterList& parameters)
  {
    parameters_ = parameters;
  }

  /**
   * @brief Reset all additional parameters.
   */
  void resetAdditionalParameters()
  {
    parameters_.reset();
  }

  /**
   * @param verbose Verbose level.
   */
  void setVerbose(size_t verbose) { verbose_ = verbose; }
  /**
   * @return Verbose level.
   */
  size_t getVerbose() const { return verbose_; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_DISTANCE_DISTANCEESTIMATION_H
