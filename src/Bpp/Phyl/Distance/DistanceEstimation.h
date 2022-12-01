//
// File: DistanceEstimation.h
// Authors:
//   Julien Dutheil
//   Vincent Ranwez
// Created: 2005-06-08 10:39:00
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

#ifndef BPP_PHYL_DISTANCE_DISTANCEESTIMATION_H
#define BPP_PHYL_DISTANCE_DISTANCEESTIMATION_H

#include <Bpp/Clonable.h>
#include <Bpp/Numeric/Function/MetaOptimizer.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Function/SimpleMultiDimensions.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

#include "../Model/SubstitutionModel.h"
#include "../PseudoNewtonOptimizer.h"
#include "../Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h"
#include "../Likelihood/RateAcrossSitesSubstitutionProcess.h"

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
  std::shared_ptr<BranchModelInterface> model_;
  std::shared_ptr<DiscreteDistribution> rateDist_;
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
   * This instance will own the model and distribution, and will take car of their recopy and destruction.
   *
   * @param model    The substitution model to use.
   * @param rateDist The discrete rate distribution to use.
   * @param verbose  The verbose level:
   *  - 0=Off,
   *  - 1=one * by row computation
   *  - 2=one * by row computation and one . by column computation
   *  - 3=2 + optimization verbose enabled
   *  - 4=3 + likelihood object verbose enabled
   */
  DistanceEstimation(
    std::shared_ptr<BranchModelInterface> model,
    std::shared_ptr<DiscreteDistribution> rateDist,
    size_t verbose = 1) :
    model_(model),
    rateDist_(rateDist),
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
   * @param model    The substitution model to use.
   * @param rateDist The discrete rate distribution to use.
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
    std::shared_ptr<BranchModelInterface> model,
    std::shared_ptr<DiscreteDistribution> rateDist,
    std::shared_ptr<const AlignmentDataInterface> sites,
    size_t verbose = 1,
    bool computeMat = true) :
    model_(model),
    rateDist_(rateDist),
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
    model_(distanceEstimation.model_),
    rateDist_(distanceEstimation.rateDist_),
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
   * @brief Assigment operator.
   *
   * Only the distance matrix is hard-copied, if there is one.
   *
   * @param distanceEstimation The object to copy.
   * @return A reference toward this object.
   */
  DistanceEstimation& operator=(const DistanceEstimation& distanceEstimation)
  {
    model_      = distanceEstimation.model_;
    rateDist_   = distanceEstimation.rateDist_;
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
  void init_()
  {
    auto desc = make_shared<MetaOptimizerInfos>();
    std::vector<std::string> name;
    name.push_back("BrLen0");
    name.push_back("BrLen1");
    desc->addOptimizer("Branch length", std::make_shared<PseudoNewtonOptimizer>(nullptr), name, 2, MetaOptimizerInfos::IT_TYPE_FULL);
    ParameterList tmp = model_->getParameters();
    tmp.addParameters(rateDist_->getParameters());
    desc->addOptimizer("substitution model and rate distribution", std::make_shared<SimpleMultiDimensions>(nullptr), tmp.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);

    defaultOptimizer_ = std::make_shared<MetaOptimizer>(nullptr, desc);
    defaultOptimizer_->setMessageHandler(nullptr);
    defaultOptimizer_->setProfiler(nullptr);
    defaultOptimizer_->getStopCondition()->setTolerance(0.0001);
    optimizer_ = dynamic_pointer_cast<OptimizerInterface>(defaultOptimizer_);
  }

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

  bool hasModel() const { return model_.get(); }

  const BranchModelInterface& model() const
  {
    if (hasModel())
      return *model_;
    else
      throw Exception("DistanceEstimation::getSubstitutionModel(). No model associated to this instance.");
  }

  std::shared_ptr<const BranchModelInterface> getModel() const
  {
    return model_;
  }

  void setModel(std::shared_ptr<BranchModelInterface> model = nullptr) { model_ = model; }

  bool hasRateDistribution() const { return rateDist_.get(); }

  const DiscreteDistribution& rateDistribution() const
  {
    if (hasRateDistribution())
      return *rateDist_;
    else
      throw Exception("DistanceEstimation::getRateDistribution(). No rate distribution associated to this instance.");
  }

  std::shared_ptr<const DiscreteDistribution> getRateDistribution() const
  {
     return rateDist_;
  }
  
  void setRateDistribution(std::shared_ptr<DiscreteDistribution> rateDist = nullptr) { rateDist_ = rateDist; }

  void setData(std::shared_ptr<const AlignmentDataInterface> sites = nullptr) { sites_ = sites; }

  std::shared_ptr<const AlignmentDataInterface> getData() const { return sites_; }

  void setOptimizer(std::shared_ptr<OptimizerInterface> optimizer)
  {
    optimizer_ = optimizer;
  }

  std::shared_ptr<const OptimizerInterface> getOptimizer() const { return optimizer_; }
  
  std::shared_ptr<OptimizerInterface> getOptimizer() { return optimizer_; }

  void resetOptimizer() { optimizer_ = dynamic_pointer_cast<OptimizerInterface>(defaultOptimizer_); }

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
