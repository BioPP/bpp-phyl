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
#include <Bpp/Seq/Container/AlignedValuesContainer.h>

// From the STL:
#include <memory>

namespace bpp
{
  class DistanceEstimation :
  public virtual Clonable
{
private:
  std::shared_ptr<BranchModel> model_;
  std::shared_ptr<DiscreteDistribution> rateDist_;
  const AlignedValuesContainer* sites_;
  std::shared_ptr<DistanceMatrix> dist_;
  Optimizer* optimizer_;
  MetaOptimizer* defaultOptimizer_;
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
    TransitionModel* model,
    DiscreteDistribution* rateDist,
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
    TransitionModel* model,
    DiscreteDistribution* rateDist,
    const AlignedValuesContainer* sites,
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
    model_(distanceEstimation.model_->clone()),
    rateDist_(distanceEstimation.rateDist_->clone()),
    sites_(distanceEstimation.sites_),
    dist_(0),
    optimizer_(dynamic_cast<Optimizer*>(distanceEstimation.optimizer_->clone())),
    defaultOptimizer_(dynamic_cast<MetaOptimizer*>(distanceEstimation.defaultOptimizer_->clone())),
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
    model_.reset(distanceEstimation.model_->clone());
    rateDist_.reset(distanceEstimation.rateDist_->clone());
    sites_      = distanceEstimation.sites_;
    if (distanceEstimation.dist_ != 0)
      dist_     = std::make_shared<DistanceMatrix>(*distanceEstimation.dist_);
    else
      dist_     = 0;
    optimizer_  = dynamic_cast<Optimizer*>(distanceEstimation.optimizer_->clone());
    // _defaultOptimizer has already been initialized since the default constructor has been called.
    verbose_    = distanceEstimation.verbose_;
    parameters_ = distanceEstimation.parameters_;
    return *this;
  }

  virtual ~DistanceEstimation()
  {
    delete defaultOptimizer_;
    delete optimizer_;
  }

  DistanceEstimation* clone() const { return new DistanceEstimation(*this); }

private:
  void init_()
  {
    MetaOptimizerInfos* desc = new MetaOptimizerInfos();
    std::vector<std::string> name;
    name.push_back("BrLen0");
    name.push_back("BrLen1");
    desc->addOptimizer("Branch length", new PseudoNewtonOptimizer(0), name, 2, MetaOptimizerInfos::IT_TYPE_FULL);
    ParameterList tmp = model_->getParameters();
    tmp.addParameters(rateDist_->getParameters());
    desc->addOptimizer("substitution model and rate distribution", new SimpleMultiDimensions(0), tmp.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);

    defaultOptimizer_ = new MetaOptimizer(0, desc);
    defaultOptimizer_->setMessageHandler(0);
    defaultOptimizer_->setProfiler(0);
    defaultOptimizer_->getStopCondition()->setTolerance(0.0001);
    optimizer_ = dynamic_cast<Optimizer*>(defaultOptimizer_->clone());
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
  
  DistanceMatrix* getMatrix() const { return dist_ == 0 ? 0 : new DistanceMatrix(*dist_); }

  bool hasModel() const { return model_.get(); }

  const BranchModel& getModel() const
  {
    if (hasModel())
      return *model_;
    else
      throw Exception("DistanceEstimation::getSubstitutionModel(). No model assciated to this instance.");
  }

  void resetSubstitutionModel(TransitionModel* model = 0) { model_.reset(model); }

  bool hasRateDistribution() const { return rateDist_.get(); }

  const DiscreteDistribution& getRateDistribution() const
  {
    if (hasRateDistribution())
      return *rateDist_;
    else
      throw Exception("DistanceEstimation::getRateDistribution(). No rate distribution assciated to this instance.");
  }

  void resetRateDistribution(DiscreteDistribution* rateDist = 0) { rateDist_.reset(rateDist); }

  void setData(const AlignedValuesContainer* sites) { sites_ = sites; }
  const AlignedValuesContainer* getData() const { return sites_; }
  void resetData() { sites_ = 0; }

  void setOptimizer(const Optimizer* optimizer)
  {
    if (optimizer_) delete optimizer_;
    optimizer_ = dynamic_cast<Optimizer*>(optimizer->clone());
  }
  const Optimizer* getOptimizer() const { return optimizer_; }
  Optimizer* getOptimizer() { return optimizer_; }
  void resetOptimizer() { optimizer_ = dynamic_cast<Optimizer*>(defaultOptimizer_->clone()); }

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
