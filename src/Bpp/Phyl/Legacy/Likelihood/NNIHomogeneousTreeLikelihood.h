// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_NNIHOMOGENEOUSTREELIKELIHOOD_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_NNIHOMOGENEOUSTREELIKELIHOOD_H

#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/VectorTools.h>

#include "../../Tree/NNISearchable.h"
#include "DRHomogeneousTreeLikelihood.h"

namespace bpp
{
/**
 * @brief Compute likelihood for a 4-tree.
 *
 * This class is used internally by DRHomogeneousTreeLikelihood to test NNI movements.
 * This function needs:
 * - two likelihood arrays corresponding to the conditional likelihoods at top and bottom nodes,
 * - a substitution model and a rate distribution, whose parameters will not be estimated but taken "as is",
 * It takes only one parameter, the branch length.
 */
class BranchLikelihood :
  public FunctionInterface,
  public AbstractParametrizable
{
protected:
  const VVVdouble* array1_, * array2_;
  std::shared_ptr<const TransitionModelInterface> model_;
  std::shared_ptr<const DiscreteDistributionInterface> rDist_;
  size_t nbStates_, nbClasses_;
  VVVdouble pxy_;
  double lnL_;
  std::vector<unsigned int> weights_;

public:
  BranchLikelihood(const std::vector<unsigned int>& weights) :
    AbstractParametrizable(""),
    array1_(0),
    array2_(0),
    model_(0),
    rDist_(0),
    nbStates_(0),
    nbClasses_(0),
    pxy_(),
    lnL_(std::log(0.)),
    weights_(weights)
  {
    addParameter_(new Parameter("BrLen", 1, 0));
  }

  BranchLikelihood(const BranchLikelihood& bl) :
    AbstractParametrizable(bl),
    array1_(bl.array1_),
    array2_(bl.array2_),
    model_(bl.model_),
    rDist_(bl.rDist_),
    nbStates_(bl.nbStates_),
    nbClasses_(bl.nbClasses_),
    pxy_(bl.pxy_),
    lnL_(bl.lnL_),
    weights_(bl.weights_)
  {}

  BranchLikelihood& operator=(const BranchLikelihood& bl)
  {
    AbstractParametrizable::operator=(bl);
    array1_ = bl.array1_;
    array2_ = bl.array2_;
    model_ = bl.model_;
    rDist_ = bl.rDist_;
    nbStates_ = bl.nbStates_;
    nbClasses_ = bl.nbClasses_;
    pxy_ = bl.pxy_;
    lnL_ = bl.lnL_;
    weights_ = bl.weights_;
    return *this;
  }

  virtual ~BranchLikelihood() {}

  BranchLikelihood* clone() const { return new BranchLikelihood(*this); }

public:
  void initModel(
      std::shared_ptr<const TransitionModelInterface> model,
      std::shared_ptr<const DiscreteDistributionInterface> rDist);

  /**
   * @warning No checking on alphabet size or number of rate classes is performed,
   * use with care!
   */
  void initLikelihoods(const VVVdouble* array1, const VVVdouble* array2)
  {
    array1_ = array1;
    array2_ = array2;
  }

  void resetLikelihoods()
  {
    array1_ = 0;
    array2_ = 0;
  }

  void setParameters(const ParameterList& parameters)
  {
    setParametersValues(parameters);
  }

  double getValue() const { return lnL_; }

  void fireParameterChanged(const ParameterList& parameters)
  {
    computeAllTransitionProbabilities();
    computeLogLikelihood();
  }

protected:
  void computeAllTransitionProbabilities();
  void computeLogLikelihood();
};


/**
 * @brief This class adds support for NNI topology estimation to the DRHomogeneousTreeLikelihood class.
 */
class NNIHomogeneousTreeLikelihood :
  public DRHomogeneousTreeLikelihood,
  public virtual NNISearchable
{
protected:
  std::shared_ptr<BranchLikelihood> brLikFunction_;

  /**
   * @brief Optimizer used for testing NNI.
   */
  std::unique_ptr<BrentOneDimension> brentOptimizer_;

  /**
   * @brief Hash used for backing up branch lengths when testing NNIs.
   */
  mutable std::map<int, double> brLenNNIValues_;

  ParameterList brLenNNIParams_;

public:
  /**
   * @brief Build a new NNIHomogeneousTreeLikelihood object.
   *
   * @param tree The tree to use.
   * @param model The substitution model to use.
   * @param rDist The rate across sites distribution to use.
   * @param checkRooted Tell if we have to check for the tree to be unrooted.
   * If true, any rooted tree will be unrooted before likelihood computation.
   * @param verbose Should I display some info?
   * @throw Exception in an error occured.
   */
  NNIHomogeneousTreeLikelihood(
    const Tree& tree,
    std::shared_ptr<TransitionModelInterface> model,
    std::shared_ptr<DiscreteDistributionInterface> rDist,
    bool checkRooted = true,
    bool verbose = true);

  /**
   * @brief Build a new NNIHomogeneousTreeLikelihood object.
   *
   * @param tree The tree to use.
   * @param data Sequences to use.
   * @param model The substitution model to use.
   * @param rDist The rate across sites distribution to use.
   * @param checkRooted Tell if we have to check for the tree to be unrooted.
   * If true, any rooted tree will be unrooted before likelihood computation.
   * @param verbose Should I display some info?
   * @throw Exception in an error occured.
   */
  NNIHomogeneousTreeLikelihood(
    const Tree& tree,
    const AlignmentDataInterface& data,
    std::shared_ptr<TransitionModelInterface> model,
    std::shared_ptr<DiscreteDistributionInterface> rDist,
    bool checkRooted = true,
    bool verbose = true);

  /**
   * @brief Copy constructor.
   */
  NNIHomogeneousTreeLikelihood(const NNIHomogeneousTreeLikelihood& lik);

  NNIHomogeneousTreeLikelihood& operator=(const NNIHomogeneousTreeLikelihood& lik);

  virtual ~NNIHomogeneousTreeLikelihood();

  NNIHomogeneousTreeLikelihood* clone() const override { return new NNIHomogeneousTreeLikelihood(*this); }

public:
  void setData(const AlignmentDataInterface& sites) override
  {
    DRHomogeneousTreeLikelihood::setData(sites);
    brLikFunction_ = std::make_unique<BranchLikelihood>(likelihoodData().getWeights());
  }

  /**
   * @name The NNISearchable interface.
   *
   * Current implementation:
   * When testing a particular NNI, only the branch length of the parent node is optimized (and roughly).
   * All other parameters (substitution model, rate distribution and other branch length are kept at there current value.
   * When performing a NNI, only the topology change is performed.
   * This is up to the user to re-initialize the underlying likelihood data to match the new topology.
   * Usually, this is achieved by calling the topologyChangePerformed() method, which call the reInit() method of the LikelihoodData object.
   * @{
   */
  const Tree& topology() const override { return tree(); }

  double getTopologyValue() const override { return getValue(); }

  double testNNI(int nodeId) const override;

  void doNNI(int nodeId) override;

  void topologyChangeTested(const TopologyChangeEvent& event) override
  {
    likelihoodData().reInit();
    // if(brLenNNIParams_.size() > 0)
    fireParameterChanged(brLenNNIParams_);
    brLenNNIParams_.reset();
  }

  void topologyChangeSuccessful(const TopologyChangeEvent& event) override
  {
    brLenNNIValues_.clear();
  }
  /** @} */
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_NNIHOMOGENEOUSTREELIKELIHOOD_H
