// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_DRHOMOGENEOUSTREELIKELIHOOD_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_DRHOMOGENEOUSTREELIKELIHOOD_H

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/VectorTools.h>

#include "AbstractHomogeneousTreeLikelihood.h"
#include "DRASDRTreeLikelihoodData.h"
#include "DRTreeLikelihood.h"

namespace bpp
{
/**
 * @brief This class implements the likelihood computation for a tree using the double-recursive
 * algorithm.
 *
 * The substitution model is the same over the tree (homogeneous model).
 * A non-uniform distribution of rates among the sites is allowed (ASRV models).</p>
 *
 * This class uses an instance of the DRASDRTreeLikelihoodData for conditionnal likelihood storage.
 *
 * All nodes share the same site patterns.
 */
class DRHomogeneousTreeLikelihood :
  public AbstractHomogeneousTreeLikelihood,
  public virtual DRTreeLikelihoodInterface
{
private:
  mutable std::unique_ptr<DRASDRTreeLikelihoodData> likelihoodData_;

protected:
  double minusLogLik_;

public:
  /**
   * @brief Build a new DRHomogeneousTreeLikelihood object without data.
   *
   * This constructor only initialize the parameters.
   * To compute a likelihood, you will need to call the setData() and the computeTreeLikelihood() methods.
   *
   * @param tree The tree to use.
   * @param model The substitution model to use.
   * @param rDist The rate across sites distribution to use.
   * @param checkRooted Tell if we have to check for the tree to be unrooted.
   * If true, any rooted tree will be unrooted before likelihood computation.
   * @param verbose Should I display some info?
   * @throw Exception in an error occured.
   */
  DRHomogeneousTreeLikelihood(
      const Tree& tree,
      std::shared_ptr<TransitionModelInterface> model,
      std::shared_ptr<DiscreteDistributionInterface> rDist,
      bool checkRooted = true,
      bool verbose = true);

  /**
   * @brief Build a new DRHomogeneousTreeLikelihood object and compute the corresponding likelihood.
   *
   * This constructor initializes all parameters, data, and likelihood arrays.
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
  DRHomogeneousTreeLikelihood(
      const Tree& tree,
      const AlignmentDataInterface& data,
      std::shared_ptr<TransitionModelInterface> model,
      std::shared_ptr<DiscreteDistributionInterface> rDist,
      bool checkRooted = true,
      bool verbose = true);

  /**
   * @brief Copy constructor.
   */
  DRHomogeneousTreeLikelihood(const DRHomogeneousTreeLikelihood& lik);

  DRHomogeneousTreeLikelihood& operator=(const DRHomogeneousTreeLikelihood& lik);

  virtual ~DRHomogeneousTreeLikelihood() {}

  DRHomogeneousTreeLikelihood* clone() const override { return new DRHomogeneousTreeLikelihood(*this); }

private:
  /**
   * @brief Method called by constructors.
   */
  void init_();

public:
  /**
   * @name The TreeLikelihood interface.
   *
   * Other methods are implemented in the AbstractTreeLikelihood class.
   *
   * @{
   */
  void setData(const AlignmentDataInterface& sites) override;
  double getLikelihood() const override;
  double getLogLikelihood() const override;
  double getLikelihoodForASite (size_t site) const override;
  double getLogLikelihoodForASite(size_t site) const override;
  size_t getSiteIndex(size_t site) const override { return likelihoodData_->getRootArrayPosition(site); }
  /** @} */

  void computeTreeLikelihood();


  /**
   * @name The DiscreteRatesAcrossSites interface implementation:
   *
   * @{
   */
  double getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const override;
  double getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const override;
  double getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const override;
  double getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const override;
  /** @} */

  /**
   * @brief Implements the Function interface.
   *
   * Update the parameter list and call the fireParameterChanged() method.
   *
   * If a subset of the whole parameter list is passed to the function,
   * only these parameters are updated and the other remain constant (i.e.
   * equal to their last value).
   *
   * @param parameters The parameter list to pass to the function.
   */
  void setParameters(const ParameterList& parameters) override;

  /**
   * @brief Function and NNISearchable interface.
   */
  double getValue() const override;

  /**
   * @name DerivableFirstOrder interface.
   *
   * @{
   */
  double getFirstOrderDerivative(const std::string& variable) const override;
  /** @{ */

  /**
   * @name DerivableSecondOrder interface.
   *
   * @{
   */
  double getSecondOrderDerivative(const std::string& variable) const override;
  double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const override { return 0; } // Not implemented for now.
  /** @} */

public:
  // Specific methods:

  DRASDRTreeLikelihoodData& likelihoodData() override { return *likelihoodData_; }
  const DRASDRTreeLikelihoodData& likelihoodData() const override { return *likelihoodData_; }

  virtual void computeLikelihoodAtNode(int nodeId, VVVdouble& likelihoodArray) const override
  {
    computeLikelihoodAtNode_(tree_->getNode(nodeId), likelihoodArray);
  }

protected:
  virtual void computeLikelihoodAtNode_(const Node* node, VVVdouble& likelihoodArray, const Node* sonNode = 0) const;

  /**
   * Initialize the arrays corresponding to each son node for the node passed as argument.
   * The method is called for each son node and the result stored in the corresponding array.
   */
  virtual void computeSubtreeLikelihoodPostfix(const Node* node); // Recursive method.
  /**
   * This method initilize the remaining likelihood arrays, corresponding to father nodes.
   * It must be called after the postfix method because it requires that the arrays for
   * son nodes to be be computed.
   */
  virtual void computeSubtreeLikelihoodPrefix(const Node* node); // Recursive method.

  virtual void computeRootLikelihood();

  virtual void computeTreeDLikelihoodAtNode(const Node* node);
  virtual void computeTreeDLikelihoods();

  virtual void computeTreeD2LikelihoodAtNode(const Node* node);
  virtual void computeTreeD2Likelihoods();

  virtual void fireParameterChanged(const ParameterList& params) override;

  virtual void resetLikelihoodArrays(const Node* node);

  /**
   * @brief This method is mainly for debugging purpose.
   *
   * @param node The node at which likelihood values must be displayed.
   */
  virtual void displayLikelihood(const Node* node);

  /**
   * @brief Compute conditional likelihoods.
   *
   * This method is the "core" likelihood computation function, performing all the product uppon all nodes, the summation for each ancestral state and each rate class.
   * It is designed for inner usage, and a maximum efficiency, so no checking is performed on the input parameters.
   * Use with care!
   *
   * @param iLik A vector of likelihood arrays, one for each conditional node.
   * @param tProb A vector of transition probabilities, one for each node.
   * @param oLik The likelihood array to store the computed likelihoods.
   * @param nbNodes The number of nodes = the size of the input vectors.
   * @param nbDistinctSites The number of distinct sites (the first dimension of the likelihood array).
   * @param nbClasses The number of rate classes (the second dimension of the likelihood array).
   * @param nbStates The number of states (the third dimension of the likelihood array).
   * @param reset Tell if the output likelihood array must be initalized prior to computation.
   * If true, the resetLikelihoodArray method will be called.
   */
  static void computeLikelihoodFromArrays(
    const std::vector<const VVVdouble*>& iLik,
    const std::vector<const VVVdouble*>& tProb,
    VVVdouble& oLik, size_t nbNodes,
    size_t nbDistinctSites,
    size_t nbClasses,
    size_t nbStates,
    bool reset = true);

  /**
   * @brief Compute conditional likelihoods.
   *
   * This method is the "core" likelihood computation function, performing all the product uppon all nodes, the summation for each ancestral state and each rate class.
   * This function is specific to non-reversible models: the subtree containing the root is specified separately.
   * It is designed for inner usage, and a maximum efficiency, so no checking is performed on the input parameters.
   * Use with care!
   *
   * @param iLik A vector of likelihood arrays, one for each conditional node.
   * @param tProb A vector of transition probabilities, one for each node.
   * @param iLikR The likelihood array for the subtree containing the root of the tree.
   * @param tProbR The transition probabilities for thr subtree containing the root of the tree.
   * @param oLik The likelihood array to store the computed likelihoods.
   * @param nbNodes The number of nodes = the size of the input vectors.
   * @param nbDistinctSites The number of distinct sites (the first dimension of the likelihood array).
   * @param nbClasses The number of rate classes (the second dimension of the likelihood array).
   * @param nbStates The number of states (the third dimension of the likelihood array).
   * @param reset Tell if the output likelihood array must be initalized prior to computation.
   * If true, the resetLikelihoodArray method will be called.
   */
  static void computeLikelihoodFromArrays(
    const std::vector<const VVVdouble*>& iLik,
    const std::vector<const VVVdouble*>& tProb,
    const VVVdouble* iLikR,
    const VVVdouble* tProbR,
    VVVdouble& oLik,
    size_t nbNodes,
    size_t nbDistinctSites,
    size_t nbClasses,
    size_t nbStates,
    bool reset = true);

  friend class DRHomogeneousMixedTreeLikelihood;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_DRHOMOGENEOUSTREELIKELIHOOD_H
