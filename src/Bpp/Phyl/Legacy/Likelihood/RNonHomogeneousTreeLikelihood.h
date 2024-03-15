// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_RNONHOMOGENEOUSTREELIKELIHOOD_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_RNONHOMOGENEOUSTREELIKELIHOOD_H

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/VectorTools.h>

#include "../Model/SubstitutionModelSet.h"
#include "AbstractNonHomogeneousTreeLikelihood.h"
#include "DRASRTreeLikelihoodData.h"

namespace bpp
{
/**
 * @brief This class implement the 'traditional' way of computing likelihood for a tree, allowing for non-homogeneous models of substitutions.
 *
 * A non uniform distribution of rates among the sites is allowed (ASRV models).</p>
 *
 * This class uses an instance of the DRASRTreeLikelihoodData for conditionnal likelihood storage.
 *
 * This class can also use a simple or recursive site compression.
 * In the simple case, computations for identical sites are not duplicated.
 * In the recursive case, computations for identical sub-sites (<i>site patterns </i>) are also not duplicated:
 * Following N. Galtier (personal communication), we define a Pattern as a distinct site
 * in a sub-dataset corresponding to the dataset with sequences associated to a particular subtree.
 * The likelihood computation is the same for a given site, hence the idea is to save time from
 * performing many times the same coputation.
 * The network between all patterns is defined by the _patternLinks double map, initialized in the
 * initLikelihoodsWithPatterns() method. This initialisation takes more time than the classic
 * initTreeLikelihood one, where all likelihoods for a given site <i>i</i> are at the <i>i</i> coordinate
 * in the likelihood tensor, but is really faster when computing the likelihood (computeLikelihoods() method).
 * Hence, if you have to compute likelihood many times while holding the tree topology unchanged,
 * you should use patterns.
 * This decreases the likelihood computation time, but at a cost: some time is spent to establish the patterns
 * relationships. Whether to use or not patterns depends on what you actllay need:
 * - The more you compute likelihoods without changing the data or topology, the more patterns are interesting
 *   (this divides the cost of computing patterns by the number of computation performed).
 *   Patterns are hence usefull when you have a high number of computation to perform, while optimizing numerical
 *   parameters for instance).
 * - Patterns are more likely to occur whith small alphabet (nucleotides).
 *
 * Important note: The input tree will be considered as rooted, since the likelihood of non-stationary models
 * depends on the position of the root. If the input tree is not rooted, it will be considered as a rotted tree
 * with a root multifurcation.
 */
class RNonHomogeneousTreeLikelihood :
  public AbstractNonHomogeneousTreeLikelihood
{
private:
  mutable std::shared_ptr<DRASRTreeLikelihoodData> likelihoodData_;
  double minusLogLik_;

public:
  /**
   * @brief Build a new NonHomogeneousTreeLikelihood object without data.
   *
   * This constructor only initialize the parameters.
   * To compute a likelihood, you will need to call the setData() and the computeTreeLikelihood() methods.
   *
   * @param tree The tree to use.
   * @param modelSet The set of substitution models to use.
   * @param rDist The rate across sites distribution to use.
   * If true, any rooted tree will be unrooted before likelihood computation.
   * @param verbose Should I display some info?
   * @param usePatterns Tell if recursive site compression should be performed.
   * @param reparametrizeRoot Should we reparametrize the branch lengths at root?
   * @throw Exception in an error occured.
   */
  RNonHomogeneousTreeLikelihood(
    const Tree& tree,
    std::shared_ptr<SubstitutionModelSet> modelSet,
    std::shared_ptr<DiscreteDistributionInterface> rDist,
    bool verbose = true,
    bool usePatterns = true,
    bool reparametrizeRoot = false);

  /**
   * @brief Build a new NonHomogeneousTreeLikelihood object and compute the corresponding likelihood.
   *
   * This constructor initializes all parameters, data, and likelihood arrays.
   *
   * @param tree The tree to use.
   * @param data Sequences to use.
   * @param modelSet The set of substitution models to use.
   * @param rDist The rate across sites distribution to use.
   * @param verbose Should I display some info?
   * @param usePatterns Tell if recursive site compression should be performed.
   * @param reparametrizeRoot Should we reparametrize the branch lengths at root?
   * @throw Exception in an error occured.
   */
  RNonHomogeneousTreeLikelihood(
    const Tree& tree,
    const AlignmentDataInterface& data,
    std::shared_ptr<SubstitutionModelSet> modelSet,
    std::shared_ptr<DiscreteDistributionInterface> rDist,
    bool verbose = true,
    bool usePatterns = true,
    bool reparametrizeRoot = false);

  RNonHomogeneousTreeLikelihood(const RNonHomogeneousTreeLikelihood& lik);

  RNonHomogeneousTreeLikelihood& operator=(const RNonHomogeneousTreeLikelihood& lik);

  virtual ~RNonHomogeneousTreeLikelihood();

  RNonHomogeneousTreeLikelihood* clone() const { return new RNonHomogeneousTreeLikelihood(*this); }

private:
  /**
   * @brief Method called by constructors.
   */
  void init_(bool usePatterns);

public:
  /**
   * @name The TreeLikelihood interface.
   *
   * Other methods are implemented in the AbstractHomogeneousTreeLikelihood class.
   *
   * @{
   */
  void setData(const AlignmentDataInterface& sites);
  double getLikelihood() const;
  double getLogLikelihood() const;
  double getLikelihoodForASite(size_t site) const;
  double getLogLikelihoodForASite(size_t site) const;
  size_t getSiteIndex(size_t site) const { return likelihoodData_->getRootArrayPosition(site); }
  /** @} */


  /**
   * @name The DiscreteRatesAcrossSites interface implementation:
   *
   * @{
   */
  double getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;
  double getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;
  double getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;
  double getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;
  /** @} */

  /**
   * @brief Implements the Function interface.
   *
   * Update the parameter list and call the applyParameters() method.
   * Then compute the likelihoods at each node (computeLikelihood() method)
   * and call the getLogLikelihood() method.
   *
   * If a subset of the whole parameter list is passed to the function,
   * only these parameters are updated and the other remain constant (i.e.
   * equal to their last value).
   *
   * @param parameters The parameter list to pass to the function.
   */
  void setParameters(const ParameterList& parameters);
  double getValue() const;

  /**
   * @name DerivableFirstOrder interface.
   *
   * @{
   */
  double getFirstOrderDerivative(const std::string& variable) const;
  /** @} */

  /**
   * @name DerivableSecondOrder interface.
   *
   * @{
   */
  double getSecondOrderDerivative(const std::string& variable) const;
  double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const { return 0; } // Not implemented for now.
  /** @} */

public:
  // Specific methods:

  DRASRTreeLikelihoodData& likelihoodData() { return *likelihoodData_; }

  const DRASRTreeLikelihoodData& likelihoodData() const { return *likelihoodData_; }

  virtual void computeTreeLikelihood();

  virtual double getDLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

  virtual double getDLikelihoodForASite(size_t site) const;

  virtual double getDLogLikelihoodForASite(size_t site) const;

  virtual double getDLogLikelihood() const;

  virtual void computeTreeDLikelihood(const std::string& variable);

  virtual double getD2LikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

  virtual double getD2LikelihoodForASite(size_t site) const;

  virtual double getD2LogLikelihoodForASite(size_t site) const;

  virtual double getD2LogLikelihood() const;

  virtual void computeTreeD2Likelihood(const std::string& variable);

protected:
  /**
   * @brief Compute the likelihood for a subtree defined by the Tree::Node <i>node</i>.
   *
   * @param node The root of the subtree.
   */
  virtual void computeSubtreeLikelihood(const Node* node); // Recursive method.

  virtual void computeDownSubtreeDLikelihood(const Node*);

  virtual void computeDownSubtreeD2Likelihood(const Node*);

  void fireParameterChanged(const ParameterList& params);

  /**
   * @brief This method is mainly for debugging purpose.
   *
   * @param node The node at which likelihood values must be displayed.
   */
  virtual void displayLikelihood(const Node* node);


  friend class RNonHomogeneousMixedTreeLikelihood;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_RNONHOMOGENEOUSTREELIKELIHOOD_H
