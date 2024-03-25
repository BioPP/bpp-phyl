// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_RHOMOGENEOUSMIXEDTREELIKELIHOOD_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_RHOMOGENEOUSMIXEDTREELIKELIHOOD_H

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/VectorTools.h>

#include "../../Model/MixedTransitionModel.h"
#include "../../Model/SubstitutionModel.h"
#include "RHomogeneousTreeLikelihood.h"

namespace bpp
{
/**
 *@ brief A class to compute the average of several
 * RHomogeneousTreeLikelihood defined from a Mixed Substitution
 * Model.
 *
 * In all the calculus, the average of the likelihoods, probabilities
 * are computed.
 **/

class RHomogeneousMixedTreeLikelihood :
  public RHomogeneousTreeLikelihood
{
private:
  std::vector< std::shared_ptr<RHomogeneousTreeLikelihood>> treeLikelihoodsContainer_;
  std::vector<double> probas_;

public:
  /**
   * @brief Build a new RHomogeneousMixedTreeLikelihood object without
   * data.
   *
   * This constructor only initialize the parameters. To compute a
   * likelihood, you will need to call the setData() and the
   * computeTreeLikelihood() methods.
   *
   * @param tree The tree to use.
   * @param model The mixed substitution model to use.
   * @param rDist The rate across sites distribution to use.
   * @param checkRooted Tell if we have to check for the tree to be unrooted.
   * If true, any rooted tree will be unrooted before likelihood computation.
   * @param verbose Should I display some info?
   * @param usePatterns Tell if recursive site compression should be performed.
   * @throw Exception in an error occurred.
   */
  RHomogeneousMixedTreeLikelihood(
      const Tree& tree,
      std::shared_ptr<TransitionModelInterface> model,
      std::shared_ptr<DiscreteDistributionInterface> rDist,
      bool checkRooted = true,
      bool verbose = true,
      bool usePatterns = true);

  /**
   * @brief Build a new RHomogeneousMixedTreeLikelihood object with data.
   *
   * This constructor initializes all parameters, data, and likelihood arrays.
   *
   * @param tree The tree to use.
   * @param data Sequences to use.
   * @param model The mixed substitution model to use.
   * @param rDist The rate across sites distribution to use.
   * @param checkRooted Tell if we have to check for the tree to be unrooted.
   * If true, any rooted tree will be unrooted before likelihood computation.
   * @param verbose Should I display some info?
   * @param usePatterns Tell if recursive site compression should be performed.
   * @throw Exception in an error occurred.
   */
  RHomogeneousMixedTreeLikelihood(
      const Tree& tree,
      const AlignmentDataInterface& data,
      std::shared_ptr<TransitionModelInterface> model,
      std::shared_ptr<DiscreteDistributionInterface> rDist,
      bool checkRooted = true,
      bool verbose = true,
      bool usePatterns = true);

  RHomogeneousMixedTreeLikelihood(const RHomogeneousMixedTreeLikelihood& lik);

  RHomogeneousMixedTreeLikelihood& operator=(const RHomogeneousMixedTreeLikelihood& lik);

  virtual ~RHomogeneousMixedTreeLikelihood();

  RHomogeneousMixedTreeLikelihood* clone() const { return new RHomogeneousMixedTreeLikelihood(*this); }

public:
  /**
   * @name The TreeLikelihood interface.
   *
   * Other methods are implemented in the RHomogeneousTreeLikelihood class.
   *
   * @{
   */
  void setData(const AlignmentDataInterface& sites);

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

public:
  // Specific methods:
  void initialize();

  void fireParameterChanged(const ParameterList& params);

  void computeTreeLikelihood();

  virtual double getDLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

  virtual void computeTreeDLikelihood(const std::string& variable);

  virtual double getD2LikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;

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

  /**
   * @brief This method is used by fireParameterChanged method.
   *
   */
  void computeAllTransitionProbabilities();
  /**
   * @brief This method is used by fireParameterChanged method.
   *
   */
  void computeTransitionProbabilitiesForNode(const Node* node);

  /**
   * @brief This method is mainly for debugging purpose.
   *
   * @param node The node at which likelihood values must be displayed.
   */
  virtual void displayLikelihood(const Node* node);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_RHOMOGENEOUSMIXEDTREELIKELIHOOD_H
