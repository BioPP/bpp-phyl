// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_NONHOMOGENEOUSTREELIKELIHOOD_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_NONHOMOGENEOUSTREELIKELIHOOD_H


#include "../Model/SubstitutionModelSet.h"
#include "TreeLikelihood.h"

namespace bpp
{
/**
 * @brief Specialization of the TreeLikelihood interface for the branch non-homogeneous and non-stationary models.
 *
 * The main difference is that the likelihood depends on the position of the root.
 * The frequencies at the root of the tree are new parameters, which are not necessarily equal to the equilibrium frequencies of the substitution model.
 *
 * This interface further assumes that the substitution model is the same for all sites, for a given branch.
 *
 * @see SubstitutionModelSet.
 */

class NonHomogeneousTreeLikelihood :
  public virtual TreeLikelihoodInterface
{
public:
  NonHomogeneousTreeLikelihood* clone() const override = 0;

public:
  std::shared_ptr<const TransitionModelInterface> getModelForSite(int nodeId, size_t siteIndex) const override
  {
    return getModelForNode(nodeId);
  }

  std::shared_ptr<TransitionModelInterface> getModelForSite(int nodeId, size_t siteIndex) override
  {
    return getModelForNode(nodeId);
  }

  /**
   * @brief Get the model associated to a given node.
   *
   * @param nodeId The id of the request node.
   * @return A pointer toward the corresponding model.
   * @throw NodeNotFoundException This exception may be thrown if the node is not found (depending on the implementation).
   */
  virtual std::shared_ptr<const TransitionModelInterface> getModelForNode(int nodeId) const = 0;

  /**
   * @brief Get the substitution model associated to a given node.
   *
   * @param nodeId The id of the request node.
   * @return A pointer toward the corresponding model.
   * @throw NodeNotFoundException This exception may be thrown if the node is not found (depending on the implementation).
   */
  virtual std::shared_ptr<TransitionModelInterface> getModelForNode(int nodeId) = 0;

  /**
   * @return The set of substitution models associated to this instance.
   */
  virtual std::shared_ptr<const SubstitutionModelSet> getSubstitutionModelSet() const = 0;

  /**
   * @return The set of substitution models associated to this instance.
   */
  virtual const SubstitutionModelSet& substitutionModelSet() const = 0;

  /**
   * @return The set of substitution models associated to this instance.
   */
  virtual std::shared_ptr<SubstitutionModelSet> getSubstitutionModelSet() = 0;

  /**
   * @return The set of substitution models associated to this instance.
   */
  virtual SubstitutionModelSet& substitutionModelSet() = 0;

  /**
   * @brief Set the substitution models for this instance.
   * @throw Exception If the model could not be set (for instance, because of a wrong alphabet type).
   */
  virtual void setSubstitutionModelSet(std::shared_ptr<SubstitutionModelSet> model) = 0;

  /**
   * @return The parameters on which the root frequencies depend.
   */
  virtual ParameterList getRootFrequenciesParameters() const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_NONHOMOGENEOUSTREELIKELIHOOD_H
