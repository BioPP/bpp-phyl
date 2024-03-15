// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_HOMOGENEOUSTREELIKELIHOOD_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_HOMOGENEOUSTREELIKELIHOOD_H


#include "../../Model/SubstitutionModel.h"
#include "TreeLikelihood.h"

namespace bpp
{
/**
 * @brief Specialization of the TreeLikelihood interface for the Homogeneous case.
 *
 * Homogeneous models assume a unique substitution model along the tree.
 * This interface further assumes that  the substitution model is the same for all sites.
 * For likelihood functions with different model per sites, see SitePartitionHomogeneousTreeLikelihood.
 *
 * @see SubstitutionModel, SitePartitionHomogeneousTreeLikelihood.
 */
class HomogeneousTreeLikelihood :
  public virtual TreeLikelihoodInterface
{
public:
  HomogeneousTreeLikelihood* clone() const = 0;

public:
  std::shared_ptr<const TransitionModelInterface> getModelForSite(int nodeId, size_t siteIndex) const
  {
    return getModel();
  }

  std::shared_ptr<TransitionModelInterface> getModelForSite(int nodeId, size_t siteIndex)
  {
    return getModel();
  }

  /**
   * @return The substitution model attached to this instance.
   */
  virtual std::shared_ptr<const TransitionModelInterface> getModel() const = 0;

  /**
   * @return The substitution model attached to this instance.
   */
  virtual std::shared_ptr<TransitionModelInterface> getModel() = 0;

  /**
   * @return Set the substitution model for this instance.
   * @throw Exception If the model could not be set (for instance, because of a wrong alphabet type).
   */
  virtual void setModel(std::shared_ptr<TransitionModelInterface> model) = 0;

  /**
   * @brief Get a SubstitutionModel pointer toward the model associated to this instance, if possible.
   *
   * Performs a cast operation on the pointer. Return NULL if cast failed.
   * @return A SubstitutionModel pointer toward the model associated to this instance.
   */
  virtual std::shared_ptr<const SubstitutionModelInterface> getSubstitutionModel() const = 0;

  /**
   * @brief Get a SubstitutionModel pointer toward the model associated to this instance, if possible.
   *
   * Performs a cast operation on the pointer. Return NULL if cast failed.
   * @return A SubstitutionModel pointer toward the model associated to this instance.
   *
   * @param nodeId Id of the node
   * @param siteIndex Position of the site
   */
  virtual std::shared_ptr<const SubstitutionModelInterface> getSubstitutionModel(int nodeId, size_t siteIndex) const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_HOMOGENEOUSTREELIKELIHOOD_H
