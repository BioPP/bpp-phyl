// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_SITEPARTITIONTREELIKELIHOOD_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_SITEPARTITIONTREELIKELIHOOD_H


#include "TreeLikelihood.h"

namespace bpp
{
/**
 * @brief Specialization of the TreeLikelihood interface for partition models, homogeneous case.
 *
 * These models allow the distinct sites of an alignment to have a different model.
 * The substitution model is however assumed to be the same along the tree.
 * Such models are hence homogeneous in time.
 */
class SitePartitionHomogeneousTreeLikelihood :
  public virtual TreeLikelihood
{
public:
  SitePartitionHomogeneousTreeLikelihood* clone() const = 0;

public:
  const TransitionModel* getModel(int nodeId, size_t siteIndex) const
  {
    return getModelForSite(siteIndex);
  }

  TransitionModel* getModel(int nodeId, size_t siteIndex)
  {
    return getModelForSite(siteIndex);
  }

  /**
   * @brief Get the substitution model associated to a given node.
   *
   * @param siteIndex The position in the alignment.
   * @return A pointer toward the corresponding model.
   */
  virtual const TransitionModel* getModelForSite(size_t siteIndex) const = 0;

  /**
   * @brief Get the substitution model associated to a given node.
   *
   * @param siteIndex The position in the alignment.
   * @return A pointer toward the corresponding model.
   */
  virtual TransitionModel* getModelForSite(size_t siteIndex) = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_SITEPARTITIONTREELIKELIHOOD_H
