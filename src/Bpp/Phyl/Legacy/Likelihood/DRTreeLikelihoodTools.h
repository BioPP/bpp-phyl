// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_DRTREELIKELIHOODTOOLS_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_DRTREELIKELIHOODTOOLS_H

#include <Bpp/Seq/Container/AlignedSequenceContainer.h>

#include "DRTreeLikelihood.h"
#include "TreeLikelihoodTools.h"

namespace bpp
{
/**
 * @brief Utilitary methods dealing with objects implementing the DRTreeLikelihood interface.
 */
class DRTreeLikelihoodTools :
  public TreeLikelihoodTools
{
public:
  /**
   * @brief Compute the posterior probabilities for each state and each rate of each distinct site.
   *
   * @param drl A DR tree likelihood object.
   * @param nodeId The id of the node at which probabilities must be computed.
   * @return A 3-dimensional array, with probabilities for each site, each rate and each state.
   */
  static VVVdouble getPosteriorProbabilitiesPerStatePerRate(
      const DRTreeLikelihoodInterface& drl,
      int nodeId);

  /**
   * @brief Compute the posterior probabilities for each state for a given node.
   *
   * This method calls the getPosteriorProbabilitiesPerStatePerRate function
   * and average the probabilities over all sites and rate classes, resulting in a
   * one-dimensionnal frequency array, with one frequency per model state.
   *
   * @param drl A DR tree likelihood object.
   * @param nodeId The id of the node at which probabilities must be computed.
   * @return vector of double with state frequencies for the given node.
   */
  static Vdouble getPosteriorStateFrequencies(
      const DRTreeLikelihoodInterface& drl,
      int nodeId);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_DRTREELIKELIHOODTOOLS_H
