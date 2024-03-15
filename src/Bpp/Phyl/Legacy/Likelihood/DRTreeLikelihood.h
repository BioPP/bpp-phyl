// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_DRTREELIKELIHOOD_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_DRTREELIKELIHOOD_H

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/VectorTools.h>

#include "AbstractNonHomogeneousTreeLikelihood.h"
#include "DRASDRTreeLikelihoodData.h"

namespace bpp
{
/**
 * @brief Interface for double-recursive (DR) implementation of the likelihood computation.
 *
 * In the DR implementation, three conditional likelihoods are stored at each node, corresponding to the three connected substrees
 * (in case of multifurcations, there may have even more).
 * The DR implementation hence uses 3x more memory than the simple recursive (R) implementation.
 * However, the likelihood of the tree can be computed independently at each node of the tree,
 * which is very convenient for topology estimation (the likelihood change of each NNI movement can be computed directly),
 * ancestral state reconstruction (all marginal ancestral states can be computed in one pass over the tree), or for
 * substitution mapping.
 *
 * This interface provides
 * - a method to access the DR likelihood data structure,
 * - a method to compute the likelihood array at each node.
 *
 * For now, this interface inherits from DiscreteRatesAcrossSitesTreeLikelihood and not TreeLikelihood,
 * since the data structure available accounts for rate across site variation.
 * This may change in the future.
 *
 * @see DRTreeLikelihoodTools
 */
class DRTreeLikelihoodInterface :
  public virtual DiscreteRatesAcrossSitesTreeLikelihoodInterface
{
public:
  DRTreeLikelihoodInterface() {}
  virtual ~DRTreeLikelihoodInterface() {}

  DRTreeLikelihoodInterface* clone() const override = 0;

public:

public:
  /**
   * @name Get the likelihood data structure associated to this class.
   *
   * @{
   */
  virtual DRASDRTreeLikelihoodData& likelihoodData() override = 0;
  virtual const DRASDRTreeLikelihoodData& likelihoodData() const override = 0;
  /** @} */

  /**
   * @brief Compute the likelihood array at a given node.
   *
   * @param nodeId The id of the node to consider.
   * @param likelihoodArray The array where to store the results.
   */
  virtual void computeLikelihoodAtNode(int nodeId, VVVdouble& likelihoodArray) const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_DRTREELIKELIHOOD_H
