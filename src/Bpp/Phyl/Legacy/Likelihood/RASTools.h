// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_RASTOOLS_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_RASTOOLS_H

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

#include "DiscreteRatesAcrossSitesTreeLikelihood.h"

namespace bpp
{
/**
 * @brief Tools to deal with Rates Across Sites (RAS) models.
 */
class RASTools
{
public:
  /**
   * @brief Get the rate distribution estimated from a dataset.
   *
   * This methods takes an objet implementing the DiscreteRatesAcroossSites
   * interface as input and use the posterior probabilities of rates for each site
   * to generate the corresponding distribution.
   *
   * @param treeLikelihood A Likelihood calculation implmenting the RAS model interface.
   * @return The posterior distribution of rate classes.
   */
  static std::unique_ptr<DiscreteDistributionInterface> getPosteriorRateDistribution(
      const DiscreteRatesAcrossSitesTreeLikelihoodInterface& treeLikelihood);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_RASTOOLS_H
