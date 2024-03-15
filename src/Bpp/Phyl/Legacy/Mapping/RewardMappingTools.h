// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _REWARDMAPPINGTOOLS_H_
#define _REWARDMAPPINGTOOLS_H_

#include "ProbabilisticRewardMapping.h"
#include "../../Mapping/Reward.h"

#include "../Likelihood/DRTreeLikelihood.h"

namespace bpp
{
/**
 * @brief Provide methods to compute reward mappings.
 *
 * For now, 4 methods are implemented, and provide reward mappings.
 *
 * See:
 * Minin, V.N. and Suchard, M.A.,
 * Fast, accurate and simulation-free stochastic mapping
 * Philosophical Transactions of the Royal Society B 2008 363:3985-95.
 *
 * @author Laurent Guéguen
 */
class LegacyRewardMappingTools
{
public:
  LegacyRewardMappingTools() {}
  virtual ~LegacyRewardMappingTools() {}

public:
  /**
   * @brief Compute the reward vectors for a particular dataset
   * using the double-recursive likelihood computation.
   *
   * @param drtl              A DRTreeLikelihood object.
   * @param nodeIds           The Ids of the nodes the reward vectors
   *                          are computed on.
   * @param reward            The Reward to use.
   * @param verbose           Print info to screen.
   * @return A vector of reward vectors (one for each site).
   * @throw Exception If the likelihood object is not initialized.
   */
  static std::unique_ptr<LegacyProbabilisticRewardMapping> computeRewardVectors(
    std::shared_ptr<const DRTreeLikelihoodInterface> drtl,
    const std::vector<int>& nodeIds,
    std::shared_ptr<Reward> reward,
    bool verbose = true);


  /**
   * @brief Write the reward vectors to a stream.
   *
   * @param rewards The reward vectors to write.
   * @param sites         The dataset associated to the vectors
   * (needed to know the position of each site in the dataset).
   * @param out           The output stream where to write the vectors.
   * @throw IOException If an output error happens.
   */
  static void writeToStream(
    const LegacyProbabilisticRewardMapping& rewards,
    const SiteContainerInterface& sites,
    std::ostream& out);


  /**
   * @brief Read the reward vectors from a stream.
   *
   * @param in            The input stream where to read the vectors.
   * @param rewards       The mapping object to fill.
   * @throw IOException If an input error happens.
   */
  static void readFromStream(std::istream& in, LegacyProbabilisticRewardMapping& rewards);


  /**
   * @brief Sum all rewards of a given branch (specified by its index).
   *
   * @param smap The reward map to use.
   * @param branchIndex The index of the reward vector for which the counts should be computed.
   * @return A vector will all rewards summed.
   */
  static double computeSumForBranch(const LegacyRewardMappingInterface& smap, size_t branchIndex);


  /**
   * @brief Sum all substitutions for each type of a given site (specified by its index).
   *
   * @param smap The substitution map to use.
   * @param siteIndex The index of the substitution vector for which the counts should be computed.
   * @return A vector will all counts summed for each types of substitutions.
   */
  static double computeSumForSite(const LegacyRewardMappingInterface& smap, size_t siteIndex);
};
} // end of namespace bpp.

#endif // _REWARDMAPPINGTOOLS_H_
