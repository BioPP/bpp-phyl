//
// File: RewardMappingToolsForASite.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: mercredi 26 septembre 2018, ÃÂ  17h 44
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef BPP_PHYL_MAPPING_REWARDMAPPINGTOOLSFORASITE_H
#define BPP_PHYL_MAPPING_REWARDMAPPINGTOOLSFORASITE_H


#include "../Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h"
#include "ProbabilisticRewardMapping.h"
#include "Reward.h"

namespace bpp
{
/**
 * @brief Provide methods to compute reward mappings.
 *
 * See:
 * Minin, V.N. and Suchard, M.A.,
 * Fast, accurate and simulation-free stochastic mapping
 * Philosophical Transactions of the Royal Society B 2008 363:3985-95.
 *
 * @author Laurent GuÃÂ©guen
 */
class RewardMappingToolsForASite
{
public:
  RewardMappingToolsForASite() {}
  virtual ~RewardMappingToolsForASite() {}

public:
  /**
   * @brief Compute the reward vectors for a particular dataset
   * using the double-recursive likelihood computation.
   *
   * @param site              The site
   * @param rltc              A LikelihoodCalculationSingleProcess object.
   * @param nodeIds           The Ids of the nodes the reward vectors
   *                          are computed on.
   * @param reward            The Reward to use.
   * @param verbose           Print info to screen.
   * @return A vector of reward vectors (one for each site).
   * @throw Exception If the likelihood object is not initialized.
   */

  static ProbabilisticRewardMapping* computeRewardVectors(
    size_t site,
    LikelihoodCalculationSingleProcess& rltc,
    const std::vector<uint>& nodeIds,
    Reward& reward,
    bool verbose = true);
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_REWARDMAPPINGTOOLSFORASITE_H
