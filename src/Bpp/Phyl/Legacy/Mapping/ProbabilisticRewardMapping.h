// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _LEGACY_PROBABILISTIC_REWARD_MAPPING_H_
#define _LEGACY_PROBABILISTIC_REWARD_MAPPING_H_

#include "RewardMapping.h"
#include "../../Mapping/Reward.h" //We use the new implementation here
#include "../../Tree/TreeExceptions.h"

#include <Bpp/Text/TextTools.h>

// From the STL:
#include <vector>

namespace bpp
{

/**
 * @brief Legacy data storage class for probabilistic rewards mappings.
 *
 * A 'probabilistic' mapping contains an single value for each branch and each site.
 * This number is an average reward.
 * Probabilistic was coined there by opposition to the'stochastic'
 * mapping, where a path (sequence of rewards along the branch) is
 * available for each branch and site.
 */
  class LegacyProbabilisticRewardMapping:
    public LegacyAbstractRewardMapping
  {
  private:
    std::shared_ptr<const Reward> reward_;

    /**
     * @brief Rewards storage.
     *
     * Rewards are stored by sites.
     */
    std::vector<std::vector<double>> mapping_;
    
  public:
    
    /**
     * @brief Build a new LegacyProbabilisticRewardMapping object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     * @param reward A pointer toward the Reward object that has been used for the mapping, if any.
     * @param numberOfSites The number of sites to map.
     */
    LegacyProbabilisticRewardMapping(const Tree& tree, std::shared_ptr<const Reward> reward, size_t numberOfSites) :
      LegacyAbstractRewardMapping(tree),
      reward_(reward),
      mapping_(0)
    {
      setNumberOfSites(numberOfSites);
    }
    
    /**
     * @brief Build a new ProbabilisticRewardMapping object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     */
    LegacyProbabilisticRewardMapping(const Tree& tree) :
      LegacyAbstractRewardMapping(tree),
      reward_(nullptr),
      mapping_(0)
    {}
    

    LegacyProbabilisticRewardMapping* clone() const override
    { 
      return new LegacyProbabilisticRewardMapping(*this);
    }

    LegacyProbabilisticRewardMapping(const LegacyProbabilisticRewardMapping& prm) = default;

    LegacyProbabilisticRewardMapping& operator=(const LegacyProbabilisticRewardMapping& prm) = default;

    virtual ~LegacyProbabilisticRewardMapping() {}

  public:

    virtual double getReward(int nodeId, size_t siteIndex) const
    {
      return mapping_[siteIndex][getNodeIndex(nodeId)];
    }
    
    /**
     * @brief (Re)-set the phylogenetic tree associated to this mapping.
     *
     * @param tree The new tree.
     */
    virtual void setTree(const Tree& tree);

    virtual void setNumberOfSites(size_t numberOfSites) override;
    
    /**
     * @brief Direct access to rewards.
     *
     * @warning No index checking is performed, use with care!
     */
    virtual double& operator()(size_t nodeIndex, size_t siteIndex) override
    {
      return mapping_[siteIndex][nodeIndex];
    }

    /**
     * @brief Direct access to rewards.
     *
     * @warning No index checking is performed, use with care!
     */
    virtual const double& operator()(size_t nodeIndex, size_t siteIndex) const override
    {
      return mapping_[siteIndex][nodeIndex];
    }
     
    /**
     * @brief Direct access to rewards.
     *
     * @warning No index checking is performed, use with care!
     */
    std::vector<double>& operator[](size_t siteIndex)
    {
      return mapping_[siteIndex];
    }

    /**
     * @brief Direct access to rewards.
     *
     * @warning No index checking is performed, use with care!
     */
    const std::vector<double>& operator[](size_t siteIndex) const
    {
      return mapping_[siteIndex];
    }
};

} //end of namespace bpp.

#endif //_LEGACY_PROBABILISTIC_REWARD_MAPPING_H_

