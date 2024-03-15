// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_REWARDMAPPING_H
#define BPP_PHYL_MAPPING_REWARDMAPPING_H

#include <Bpp/Clonable.h>

#include "Mapping.h"

// From the STL:
#include <vector>
#include <memory>

namespace bpp
{
/**
 * @brief General interface for storing reward mapping data.
 *
 * Since only probabilistic reward mapping is implemented for now, the basal
 * interface only contains a few methods.
 * More methods are expected to be added later.
 */
class RewardMappingInterface :
  public virtual MappingInterface
{
public:
  RewardMappingInterface() {}
  virtual ~RewardMappingInterface() {}

  RewardMappingInterface* clone() const override = 0;

public:
  virtual double& operator()(uint branchId, size_t siteIndex) = 0;
  virtual double operator()(uint branchId, size_t siteIndex) const = 0;
};


/**
 * @brief Partial implementation of the substitution mapping interface.
 *
 * This implementation copies the input tree in a TreeTemplate<Node> object.
 */

class AbstractRewardMapping :
  public virtual RewardMappingInterface,
  public AbstractMapping
{
public:
  AbstractRewardMapping() : AbstractMapping() {}

  AbstractRewardMapping(const AbstractRewardMapping& absm) :
    AbstractMapping(absm) {}

  AbstractRewardMapping* clone() const override = 0;

  AbstractRewardMapping& operator=(const AbstractRewardMapping& absm)
  {
    AbstractMapping::operator=(absm);
    return *this;
  }

  virtual ~AbstractRewardMapping() {}
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_REWARDMAPPING_H
