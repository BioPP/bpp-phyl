// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _LEGACY_REWARD_MAPPING_H_
#define _LEGACY_REWARD_MAPPING_H_

#include "Mapping.h"

#include <Bpp/Clonable.h>

//From the STL:
#include <vector>
#include <memory>

namespace bpp
{

/**
 * @brief Legacy interface for storing reward mapping data.
 *
 * Since only probabilistic reward mapping is implemented for now, the basal 
 * interface only contains a few methods.
 * More methods are expected to be added later.
 */
  
  class LegacyRewardMappingInterface:
    virtual public LegacyMappingInterface
  {
    
  public:
    LegacyRewardMappingInterface() {}
    virtual ~LegacyRewardMappingInterface() {}

    LegacyRewardMappingInterface* clone() const override = 0;

  public:
    
    virtual double& operator()(size_t nodeIndex, size_t siteIndex) = 0;
    virtual const double& operator()(size_t nodeIndex, size_t siteIndex) const = 0;
};



/**
 * @brief Partial implementation of the reward mapping interface.
 *
 * This implementation copies the input tree in a TreeTemplate<Node> object.
 */
class LegacyAbstractRewardMapping:
    public virtual LegacyRewardMappingInterface,
    public LegacyAbstractMapping
{
public:
  LegacyAbstractRewardMapping(const Tree& tree) : LegacyAbstractMapping(tree){}

  LegacyAbstractRewardMapping(const LegacyAbstractRewardMapping& arm) = default;

  LegacyAbstractRewardMapping* clone() const override = 0;

  LegacyAbstractRewardMapping& operator=(const LegacyAbstractRewardMapping& arm) = default;

  virtual ~LegacyAbstractRewardMapping() {}

};

} //end of namespace bpp.

#endif //_LEGACY_REWARD_MAPPING_H_

