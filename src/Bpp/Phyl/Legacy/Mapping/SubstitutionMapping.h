// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _LEGACY_SUBSTITUTION_MAPPING_H_
#define _LEGACY_SUBSTITUTION_MAPPING_H_

#include "Mapping.h"

#include <Bpp/Clonable.h>

// From the STL:
#include <vector>
#include <memory>

namespace bpp
{
/**
 * @brief Legacy interface for storing mapping data.
 *
 * There are several kinds of mapping:
 * - Exact mapping, storing the positions of each substitution onto each branch,
 * - Probabilistic mapping, storing the number of substitutions onto each branch.
 *
 * Since only probabilistic substitution mapping is implemented for now, the basal
 * interface only contains a few methods.
 * More methods are expected to be added later.
 */
class LegacySubstitutionMappingInterface :
  public virtual LegacyMappingInterface
{
public:
  LegacySubstitutionMappingInterface() {}
  virtual ~LegacySubstitutionMappingInterface() {}

  LegacySubstitutionMappingInterface* clone() const override = 0;

public:
  /**
   * @return The number of distinct types of substitutions mapped.
   */
  virtual size_t getNumberOfSubstitutionTypes() const = 0;

  virtual double& operator()(size_t nodeIndex, size_t siteIndex, size_t type) = 0;
  virtual const double& operator()(size_t nodeIndex, size_t siteIndex, size_t type) const = 0;
};


/**
 * @brief Partial implementation of the substitution mapping interface.
 *
 * This implementation copies the input tree in a TreeTemplate<Node> object.
 */
class LegacyAbstractSubstitutionMapping :
  public virtual LegacySubstitutionMappingInterface,
  public LegacyAbstractMapping
{
public:
  LegacyAbstractSubstitutionMapping(const Tree& tree) : LegacyAbstractMapping(tree){}

  LegacyAbstractSubstitutionMapping(const LegacyAbstractSubstitutionMapping& absm) = default;

  LegacyAbstractSubstitutionMapping* clone() const override = 0;

  LegacyAbstractSubstitutionMapping& operator=(const LegacyAbstractSubstitutionMapping& absm) = default;

  virtual ~LegacyAbstractSubstitutionMapping() {}
};
} // end of namespace bpp.

#endif // _LEGACY_SUBSTITUTION_MAPPING_H_
