// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_SUBSTITUTIONMAPPING_H
#define BPP_PHYL_MAPPING_SUBSTITUTIONMAPPING_H

#include <Bpp/Clonable.h>

#include "Mapping.h"

// From the STL:
#include <vector>
#include <memory>

namespace bpp
{
/**
 * @brief General interface for storing mapping data.
 *
 * There are several kinds of mapping:
 * - Exact mapping, storing the positions of each substitution onto each branch,
 * - Probabilistic mapping, storing the number of substitutions onto each branch.
 *
 * Since only probabilistic substitution mapping is implemented for now, the basal
 * interface only contains a few methods.
 * More methods are expected to be added later.
 */
class SubstitutionMapping :
  public virtual MappingInterface
{
public:
  SubstitutionMapping() {}
  virtual ~SubstitutionMapping() {}

  SubstitutionMapping* clone() const override = 0;

public:
  /**
   * @return The number of distinct types of substitutions mapped.
   */
  virtual size_t getNumberOfSubstitutionTypes() const = 0;
  virtual void setNumberOfSubstitutionTypes(size_t numberOfTypes) = 0;

  virtual double& operator()(unsigned int branchId, size_t siteIndex, size_t type) = 0;
  virtual const double& operator()(unsigned int branchId, size_t siteIndex, size_t type) const = 0;
};


/**
 * @brief Partial implementation of the substitution mapping interface.
 *
 * This implementation copies the input tree in a TreeTemplate<Node> object.
 */
class AbstractSubstitutionMapping :
  virtual public SubstitutionMapping,
  virtual public AbstractMapping
{
private:
  size_t numberOfTypes_;

public:
  AbstractSubstitutionMapping() : AbstractMapping(), numberOfTypes_(0){}

  AbstractSubstitutionMapping(const AbstractSubstitutionMapping& absm) :
    AbstractMapping(absm), numberOfTypes_(absm.numberOfTypes_) {}

  AbstractSubstitutionMapping* clone() const = 0;

  AbstractSubstitutionMapping& operator=(const AbstractSubstitutionMapping& absm)
  {
    AbstractMapping::operator=(absm);
    numberOfTypes_ = absm.numberOfTypes_;

    return *this;
  }

  virtual ~AbstractSubstitutionMapping() {}

  size_t getNumberOfSubstitutionTypes() const
  {
    return numberOfTypes_;
  }

  virtual void setNumberOfSubstitutionTypes(size_t numberOfTypes)
  {
    numberOfTypes_ = numberOfTypes;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_SUBSTITUTIONMAPPING_H
