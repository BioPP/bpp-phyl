// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_SUBSTITUTIONDISTANCE_H
#define BPP_PHYL_MAPPING_SUBSTITUTIONDISTANCE_H


#include "SubstitutionCount.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Interface allowing for using distances between states in
 * substitution counts.
 *
 * These distances are used for integration of substitution distances
 * before the mapping process (see WeightedSubstitutionCount for
 * weight on the final counts, after the mapping process).
 *
 */

class SubstitutionDistance :
  public virtual SubstitutionCountInterface
{
public:
  virtual void setDistances(std::shared_ptr<const AlphabetIndex2> index) = 0;
  virtual bool hasDistances() const = 0;
  virtual std::shared_ptr<const AlphabetIndex2> getDistances() const = 0;
};

/**
 * @brief Partial implementation of the SubstitutionDistance interface.
 */
class AbstractSubstitutionDistance :
  public virtual SubstitutionDistance
{
protected:
  std::shared_ptr<const AlphabetIndex2> distances_;

public:
  AbstractSubstitutionDistance(std::shared_ptr<const AlphabetIndex2> distances) :
    distances_(distances)
  {}

  AbstractSubstitutionDistance(const AbstractSubstitutionDistance& index) :
    distances_(index.distances_)
  {}

  AbstractSubstitutionDistance& operator=(const AbstractSubstitutionDistance& index)
  {
    distances_ = index.distances_;

    return *this;
  }

  virtual ~AbstractSubstitutionDistance()
  {}

public:
  void setDistances(std::shared_ptr<const AlphabetIndex2> distances);
  bool hasDistances() const { return distances_.get() != 0; }
  std::shared_ptr<const AlphabetIndex2> getDistances() const { return distances_; }

protected:
  virtual void distancesHaveChanged() = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_SUBSTITUTIONDISTANCE_H
