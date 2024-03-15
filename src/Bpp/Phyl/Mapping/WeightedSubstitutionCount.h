// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_WEIGHTEDSUBSTITUTIONCOUNT_H
#define BPP_PHYL_MAPPING_WEIGHTEDSUBSTITUTIONCOUNT_H


#include "SubstitutionCount.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Interface allowing for weighting of substitution counts
 * according to state properties.
 *
 * These weights are used for the final counts, after the mapping
 * process (see SubstitutionDistances for integration of substitution
 * distances before the mapping process).
 */

class WeightedSubstitutionCount :
  public virtual SubstitutionCountInterface
{
public:
  virtual void setWeights(std::shared_ptr<const AlphabetIndex2> index) = 0;
  virtual bool hasWeights() const = 0;
  virtual std::shared_ptr<const AlphabetIndex2> getWeights() const = 0;
};

/**
 * @brief Partial implementation of the WeightedSubstitutionCount interface.
 */
class AbstractWeightedSubstitutionCount :
  public virtual WeightedSubstitutionCount
{
protected:
  std::shared_ptr<const AlphabetIndex2> weights_;

public:
  AbstractWeightedSubstitutionCount(std::shared_ptr<const AlphabetIndex2> weights) :
    weights_(weights)
  {}

  AbstractWeightedSubstitutionCount(const AbstractWeightedSubstitutionCount& index) :
    weights_(index.weights_)
  {}

  AbstractWeightedSubstitutionCount& operator=(const AbstractWeightedSubstitutionCount& index)
  {
    weights_ = index.weights_;

    return *this;
  }

  virtual ~AbstractWeightedSubstitutionCount()
  {}

public:
  void setWeights(std::shared_ptr<const AlphabetIndex2> weights);
  bool hasWeights() const { return weights_.get() != 0; }
  std::shared_ptr<const AlphabetIndex2> getWeights() const { return weights_; }

protected:
  virtual void weightsHaveChanged() = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_WEIGHTEDSUBSTITUTIONCOUNT_H
