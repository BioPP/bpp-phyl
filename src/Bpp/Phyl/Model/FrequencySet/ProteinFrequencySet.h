// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_FREQUENCYSET_PROTEINFREQUENCYSET_H
#define BPP_PHYL_MODEL_FREQUENCYSET_PROTEINFREQUENCYSET_H

#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>

#include "FrequencySet.h"

namespace bpp
{
/**
 * @brief Parametrize a set of state frequencies for proteins.
 */
class ProteinFrequencySetInterface :
  public virtual FrequencySetInterface
{
public:
  ProteinFrequencySetInterface* clone() const override = 0;

  virtual std::shared_ptr<const ProteicAlphabet> getProteicAlphabet() const = 0;
};

/**
 * @brief Protein FrequencySet using 19 independent parameters to
 * model the 20 frequencies.
 *
 * The parameters are called @f$ \theta_{i \in 1..19} @f$, and are
 * initialized so that all frequencies are equal to 0.005. The
 * parametrization depends on the method used. Default
 * method is 1 (ie global ratio).
 *
 * @see Simplex
 */
class FullProteinFrequencySet :
  public virtual ProteinFrequencySetInterface,
  public FullFrequencySet
{
public:
  FullProteinFrequencySet(
      std::shared_ptr<const ProteicAlphabet> alphabet,
      bool allowNullFreqs = false,
      unsigned short method = 1,
      const std::string& name = "Full") :
    FullFrequencySet(
      std::make_shared<const CanonicalStateMap>(alphabet, false),
      allowNullFreqs,
      method,
      name)
  {}

  FullProteinFrequencySet(
      std::shared_ptr<const ProteicAlphabet> alphabet,
      const std::vector<double>& initFreqs,
      bool allowNullFreqs = false,
      unsigned short method = 1,
      const std::string& name = "Full") :
    FullFrequencySet(
      std::make_shared<const CanonicalStateMap>(alphabet, false),
      initFreqs,
      allowNullFreqs,
      method,
      name)
  {}

  FullProteinFrequencySet* clone() const override { return new FullProteinFrequencySet(*this); }

public:
  std::shared_ptr<const ProteicAlphabet> getProteicAlphabet() const override
  {
    return std::dynamic_pointer_cast<const ProteicAlphabet>(getAlphabet());
  }
};

/**
 * @brief FrequencySet useful for homogeneous and stationary models, protein implementation
 *
 * This set contains no parameter.
 */
class FixedProteinFrequencySet :
  public virtual ProteinFrequencySetInterface,
  public FixedFrequencySet
{
public:
  FixedProteinFrequencySet(
      std::shared_ptr<const ProteicAlphabet> alphabet,
      const std::vector<double>& initFreqs,
      const std::string& name = "Fixed") :
    FixedFrequencySet(
      std::make_shared<const CanonicalStateMap>(alphabet, false),
      initFreqs,
      name)
  {}

  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet.
   */
  FixedProteinFrequencySet(
      std::shared_ptr<const ProteicAlphabet> alphabet,
      const std::string& name = "Fixed") :
    FixedFrequencySet(
      std::make_shared<const CanonicalStateMap>(alphabet, false),
      name)
  {}

  FixedProteinFrequencySet* clone() const override { return new FixedProteinFrequencySet(*this); }

  std::shared_ptr<const ProteicAlphabet> getProteicAlphabet() const override
  {
    return std::dynamic_pointer_cast<const ProteicAlphabet>(getAlphabet());
  }
};

/**
 * @brief FrequencySet from file
 *
 * This set contains no parameter.
 */

class UserProteinFrequencySet :
  public virtual ProteinFrequencySetInterface,
  public UserFrequencySet
{
public:
  UserProteinFrequencySet(
      std::shared_ptr<const ProteicAlphabet> alphabet,
      const std::string& path,
      size_t nCol = 1) :
    UserFrequencySet(
      std::make_shared<const CanonicalStateMap>(alphabet, false),
      path,
      nCol)
  {}

  UserProteinFrequencySet* clone() const override { return new UserProteinFrequencySet(*this); }

  std::shared_ptr<const ProteicAlphabet> getProteicAlphabet() const override
  {
    return std::dynamic_pointer_cast<const ProteicAlphabet>(getAlphabet());
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_FREQUENCYSET_PROTEINFREQUENCYSET_H
