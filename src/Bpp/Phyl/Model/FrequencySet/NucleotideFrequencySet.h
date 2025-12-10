// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_FREQUENCYSET_NUCLEOTIDEFREQUENCYSET_H
#define BPP_PHYL_MODEL_FREQUENCYSET_NUCLEOTIDEFREQUENCYSET_H

#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

#include "FrequencySet.h"

namespace bpp
{
/**
 * @brief Parametrize a set of state frequencies for nucleotides.
 */
class NucleotideFrequencySetInterface :
  public virtual FrequencySetInterface
{
public:
  NucleotideFrequencySetInterface* clone() const override = 0;

  virtual std::shared_ptr<const NucleicAlphabet> getNucleicAlphabet() const = 0;
};

/**
 * @brief Nucleotide FrequencySet using only one parameter, the GC
 *        content (denoted as 'GC.theta')
 */

class GCFrequencySet :
  public virtual NucleotideFrequencySetInterface,
  public AbstractFrequencySet
{
public:
  GCFrequencySet(std::shared_ptr<const NucleicAlphabet> alphabet) :
    AbstractFrequencySet(std::make_shared<const CanonicalStateMap>(alphabet, false), "GC.", "GC")
  {
    addParameter_(new Parameter("GC.theta", 0.5, Parameter::PROP_CONSTRAINT_IN));
    getFreq_(0) = getFreq_(1) = getFreq_(2) = getFreq_(3) = 0.25;
  }

  GCFrequencySet(std::shared_ptr<const NucleicAlphabet> alphabet, double theta) :
    AbstractFrequencySet(std::make_shared<const CanonicalStateMap>(alphabet, false), "GC.", "GC")
  {
    addParameter_(new Parameter("GC.theta", theta, Parameter::PROP_CONSTRAINT_IN));
    getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
    getFreq_(1) = getFreq_(2) = theta / 2.;
  }

  GCFrequencySet* clone() const override
  {
    return new GCFrequencySet(*this);
  }

  GCFrequencySet(const GCFrequencySet& gcf) :
    AbstractFrequencySet(gcf)
  {}

public:
  std::shared_ptr<const NucleicAlphabet> getNucleicAlphabet() const override
  {
    return std::dynamic_pointer_cast<const NucleicAlphabet>(getAlphabet());
  }

  void setFrequencies(const std::vector<double>& frequencies) override;

protected:
  void fireParameterChanged(const ParameterList& parameters) override;
};

/**
 * @brief Nucleotide FrequencySet using three independent parameters
 * (theta, theta1, theta2) to modelize the four frequencies:
 *
 * \f[
 * \begin{cases}
 * \theta = \pi_C + \pi_G\\
 * \theta_1 = \frac{\pi_A}{1 - \theta} = \frac{\pi_A}{\pi_A + \pi_T}\\
 * \theta_2 = \frac{\pi_G}{\theta} = \frac{\pi_G}{\pi_C + \pi_G}\\
 * \end{cases}
 * \Longleftrightarrow
 * \begin{cases}
 * \pi_A = \theta_1 (1 - \theta)\\
 * \pi_C = (1 - \theta_2) \theta\\
 * \pi_G = \theta_2 \theta\\
 * \pi_T = (1 - \theta_1)(1 - \theta).
 * \end{cases}
 * \f]
 *
 * with \f$\pi_x\f$ the frequency of nucleotide \f$x\f$.
 */
class FullNucleotideFrequencySet :
  public virtual NucleotideFrequencySetInterface,
  public AbstractFrequencySet
{
public:
  FullNucleotideFrequencySet(
      std::shared_ptr<const NucleicAlphabet> alphabet,
      bool allowNullFreqs = false,
      const std::string& name = "Full");

  FullNucleotideFrequencySet(
      std::shared_ptr<const NucleicAlphabet> alphabet,
      double theta, double theta1, double theta2,
      bool allowNullFreqs = false,
      const std::string& name = "Full");

  FullNucleotideFrequencySet* clone() const override { return new FullNucleotideFrequencySet(*this); }

public:
  std::shared_ptr<const NucleicAlphabet> getNucleicAlphabet() const override
  {
    return std::dynamic_pointer_cast<const NucleicAlphabet>(getAlphabet());
  }

  void setFrequencies(const std::vector<double>& frequencies) override;

protected:
  void fireParameterChanged(const ParameterList& parameters) override;
};


/**
 * @brief FrequencySet useful for homogeneous and stationary models, nucleotide implementation
 *
 * This set contains no parameter.
 */
class FixedNucleotideFrequencySet :
  public virtual NucleotideFrequencySetInterface,
  public FixedFrequencySet
{
public:
  FixedNucleotideFrequencySet(
      std::shared_ptr<const NucleicAlphabet> alphabet,
      const std::vector<double>& initFreqs,
      const std::string& name = "Fixed") :
    FixedFrequencySet(std::make_shared<const CanonicalStateMap>(alphabet, false), initFreqs, name) {}

  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet.
   */
  FixedNucleotideFrequencySet(
      std::shared_ptr<const NucleicAlphabet> alphabet,
      const std::string& name = "Fixed") :
    FixedFrequencySet(std::make_shared<const CanonicalStateMap>(alphabet, false), name) {}

  FixedNucleotideFrequencySet* clone() const override { return new FixedNucleotideFrequencySet(*this); }

  std::shared_ptr<const NucleicAlphabet> getNucleicAlphabet() const override
  {
    return std::dynamic_pointer_cast<const NucleicAlphabet>(getAlphabet());
  }
};

/**
 * @brief FrequencySet useful for homogeneous and stationary models, nucleotide implementation
 *
 * This set contains no parameter.
 */
class UserNucleotideFrequencySet :
  public virtual NucleotideFrequencySetInterface,
  public UserFrequencySet
{
public:
  UserNucleotideFrequencySet(
      std::shared_ptr<const NucleicAlphabet> alphabet,
      const std::string& path,
      size_t nCol = 1) :
    UserFrequencySet(std::make_shared<const CanonicalStateMap>(alphabet, false), path, nCol) {}

  UserNucleotideFrequencySet* clone() const override { return new UserNucleotideFrequencySet(*this); }

  std::shared_ptr<const NucleicAlphabet> getNucleicAlphabet() const override
  {
    return std::dynamic_pointer_cast<const NucleicAlphabet>(getAlphabet());
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_FREQUENCYSET_NUCLEOTIDEFREQUENCYSET_H
