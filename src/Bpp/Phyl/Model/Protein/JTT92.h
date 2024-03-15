// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_PROTEIN_JTT92_H
#define BPP_PHYL_MODEL_PROTEIN_JTT92_H


#include "../AbstractSubstitutionModel.h"
#include "../FrequencySet/ProteinFrequencySet.h"
#include "ProteinSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>

namespace bpp
{
/**
 * @brief The Jones, Taylor and Thornton substitution model for proteins.
 *
 * Exchangeabilities have been computed using the DCMut method of Kosiol and Goldman.
 * The exchangability matrix is normalized so that \f$Q = S . \pi\f$ and \f$\sum_i Q_{i,i}\pi_i = -1\f$.
 * The original frequencies can be used, or alternatively a parametrized version, corresponding to the
 * so-called JTT92+F model.
 * Eigen values and vectors are obtained numerically.
 *
 * References:
 * - Jones DT, Taylor WR and Thornton JM (1992), _Computer Applications In The Biosciences_, 8(3) 275-82.
 * - Kosiol C and Goldman N (2005), _Molecular Biology And Evolution_ 22(2) 193-9.
 */
class JTT92 :
  public AbstractReversibleProteinSubstitutionModel
{
private:
  std::shared_ptr<ProteinFrequencySetInterface> freqSet_;

public:
  /**
   * @brief Build a simple JTT92 model, with original equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   */
  JTT92(std::shared_ptr<const ProteicAlphabet> alpha);

  /**
   * @brief Build a JTT92 model with special equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   * @param freqSet A pointer toward a protein frequencies set, which will be owned by this instance.
   * @param initFreqs Tell if the frequency set should be initialized with the original JTT92 values.
   * Otherwise, the values of the set will be used.
   */
  JTT92(
      std::shared_ptr<const ProteicAlphabet> alpha,
      std::unique_ptr<ProteinFrequencySetInterface> freqSet,
      bool initFreqs = false);

  JTT92(const JTT92& model) :
    AbstractParameterAliasable(model),
    AbstractReversibleProteinSubstitutionModel(model),
    freqSet_(model.freqSet_->clone())
  {}

  JTT92& operator=(const JTT92& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractReversibleProteinSubstitutionModel::operator=(model);
    freqSet_.reset(model.freqSet_->clone());
    return *this;
  }

  virtual ~JTT92() {}

  JTT92* clone() const override { return new JTT92(*this); }

public:
  std::string getName() const override
  {
    if (freqSet_->getNamespace().find("JTT92+F.") != std::string::npos)
      return "JTT92+F";
    else
      return "JTT92";
  }

  void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix);
    freqSet_->setNamespace(prefix + freqSet_->getName() + ".");
  }


  void fireParameterChanged(const ParameterList& parameters) override
  {
    freqSet_->matchParametersValues(parameters);
    freq_ = freqSet_->getFrequencies();
    AbstractReversibleSubstitutionModel::fireParameterChanged(parameters);
  }

  void setFrequencySet(const ProteinFrequencySetInterface& freqSet)
  {
    freqSet_.reset(freqSet.clone());
    resetParameters_();
    addParameters_(freqSet_->getParameters());
  }

  const FrequencySetInterface& frequencySet() const override
  {
    if (freqSet_)
      return *freqSet_;
    throw NullPointerException("JTT92::frequencySet(). No associated FrequencySet.");
  }
    
  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_JTT92_H
