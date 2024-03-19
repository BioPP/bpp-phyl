// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_PROTEIN_LG08_H
#define BPP_PHYL_MODEL_PROTEIN_LG08_H


#include "../AbstractSubstitutionModel.h"
#include "../FrequencySet/ProteinFrequencySet.h"
#include "ProteinSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>

namespace bpp
{
/**
 * @brief The Le and Gascuel substitution model for proteins.
 *
 * The exchangability matrix is normalized so that \f$Q = S . \pi\f$ and \f$\sum_i Q_{i,i}\pi_i = -1\f$.
 * The original frequencies can be used, or alternatively a parametrized version, corresponding to the
 * so-called LG08+F model.
 * Eigen values and vectors are obtained numerically.
 *
 * Reference:
 * - Le, Si Q. and Gascuel, O. (2008), _Molecular Biology And Evolution_ 25, 1307--1320.
 */
class LG08 :
  public AbstractReversibleProteinSubstitutionModel
{
private:
  std::unique_ptr<ProteinFrequencySetInterface> freqSet_;

public:
  /**
   * @brief Build a simple LG08 model, with original equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   */
  LG08(std::shared_ptr<const ProteicAlphabet> alpha);

  /**
   * @brief Build a LG08 model with special equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   * @param freqSet A pointer toward a protein frequencies set, which will be owned by this instance.
   * @param initFreqs Tell if the frequency set should be initialized with the original LG08 values.
   * Otherwise, the values of the set will be used.
   */
  LG08(
      std::shared_ptr<const ProteicAlphabet> alpha,
      std::unique_ptr<ProteinFrequencySetInterface> freqSet,
      bool initFreqs = false);

  LG08(const LG08& model) :
    AbstractParameterAliasable(model),
    AbstractReversibleProteinSubstitutionModel(model),
    freqSet_(model.freqSet_->clone())
  {}

  LG08& operator=(const LG08& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractReversibleProteinSubstitutionModel::operator=(model);
    freqSet_.reset(model.freqSet_->clone());
    return *this;
  }

  virtual ~LG08() {}

  LG08* clone() const override { return new LG08(*this); }

public:
  std::string getName() const override
  {
    if (freqSet_->getNamespace().find("LG08+F.") != std::string::npos)
      return "LG08+F";
    else
      return "LG08";
  }

  void fireParameterChanged(const ParameterList& parameters) override
  {
    freqSet_->matchParametersValues(parameters);
    freq_ = freqSet_->getFrequencies();
    AbstractReversibleSubstitutionModel::fireParameterChanged(parameters);
  }

  void setFrequencySet(const ProteinFrequencySetInterface& freqSet)
  {
    freqSet_.reset(dynamic_cast<ProteinFrequencySetInterface*>(freqSet.clone()));
    resetParameters_();
    addParameters_(freqSet_->getParameters());
  }

  void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix);
    freqSet_->setNamespace(prefix + freqSet_->getName() + ".");
  }

  const FrequencySetInterface& frequencySet() const override
  {
    if (freqSet_)
      return *freqSet_;
    throw NullPointerException("LG08::frequencySet(). No associated FrequencySet.");
  }

  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_LG08_H
