// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_PROTEIN_WAG01_H
#define BPP_PHYL_MODEL_PROTEIN_WAG01_H


#include "../AbstractSubstitutionModel.h"
#include "../FrequencySet/ProteinFrequencySet.h"
#include "ProteinSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>

namespace bpp
{
/**
 * @brief The Whelan and Goldman substitution model for proteins.
 *
 * The exchangability matrix is normalized so that \f$Q = S . \pi\f$ and
 * \f$\sum_i Q_{i,i}\pi_i = -1\f$.
 *
 * The original frequencies can be used, or alternatively a
 * parametrized version, corresponding to the so-called WAG01+F
 * model. Eigen values and vectors are obtained numerically.
 *
 * Reference:
 *
 * Whelan, S. and N. Goldman. 2001. A general empirical model of
 * protein evolution derived from multiple protein families using a
 * maximum likelihood approach. Molecular Biology and Evolution 18:691-699.
 */
class WAG01 :
  public AbstractReversibleProteinSubstitutionModel
{
private:
  std::unique_ptr<ProteinFrequencySetInterface> freqSet_;

public:
  /**
   * @brief Build a simple WAG01 model, with original equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   */
  WAG01(std::shared_ptr<const ProteicAlphabet> alpha);

  /**
   * @brief Build a WAG01 model with special equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   * @param freqSet A pointer toward a protein frequencies set, which will be owned by this instance.
   * @param initFreqs Tell if the frequency set should be initialized with the original WAG01 values.
   * Otherwise, the values of the set will be used.
   */
  WAG01(
      std::shared_ptr<const ProteicAlphabet> alpha,
      std::unique_ptr<ProteinFrequencySetInterface> freqSet,
      bool initFreqs = false);

  WAG01(const WAG01& model) :
    AbstractParameterAliasable(model),
    AbstractReversibleProteinSubstitutionModel(model),
    freqSet_(model.freqSet_->clone())
  {}

  WAG01& operator=(const WAG01& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractReversibleProteinSubstitutionModel::operator=(model);
    freqSet_.reset(model.freqSet_->clone());
    return *this;
  }

  virtual ~WAG01() {}

  WAG01* clone() const override { return new WAG01(*this); }

public:
  std::string getName() const override
  {
    if (freqSet_->getNamespace().find("WAG01+F.") != std::string::npos)
      return "WAG01+F";
    else
      return "WAG01";
  }

  void fireParameterChanged(const ParameterList& parameters) override
  {
    freqSet_->matchParametersValues(parameters);
    freq_ = freqSet_->getFrequencies();
    AbstractReversibleSubstitutionModel::fireParameterChanged(parameters);
  }

  void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix);
    freqSet_->setNamespace(prefix + freqSet_->getName() + ".");
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
    throw NullPointerException("WAG01::frequencySet(). No associated FrequencySet.");
  }
    
  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override;

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_WAG01_H
