// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_PROTEIN_USERPROTEINSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_PROTEIN_USERPROTEINSUBSTITUTIONMODEL_H


#include "../AbstractSubstitutionModel.h"
#include "../FrequencySet/ProteinFrequencySet.h"
#include "ProteinSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>

// From the STL:
#include <string>

namespace bpp
{
/**
 * @brief Build an empirical protein substitution model from a file.
 *
 * The file must follow PAML's format, and contain the exchangeabilities components (\f$S_{i,j}\f$)
 * and all equilibrium frequencies (\f$\pi_{i}\f$).
 * The generator is build so that \f$Q_{i,j} = \pi_i . S_{i,j}\f$, and is normalized
 * so that \f$\sum_i Q_{i,i} \times \pi_i = -1\f$.
 */
class UserProteinSubstitutionModel :
  public AbstractReversibleProteinSubstitutionModel
{
private:
  std::string path_;
  std::unique_ptr<ProteinFrequencySetInterface> freqSet_;

public:
  /**
   * @brief Build a protein model from a PAML file, with original equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   * @param path The path toward the file to parse.
   * @param prefix The parameter namespace to use.
   */
  UserProteinSubstitutionModel(
    std::shared_ptr<const ProteicAlphabet> alpha,
    const std::string& path,
    const std::string& prefix);

  /**
   * @brief Build a protein model from a PAML file, with special equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   * @param freqSet A pointer toward a protein frequencies set, which will be owned by this instance.
   * @param path The path toward the file to parse.
   * @param prefix The parameter namespace to use.
   * @param initFreqs Tell if the frequency set should be initialized with the original JTT92 values.
   * Otherwise, the values of the set will be used.
   */
  UserProteinSubstitutionModel(
    std::shared_ptr<const ProteicAlphabet> alpha,
    const std::string& path,
    std::unique_ptr<ProteinFrequencySetInterface> freqSet,
    const std::string& prefix,
    bool initFreqs = false
    );

  UserProteinSubstitutionModel(const UserProteinSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    AbstractReversibleProteinSubstitutionModel(model),
    path_(model.path_),
    freqSet_(model.freqSet_->clone())
  {}

  UserProteinSubstitutionModel& operator=(const UserProteinSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractReversibleProteinSubstitutionModel::operator=(model);
    path_ = model.path_;
    freqSet_.reset(model.freqSet_->clone());
    return *this;
  }

  virtual ~UserProteinSubstitutionModel() {}

  UserProteinSubstitutionModel* clone() const override { return new UserProteinSubstitutionModel(*this); }

public:
  std::string getName() const override;
  
  const std::string& getPath() const { return path_; }

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
    throw NullPointerException("UserProteinSubstitutionModel::frequencySet(). No associated FrequencySet.");
  }
    
  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override;

protected:
  void readFromFile();
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_USERPROTEINSUBSTITUTIONMODEL_H
