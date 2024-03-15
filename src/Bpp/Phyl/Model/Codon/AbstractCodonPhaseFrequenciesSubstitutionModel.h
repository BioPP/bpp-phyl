// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_ABSTRACTCODONPHASEFREQUENCIESSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_ABSTRACTCODONPHASEFREQUENCIESSUBSTITUTIONMODEL_H


#include "../FrequencySet/CodonFrequencySet.h"
#include "CodonSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Abstract Class for substitution models on codons
 *  parametrized by a frequency.
 *
 * @author Laurent GuÃÂ©guen
 *
 * This class should be used with models which equilibrium
 * distribution is fixed, ans does not depend on the parameters.
 * Otherwise there may be problems of identifiability of the
 * parameters.
 *
 * If we denote @f$F@f$ the given frequencies for codons,
 * @f$F_{j_k}@f$ is the frequency of letter @f$j@f$ in phase
 * @f$k@f$.
 *
 * For codons @f$i=i_1i_2i_3@f$ and @f$j=j_1j_2j_3@f$, the generator
 * term defined from inherited and inheriting classes,
 * @f$Q_{ij})@f$, is multiplied by the product of the @f$F_{j_k}@f$
 * for each @f$k \in 1, 2, 3@f$ such that @f$i_k \neq j_k@f$.
 */
class AbstractCodonPhaseFrequenciesSubstitutionModel :
  public virtual CoreCodonSubstitutionModelInterface,
  public virtual AbstractParameterAliasable
{
private:
  /**
   * @brief Position dependent version of Codon Frequencies Set
   */
  std::unique_ptr<CodonFrequencySetInterface> posFreqSet_;
  std::string freqName_;

public:
  /**
   * @brief Build a AbstractCodonPhaseFrequenciesSubstitutionModel instance
   *
   * @param pfreq pointer to the AbstractFrequencySet equilibrium frequencies.
   *        It is owned by the instance.
   * @param prefix the Namespace
   */
  AbstractCodonPhaseFrequenciesSubstitutionModel(
    std::unique_ptr<CodonFrequencySetInterface> pfreq,
    const std::string& prefix);

  AbstractCodonPhaseFrequenciesSubstitutionModel(const AbstractCodonPhaseFrequenciesSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    posFreqSet_(model.posFreqSet_->clone()),
    freqName_(model.freqName_)
  {}

  AbstractCodonPhaseFrequenciesSubstitutionModel& operator=(const AbstractCodonPhaseFrequenciesSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    posFreqSet_.reset(model.posFreqSet_->clone());
    freqName_ = model.freqName_;

    return *this;
  }

  AbstractCodonPhaseFrequenciesSubstitutionModel* clone() const override
  {
    return new AbstractCodonPhaseFrequenciesSubstitutionModel(*this);
  }

  virtual ~AbstractCodonPhaseFrequenciesSubstitutionModel();

  void fireParameterChanged(const ParameterList& parameters) override;

  void setFreq(std::map<int, double>& frequencies) override;

  void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix);
    posFreqSet_->setNamespace(prefix + freqName_);
  }

  double getCodonsMulRate(size_t, size_t) const override;

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    return *posFreqSet_;
  }

  const FrequencySetInterface& frequencySet() const 
  {
    return *posFreqSet_;
  }

  bool hasCodonFrequencySet() const override
  {
    return (posFreqSet_ != nullptr);
  }

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTCODONPHASEFREQUENCIESSUBSTITUTIONMODEL_H
