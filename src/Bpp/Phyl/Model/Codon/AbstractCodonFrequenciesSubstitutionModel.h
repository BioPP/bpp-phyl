// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_ABSTRACTCODONFREQUENCIESSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_ABSTRACTCODONFREQUENCIESSUBSTITUTIONMODEL_H


#include "../FrequencySet/CodonFrequencySet.h"
#include "CodonSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Abstract Class for substitution models on codons
 *  parametrized by frequencies.
 *
 * This class should be used with models which equilibrium
 * distribution is fixed, and does not depend on the parameters.
 * Otherwise there may be problems of identifiability of the
 * parameters.
 *
 * @author Laurent Guéguen
 *
 * If we denote @f$F@f$ the equilibrium frequency, the generator term
 * defined from inherited and inheriting classes, @f$Q_{ij})@f$, is
 * multiplied by @f$F_{j}@f$.
 */
class AbstractCodonFrequenciesSubstitutionModel :
  public virtual CoreCodonSubstitutionModelInterface,
  public virtual AbstractParameterAliasable
{
private:
  std::unique_ptr<CodonFrequencySetInterface> pfreqset_;
  std::string freqName_;

public:
  /**
   * @brief Build a AbstractCodonFrequenciesSubstitutionModel instance
   *
   * @param pfreq pointer to the AbstractFrequencySet equilibrium frequencies.
   *        It is owned by the instance.
   * @param prefix the Namespace
   */
  AbstractCodonFrequenciesSubstitutionModel(
      std::unique_ptr<CodonFrequencySetInterface> pfreq,
      const std::string& prefix);

  AbstractCodonFrequenciesSubstitutionModel(const AbstractCodonFrequenciesSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    pfreqset_(model.pfreqset_->clone()),
    freqName_(model.freqName_)
  {}

  AbstractCodonFrequenciesSubstitutionModel& operator=(const AbstractCodonFrequenciesSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    pfreqset_   = std::unique_ptr<CodonFrequencySetInterface>(model.pfreqset_->clone());
    freqName_   = model.freqName_;
    return *this;
  }

  AbstractCodonFrequenciesSubstitutionModel* clone() const override
  {
    return new AbstractCodonFrequenciesSubstitutionModel(*this);
  }

  virtual ~AbstractCodonFrequenciesSubstitutionModel();

  void fireParameterChanged(const ParameterList& parameters) override;

  void setFreq(std::map<int, double>& frequencies) override;

  void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix);
    pfreqset_->setNamespace(prefix + freqName_);
  }

  double getCodonsMulRate(size_t, size_t) const override;

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    return *pfreqset_;
  }

  const FrequencySetInterface& frequencySet() const
  {
    return *pfreqset_;
  }

  bool hasCodonFrequencySet() const override
  {
    return pfreqset_ != nullptr;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTCODONFREQUENCIESSUBSTITUTIONMODEL_H
