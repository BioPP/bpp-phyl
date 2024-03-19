// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_ABSTRACTCODONCPGSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_ABSTRACTCODONCPGSUBSTITUTIONMODEL_H

#include <Bpp/Numeric/AbstractParameterAliasable.h>

#include "CodonSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Abstract class for modelling of CpG -> CpA or TpG (symetric)
 *  hypermutability substitution rate inside codons. Note that the
 *  neihbouring effects between codons are not considered.
 *
 * @author Laurent GuÃÂ©guen
 *
 * Substitution rate from C to T (resp. from G to A) is multiplied by
 *  a factor @f$\rho@f$ if C is followed by a G (resp. if G is
 *  following a C).
 *
 * Hypermutability parameter is named \c "rho".
 */
class AbstractCodonCpGSubstitutionModel :
  public virtual CoreCodonSubstitutionModelInterface,
  public virtual AbstractParameterAliasable
{
private:
  double rho_;

  std::shared_ptr<const StateMapInterface> stateMap_;

public:
  /**
   * @brief Build a new AbstractCodonCpGSubstitutionModel object from
   *  a pointer to NucleotideSubstitutionModel.
   *
   * @param alphabet the codon alphabet
   * @param prefix the Namespace
   */
  AbstractCodonCpGSubstitutionModel(
      std::shared_ptr<const CodonAlphabet> alphabet,
      const std::string& prefix);

  AbstractCodonCpGSubstitutionModel(const AbstractCodonCpGSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    rho_(model.rho_),
    stateMap_(model.stateMap_)
  {}

  AbstractCodonCpGSubstitutionModel& operator=(
      const AbstractCodonCpGSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    rho_ = model.rho_;
    stateMap_ = model.stateMap_;
    return *this;
  }

  AbstractCodonCpGSubstitutionModel* clone() const override
  {
    return new AbstractCodonCpGSubstitutionModel(*this);
  }

  virtual ~AbstractCodonCpGSubstitutionModel() {}

public:
  void fireParameterChanged(const ParameterList& parameters) override;

  double getCodonsMulRate(size_t i, size_t j) const override;

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    throw NullPointerException("AbstractCodonCpGSubstitutionModel::frequencySet. No associated FrequencySet.");
  }

  bool hasCodonFrequencySet() const override
  {
    return false;
  }

  void setFreq(std::map<int, double>& frequencies) override {}
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTCODONCPGSUBSTITUTIONMODEL_H
