// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_ABSTRACTCODONBGCSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_ABSTRACTCODONBGCSUBSTITUTIONMODEL_H

#include <Bpp/Numeric/AbstractParameterAliasable.h>

#include "CodonSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Abstract class for modelling of non-synonymous and
 * synonymous substitution rates in codon models, with gBGC.
 *
 * @author Laurent GuÃÂ©guen
 *
 * The non-synonymous substitution rate is multiplied with
 * @f$\frac{\epsilon B+S}{1-e^{-(\epsilon B+S)}}@f$.
 *
 * The synonymous substitution rate is multiplied with @f$\frac{\epsilon
 * B}{1-e^{-\epsilon B}}@f$.
 *
 *
 * with positive parameter @f$S@f$ that stands for selection, and real
 * parameter @f$B@f$ for biased gene conversion. In the formula,
 * @f$\epsilon = 1@f$ for AT->GC substitutions, @f$\epsilon = -1@f$
 * for GC->AT substitution, and  @f$\epsilon = 0@f$ otherwise.
 *
 *
 * References:
 * - Galtier N, Duret L, GlÃÂ©min S, Ranwez V (2009) GC-biased gene
 * conversion promotes the fixation of deleterious amino acid changes
 * in primates, Trends in Genetics, vol. 25(1) pp.1-5.
 */
class AbstractCodonBGCSubstitutionModel :
  public virtual CoreCodonSubstitutionModelInterface,
  public virtual AbstractParameterAliasable
{
private:
  std::shared_ptr<const GeneticCode> pgencode_;

  double B_, S_;

  std::shared_ptr<const StateMapInterface> stateMap_;

public:
  /**
   * @brief Build a new AbstractCodonBGCSubstitutionModel object.
   *
   * @param pgencode the genetic code
   * @param prefix the Namespace
   */
  AbstractCodonBGCSubstitutionModel(
      std::shared_ptr<const GeneticCode> pgencode,
      const std::string& prefix);

  AbstractCodonBGCSubstitutionModel(const AbstractCodonBGCSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    pgencode_(model.pgencode_),
    B_(model.B_),
    S_(model.S_),
    stateMap_(model.stateMap_)
  {}

  AbstractCodonBGCSubstitutionModel& operator=(
      const AbstractCodonBGCSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    pgencode_ = model.pgencode_;
    B_ = model.B_;
    S_ = model.S_;
    stateMap_ = model.stateMap_;

    return *this;
  }

  AbstractCodonBGCSubstitutionModel* clone() const override
  {
    return new AbstractCodonBGCSubstitutionModel(*this);
  }

  virtual ~AbstractCodonBGCSubstitutionModel() {}

public:
  void fireParameterChanged(const ParameterList& parameters) override;

  double getCodonsMulRate(size_t i, size_t j) const override;

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    throw NullPointerException("AbstractCodonBGCSubstitutionModel::frequencySet. No associated FrequencySet.");
  }

  bool hasCodonFrequencySet() const override
  {
    return false;
  }

  void setFreq(std::map<int, double>& frequencies) override {}
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTCODONBGCSUBSTITUTIONMODEL_H
