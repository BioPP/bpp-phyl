// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_ABSTRACTCODONDISTANCESUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_ABSTRACTCODONDISTANCESUBSTITUTIONMODEL_H

#include <Bpp/Numeric/AbstractParameterAliasable.h>

#include "CodonSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Abstract class for modelling of non-synonymous and
 *  synonymous substitution rates in codon models.
 *
 * @author Laurent Guéguen
 *
 * If a distance @f$d@f$ between amino-acids is defined, the
 *  non-synonymous rate is multiplied with, if the coded amino-acids
 *  are @f$x@f$ and @f$y@f$, @f$\beta*\exp(-\alpha.d(x,y))@f$ with
 *  non-negative parameter \c "alpha" and positive parameter \c
 *  "beta".
 *
 * If such a distance is not defined, the non-synonymous substitution
 *  rate is multiplied with @f$\beta@f$ with positive parameter \c
 *  "beta" (ie @f$d=0@f$).
 *
 * If paramSynRate is true, the synonymous substitution rate is
 *  multiplied with @f$\gamma@f$ (with optional positive parameter \c
 *  "gamma"), else it is multiplied with 1.
 *
 * References:
 * - Goldman N. and Yang Z. (1994), _Molecular Biology And Evolution_ 11(5) 725--736.
 * - Kosakovsky Pond, S. and Muse, S.V. (2005), _Molecular Biology And Evolution_,
 *   22(12), 2375--2385.
 * - Mayrose, I. and Doron-Faigenboim, A. and Bacharach, E. and Pupko T.
 *   (2007), Bioinformatics, 23, i319--i327.
 */
class AbstractCodonDistanceSubstitutionModel :
  public virtual CoreCodonSubstitutionModelInterface,
  public virtual AbstractParameterAliasable
{
private:
  std::shared_ptr<const AlphabetIndex2> pdistance_;

  std::shared_ptr<const GeneticCode> pgencode_;

  double alpha_, beta_;

  double gamma_;

  std::shared_ptr<const StateMapInterface> stateMap_;

public:
  /**
   * @brief Build a new AbstractCodonDistanceSubstitutionModel object.
   *
   * @param pdist optional pointer to a distance between amino-acids
   * @param pgencode the genetic code
   * @param prefix the Namespace
   * @param paramSynRate is true iff synonymous rate is parametrised
   *       (default=false).
   */
  AbstractCodonDistanceSubstitutionModel(
      std::shared_ptr<const AlphabetIndex2> pdist,
      std::shared_ptr<const GeneticCode> pgencode,
      const std::string& prefix,
      bool paramSynRate = false);

  AbstractCodonDistanceSubstitutionModel(const AbstractCodonDistanceSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    pdistance_(model.pdistance_),
    pgencode_(model.pgencode_),
    alpha_(model.alpha_),
    beta_(model.beta_),
    gamma_(model.gamma_),
    stateMap_(model.stateMap_)
  {}

  AbstractCodonDistanceSubstitutionModel& operator=(
      const AbstractCodonDistanceSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    pdistance_ = model.pdistance_;
    pgencode_ = model.pgencode_;
    alpha_ = model.alpha_;
    beta_ = model.beta_;
    gamma_ = model.gamma_;
    stateMap_ = model.stateMap_;

    return *this;
  }

  AbstractCodonDistanceSubstitutionModel* clone() const override
  {
    return new AbstractCodonDistanceSubstitutionModel(*this);
  }

  virtual ~AbstractCodonDistanceSubstitutionModel() {}

public:
  void fireParameterChanged(const ParameterList& parameters) override;

  double getCodonsMulRate(size_t i, size_t j) const override;

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    throw NullPointerException("AbstractCodonDistanceSubstitutionModel::frequencySet. No associated FrequencySet.");
  }

  bool hasCodonFrequencySet() const override
  {
    return false;
  }

  void setFreq(std::map<int, double>& frequencies) override {}
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTCODONDISTANCESUBSTITUTIONMODEL_H
