// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_ABSTRACTDFPSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_ABSTRACTDFPSUBSTITUTIONMODEL_H


#include "../AbstractSubstitutionModel.h"
#include "CodonSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>

namespace bpp
{
/**
 * @brief Class for neutral substitution models on triplets, following
 * the mutation process proposed in Doron-Fagenboim & Pupko, 2006, but
 * without equilibrium frequencies. This model is an extension of
 * Kimura 2-rates substitution model to codons.
 *
 * There are 5 five free parameters:
 *
 * \c "tr" : rate between codons that differ by 1 transition
 * \c "tv" : rate between codons that differ by 1 transversion. This
 *           rate is set to 1, and other parameters are ratio on this one (so
 *           'tr' is similar to 'kappa' in YN98 model).
 * \c "trr" : rate between codons that differ by 2 transversions
 * \c "tvv" : rate between codons that differ by 2 transversions
 * \c "trv" : rate between codons that differ by 1 transition & 1 transversion
 * \c "tsub" : rate between codons that differ by 3 substitutions
 *
 * Reference: Adi Doron-Faigenboim, Tal Pupko, 2007, A Combined
 * Empirical and Mechanistic Codon Model, Molecular Biology and
 * Evolution, Volume 24, Issue 2, Pages 388Ã¢ÂÂ397,
 * https://doi.org/10.1093/molbev/msl175
 *
 */

class AbstractDFPSubstitutionModel :
  public virtual CodonSubstitutionModelInterface,
  public AbstractSubstitutionModel
{
private:
  std::shared_ptr<const GeneticCode> gCode_;

  double tr_, trr_, tvv_, trv_, tsub_;

public:
  /**
   * @brief Build a new AbstractDFPSubstitutionModel object
   */
  AbstractDFPSubstitutionModel(
      std::shared_ptr<const GeneticCode> gCode,
      const std::string& prefix = "AbstractDFP. ");

  AbstractDFPSubstitutionModel(const AbstractDFPSubstitutionModel& mod) :
    AbstractParameterAliasable(mod),
    AbstractSubstitutionModel(mod),
    gCode_(mod.gCode_),
    tr_(mod.tr_), trr_(mod.trr_), tvv_(mod.tvv_), trv_(mod.trv_), tsub_(mod.tsub_)
  {}

  AbstractDFPSubstitutionModel& operator=(const AbstractDFPSubstitutionModel& mod)
  {
    AbstractSubstitutionModel::operator=(mod);
    gCode_ = mod.gCode_;
    tr_ = mod.tr_;
    trr_ = mod.trr_;
    tvv_ =  mod.tvv_;
    trv_ = mod.trv_;
    tsub_ = mod.tsub_;

    return *this;
  }

  virtual ~AbstractDFPSubstitutionModel() {}

  AbstractDFPSubstitutionModel* clone() const override = 0;

public:
  std::shared_ptr<const GeneticCode> getGeneticCode() const override { return gCode_; }

  void fireParameterChanged(const ParameterList& parameters) override;

  using BranchModelInterface::getNumberOfStates;

  size_t getNumberOfStates() { return 64; }

  /**
   * @brief Calls  the multiplication by the specific codon-codon rate.
   */
  double getCodonsMulRate(size_t i, size_t j) const override;

protected:
  /**
   * @brief Method inherited from AbstractSubstitutionModel
   *
   * This method sets the rates to/from stop codons to zero and
   * set the generator given parameters.
   */
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTDFPSUBSTITUTIONMODEL_H
