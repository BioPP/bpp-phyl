// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_ABSTRACTCODONAARATESUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_ABSTRACTCODONAARATESUBSTITUTIONMODEL_H


#include "../Protein/ProteinSubstitutionModel.h"
#include "CodonSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Abstract class for modelling of non-synonymous and
 *  synonymous substitution rates in codon models, given an amino acid
 *  rate matrix (from a shared_ptr model).
 *
 * @author Laurent GuÃÂ©guen
 *
 * From the generator @f$g@f$ between amino-acids, the non-synonymous
 *  rate is multiplied with, if the coded amino-acids are @f$x@f$ and
 *  @f$y@f$, @f$\beta*g(x,y)@f$ with positive parameter \c "beta".
 *
 * If paramSynRate is true, the synonymous substitution rate is
 *  multiplied with @f$\gamma@f$ (with optional positive parameter \c
 *  "gamma"), else it is multiplied with 1.
 *
 *
 */

class AbstractCodonAARateSubstitutionModel :
  public virtual CoreCodonSubstitutionModelInterface,
  public virtual AbstractParameterAliasable
{
private:
  std::shared_ptr<ProteinSubstitutionModelInterface> pAAmodel_;

  std::shared_ptr<const GeneticCode> pgencode_;

  double beta_;

  double gamma_;

  std::shared_ptr<const StateMapInterface> stateMap_;

public:
  /**
   * @brief Build a new AbstractCodonAARateSubstitutionModel object from
   *  a pointer to NucleotideSubstitutionModel.
   *
   * @param pmodel shared_ptr to an amino_acid generator
   * @param pgencode the genetic code
   * @param prefix the Namespace
   * @param paramSynRate is true iff synonymous rate is parameterised
   *       (default=false).
   */
  AbstractCodonAARateSubstitutionModel(
    std::shared_ptr<ProteinSubstitutionModelInterface> pmodel,
    std::shared_ptr<const GeneticCode> pgencode,
    const std::string& prefix,
    bool paramSynRate = false);


  AbstractCodonAARateSubstitutionModel(const AbstractCodonAARateSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    pAAmodel_(model.pAAmodel_),
    pgencode_(model.pgencode_),
    beta_(model.beta_),
    gamma_(model.gamma_),
    stateMap_(model.stateMap_)
  {}

  AbstractCodonAARateSubstitutionModel& operator=(
    const AbstractCodonAARateSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    pAAmodel_ = model.pAAmodel_;
    pgencode_ = model.pgencode_;
    beta_ = model.beta_;
    gamma_ = model.gamma_;
    stateMap_ = model.stateMap_;

    return *this;
  }

  AbstractCodonAARateSubstitutionModel* clone() const override
  {
    return new AbstractCodonAARateSubstitutionModel(*this);
  }

  virtual ~AbstractCodonAARateSubstitutionModel() {}

public:
  void fireParameterChanged(const ParameterList& parameters) override;

  double getCodonsMulRate(size_t i, size_t j) const override;

  void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix);
    pAAmodel_->setNamespace(prefix + pAAmodel_->getNamespace());
  }

  /*
   * @brief links to a new AA model
   *
   */
  void setAAModel(std::shared_ptr<ProteinSubstitutionModelInterface> model)
  {
    pAAmodel_ = model;
  }

  const ProteinSubstitutionModelInterface& aaModel() const
  {
    return *pAAmodel_;
  }

  const std::shared_ptr<ProteinSubstitutionModelInterface> getAAModel() const
  {
    return pAAmodel_;
  }

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    throw NullPointerException("AbstractCodonAARateSubstitutionModel::frequencySet. No associated FrequencySet.");
  }

  bool hasCodonFrequencySet() const override
  {
    return false;
  }

  void setFreq(std::map<int, double>& frequencies) override {}
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTCODONAARATESUBSTITUTIONMODEL_H
