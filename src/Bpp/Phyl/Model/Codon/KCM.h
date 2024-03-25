// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_KCM_H
#define BPP_PHYL_MODEL_CODON_KCM_H


#include "../AbstractBiblioSubstitutionModel.h"
#include "KroneckerCodonDistanceSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The general multiple substitution model for codons, from
 * Zaheri & al, 2014.
 *
 * @author Laurent Guéguen
 *
 * This model is built from one or several nucleotide substitution
 * models. It also allows distinct equilibrium frequencies between
 * codons. A multiplicative factor accounts for the selective
 * restraints at the amino acid level, depending on the synonymy of
 * the amino acids.
 *
 *
 * This model includes :
 *
 * - parameters of the nucleotide models,
 * - @f$\omega@f$ for the ratio of non-synonymous vs synonymous
 * substitution rates,
 * - parameters of the equilibrium frequencies

 * The codon frequencies @f$\pi_j@f$ are either observed or inferred.
 *
 * Reference:
 * -  Zaheri, M. and Dib, L. and Salamin, N. A Generalized
 * Mechanistic Codon Model,  Mol Biol Evol. 2014 Sep;31(9):2528-41.
 * doi: 10.1093/molbev/msu196. Epub 2014 Jun 23.
 */
class KCM :
  public AbstractBiblioSubstitutionModel,
  public virtual CodonReversibleSubstitutionModelInterface
{
private:
  std::unique_ptr<KroneckerCodonDistanceSubstitutionModel> pmodel_;
  bool oneModel_;

public:
  /**
   * @brief constructor.
   *
   * If onemod, a unique GTR model is used, otherwise three
   * different GTR models are used.
   */
  KCM(std::shared_ptr<const GeneticCode> gc, bool oneModel);

  KCM(const KCM& kcm);

  KCM& operator=(const KCM&);

  virtual ~KCM() {}

  KCM* clone() const override { return new KCM(*this); }

public:
  std::string getName() const override { return "KCM" + std::string(oneModel_ ? "7" : "19") + ".";}

  const SubstitutionModelInterface& substitutionModel() const override { return *pmodel_; }

  std::shared_ptr<const GeneticCode> getGeneticCode() const override { return pmodel_->getGeneticCode(); }

  double getCodonsMulRate(size_t i, size_t j) const override { return pmodel_->getCodonsMulRate(i, j); }

protected:
  SubstitutionModelInterface& substitutionModel_() override { return *pmodel_; }

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    throw NullPointerException("KCM::codonFrequencySet. No FrequencySet available for this model.");
  }

  bool hasCodonFrequencySet() const override
  {
    return false;
  }

  void setFreq(std::map<int, double>& frequencies) override
  {
    AbstractBiblioSubstitutionModel::setFreq(frequencies);
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_KCM_H
