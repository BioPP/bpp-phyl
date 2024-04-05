// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_MG94_H
#define BPP_PHYL_MODEL_CODON_MG94_H


#include "../AbstractBiblioSubstitutionModel.h"
#include "CodonDistancePhaseFrequenciesSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The Muse and Gaut (1994) substitution model for codons.
 * @author Laurent GuÃÂ©guen
 *
 * This model has one ratio @f$\rho@f$ of synonymous substitution rate
 * over non-synonymous substitution rate. It allows distinct
 * equilibrium frequencies between nucleotides.
 *
 * For codons @f$i=i_1i_2i_3@f$ and @f$j=j_1j_2j_3@f$, the generator term
 * @f$Q_{ij} (i \neq j)@f$ is:
 *
 * 0 if 2 or 3 of the pair @f$(i_1,j_1)(i_2,j_2) (i_3,j_3) @f$ are different.
 *
 * @f$\mu \rho \pi_{j_k} @f$  if exactly 1 of the pairs
 * @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ is different (@f$k@f$), and that
 * difference is non-synonymous.
 *
 * @f$\mu \pi_{j_k} @f$  if exactly 1 of the pairs
 * @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ is different (@f$k@f$), and that
 * difference is synonymous.
 *
 * @f$\rho@f$ corresponds to ratio @f$\frac{\beta}{\alpha}@f$ in
 * original publication.
 *
 * @f$\mu@f$ is a normalization factor.
 *
 * This model includes one parameter (@f$\rho@f$). The codon
 * frequencies @f$\pi_j@f$ are either observed or inferred.
 *
 * Reference:
 * - Muse S.V. and Gaut B.S. (1994), Molecular_ Biology And Evolution_ 11(5) 715--724.
 */
class MG94 :
  public AbstractBiblioSubstitutionModel,
  public virtual CodonReversibleSubstitutionModelInterface
{
private:
  std::unique_ptr<CodonDistancePhaseFrequenciesSubstitutionModel> pmodel_;

public:
  MG94(
      std::shared_ptr<const GeneticCode> gc,
      std::unique_ptr<CodonFrequencySetInterface> codonFreqs);

  MG94(const MG94& mg94);

  MG94& operator=(const MG94& mg94);

  virtual ~MG94();

  MG94* clone() const override { return new MG94(*this); }

public:
  std::string getName() const override { return "MG94"; }

  const SubstitutionModelInterface& substitutionModel() const override { return *pmodel_; }

  std::shared_ptr<const GeneticCode> getGeneticCode() const override { return pmodel_->getGeneticCode(); }

  double getCodonsMulRate(size_t i, size_t j) const override { return pmodel_->getCodonsMulRate(i, j); }

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    return pmodel_->codonFrequencySet();
  }

  bool hasCodonFrequencySet() const override
  {
    return pmodel_->hasCodonFrequencySet();
  }

  void setFreq(std::map<int, double>& frequencies) override
  {
    AbstractBiblioSubstitutionModel::setFreq(frequencies);
  }

protected:
  SubstitutionModelInterface& substitutionModel_() override { return *pmodel_; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_MG94_H
