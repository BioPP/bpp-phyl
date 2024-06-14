// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_YN98_H
#define BPP_PHYL_MODEL_CODON_YN98_H


#include "../AbstractBiblioSubstitutionModel.h"
#include "CodonDistanceFrequenciesSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The Yang and Nielsen (1998) substitution model for codons.
 * @author Laurent Guéguen
 *
 * This model has one rate of transitions and one rate of
 * transversion. It also allows distinct equilibrium frequencies
 * between codons. A multiplicative factor accounts for the selective
 * restraints at the amino acid level, depending on the synonymy of
 * the amino acids.
 *
 * For codons @f$i=i_1i_2i_3@f$ and @f$j=j_1j_2j_3@f$, the generator
 * term @f$Q_{ij} (i \neq j)@f$ is:
 *
 * 0 if 2 or 3 of the pair @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ are different.
 *
 * @f$\mu \pi_j \omega@f$ if exactly 1 of the pairs @f$(i_1,j_1)
 * (i_2,j_2) (i_3,j_3) @f$ is different, that difference is a
 * transversion and amino acids coded by i and j are different.
 *
 * @f$\mu \pi_j @f$ if exactly 1 of the pairs @f$(i_1,j_1) (i_2,j_2)
 * (i_3,j_3) @f$ is different, that difference is a transversion and
 * amino acids coded by i and j are the same.
 *
 * @f$\mu \kappa \pi_j \omega@f$ if exactly 1 of the pairs
 * @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ is different, that difference
 * is a transition and amino acids coded by i and j are different.
 *
 * @f$\mu \kappa \pi_j @f$ if exactly 1 of the pairs @f$(i_1,j_1)
 * (i_2,j_2) (i_3,j_3) @f$ is different, that difference is a
 * transition and amino acids coded by @f$i@f$ and @f$j@f$ are the
 * same.
 *
 * @f$\mu@f$ is a normalization factor.
 *
 * This model includes 2 parameters (@f$\kappa@f$ and @f$\omega@f$).
 * The codon frequencies @f$\pi_j@f$ are either observed or inferred.
 *
 * Reference:
 * -  Yang Z. and Nielsen R. (1998), _Journal of Molecular Evolution_ 46:409--418.
 */
class YN98 :
  public AbstractBiblioSubstitutionModel,
  public virtual CodonReversibleSubstitutionModelInterface
{
private:
  std::unique_ptr<CodonDistanceFrequenciesSubstitutionModel> pmodel_;

public:
  YN98(
      std::shared_ptr<const GeneticCode> gc,
      std::unique_ptr<CodonFrequencySetInterface> codonFreqs);

  YN98(const YN98& yn98);

  YN98& operator=(const YN98&);

  virtual ~YN98() {}

  YN98* clone() const override { return new YN98(*this); }

public:
  std::string getName() const override { return "YN98"; }

  const SubstitutionModelInterface& substitutionModel() const override { return *pmodel_; }

  std::shared_ptr<const GeneticCode> getGeneticCode() const override
  {
    return pmodel_->getGeneticCode();
  }

  const FrequencySetInterface& frequencySet() const override
  {
    return AbstractWrappedModel::frequencySet();
  }

  double getCodonsMulRate(size_t i, size_t j) const override
  {
    return pmodel_->getCodonsMulRate(i, j);
  }

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

  using AbstractBiblioSubstitutionModel::frequencySet;

protected:
  SubstitutionModelInterface& substitutionModel_() override { return *pmodel_; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_YN98_H
