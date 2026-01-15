// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_MG94xREV_H
#define BPP_PHYL_MODEL_CODON_MG94xREV_H


#include "../AbstractBiblioSubstitutionModel.h"
#include "CodonDistancePhaseFrequenciesSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The Muse and Gaut (1994) substitution model for codons with GTR reversibility matrix, as in HyPhy.
 * @author Laurent GuÃÂ©guen
 *
 * This model has one rate @f$\alpha@f$ of synonymous substitution
 * rate and one rate @f$\beta" for non-synonymous substitution rate.
 * It is based on a GTR mutation exchangeability matrix ($\Theta$)
 * with codon specific equilibrium frequencies.
 *
 * For codons @f$i=i_1i_2i_3@f$ and @f$j=j_1j_2j_3@f$, the generator term
 * @f$Q_{ij} (i \neq j)@f$ is:
 *
 * 0 if 2 or 3 of the pair @f$(i_1,j_1)(i_2,j_2) (i_3,j_3) @f$ are different.
 *
 * @f$\mu \alpha \Theta_{i_k,j_k} \pi_{j_k} @f$ if exactly 1 of the
 * pairs @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ is different (@f$k@f$),
 * and that difference is synonymous.
 *
 * @f$\mu \rho \Theta_{i_k,j_k} \pi_{j_k} @f$ if exactly 1 of the pairs
 * @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ is different (@f$k@f$), and
 * that difference is non-synonymous.
 *
 * This model includes parameter (@f$\omega@f$) and 5 non equilibrium
 * parameters of GTR model. The codon frequencies @f$\pi_j@f$ are
 * either observed or inferred.
 *
 * @see GTR
 *
 * Reference:
 * - Kosakovsky Pond & Frost, 2005, MBE
 */

  class MG94xREV :
  public AbstractBiblioSubstitutionModel,
  public virtual CodonReversibleSubstitutionModelInterface
{
private:
  std::unique_ptr<CodonDistancePhaseFrequenciesSubstitutionModel> pmodel_;

public:
  MG94xREV(
      std::shared_ptr<const GeneticCode> gc,
      std::unique_ptr<CodonFrequencySetInterface> codonFreqs);

  MG94xREV(const MG94xREV& mg94);

  MG94xREV& operator=(const MG94xREV& mg94);

  virtual ~MG94xREV();

  MG94xREV* clone() const override { return new MG94xREV(*this); }

public:
  std::string getName() const override { return "MG94xREV"; }

  const SubstitutionModelInterface& substitutionModel() const override { return *pmodel_; }

  std::shared_ptr<const GeneticCode> getGeneticCode() const override { return pmodel_->getGeneticCode(); }

  const FrequencySetInterface& frequencySet() const override
  {
    return AbstractWrappedModel::frequencySet();
  }

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
#endif // BPP_PHYL_MODEL_CODON_MG94xREV_H
