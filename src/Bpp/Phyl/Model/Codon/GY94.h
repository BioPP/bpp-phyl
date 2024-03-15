// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_GY94_H
#define BPP_PHYL_MODEL_CODON_GY94_H

#include <Bpp/Seq/AlphabetIndex/GranthamAAChemicalDistance.h>

#include "../AbstractBiblioSubstitutionModel.h"
#include "CodonDistanceFrequenciesSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The Goldman and Yang (1994) substitution model for codons.
 * @author Laurent GuÃÂ©guen
 *
 * This model has one rate of transitions and one rate of
 * transversion. It also allows distinct equilibrium frequencies
 * between codons. A multiplicative factor accounts for the selective
 * restraints at the amino acid level. This factor applies on the
 * distance @f$d@f$ between amino acids given by Grantham (1974).
 *
 * For codons @f$i=i_1i_2i_3@f$ and @f$j=j_1j_2j_3@f$, the generator
 * term @f$Q_{ij} (i \neq j)@f$ is:
 *
 * 0 if 2 or 3 of the pair @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ are
 * different.
 *
 * @f$\mu \pi_j \exp(-d_{aa_i,aa_j}/V)@f$ if exactly 1 of the pairs
 * @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ is different, and that
 * difference is a transversion.
 *
 * @f$\mu \kappa \pi_j \exp(-d_{aa_i,aa_j}/V)@f$ if exactly 1 of the
 * pairs @f$(i_1,j_1) (i_2,j_2) (i_3,j_3) @f$ is different, and that
 * difference is a transition.
 *
 * @f$\mu@f$ is a normalization factor.
 *
 * This model includes 2 parameters (@f$\kappa@f$ and @f$V@f$). The
 * codon frequencies @f$\pi_j@f$ are either observed or infered.
 *
 * Reference:
 * - Goldman N. and Yang Z. (1994), _Molecular Biology And Evolution_ 11(5) 725--736.
 */
class GY94 :
  public AbstractBiblioSubstitutionModel,
  public virtual ReversibleSubstitutionModelInterface
{

private:

  std::shared_ptr<const GranthamAAChemicalDistance> gacd_;
  std::unique_ptr<CodonDistanceFrequenciesSubstitutionModel> pmodel_;

public:

  GY94(
      std::shared_ptr<const GeneticCode> gc,
      std::unique_ptr<CodonFrequencySetInterface> codonFreqs);

  virtual ~GY94();

  GY94(const GY94& gy94);

  GY94& operator=(const GY94& gy94);

  GY94* clone() const override { return new GY94(*this); }

public:

  std::string getName() const override { return "GY94"; }

  const SubstitutionModelInterface& substitutionModel() const override { return *pmodel_; }

  std::shared_ptr<const GeneticCode> getGeneticCode() const { return pmodel_->getGeneticCode(); }

  double getCodonsMulRate(size_t i, size_t j) const { return pmodel_->getCodonsMulRate(i, j); }

protected:
  
  SubstitutionModelInterface& substitutionModel_() override { return *pmodel_; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_GY94_H
