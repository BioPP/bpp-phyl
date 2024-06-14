// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_RELAX_H
#define BPP_PHYL_MODEL_CODON_RELAX_H

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

#include "YNGP_M.h"

namespace bpp
{
/**
 * @brief The RELAX (2014) branch-site model for codons
 *
 * @author Keren Halabi
 *
 * This model consists of a mixture of YN98 models, and allows for
 * mixture in two levels: the site level and the branch level.
 *
 * @f$\omega_0 = omega_1 * p < 1 @f$ (with probability @f$p_0 @f$)
 *
 * @f$\omega_1 <= 1 @f$ (with probability @f$p_1 @f$)
 *
 * @f$\omega_2 > 1 @f$ (with probability @f$1-p_1-p_0 @f$)
 *
 * Each omega is raised to the power of a selection intensity
 * parameter k.
 *
 * Used alone, this model is over-parameterized, and should be used in
 * a branch heterogeneous modeling with shared parameters with other
 * RELAX models.
 *
 * References:
 *
 * Joel O., et al. "RELAX: detecting relaxed selection in a phylogenetic framework."
 * Molecular biology and evolution 32.3 (2014): 820-832.Ã¢ÂÂ
 */
class RELAX :
  public YNGP_M
{
public:
  RELAX(
      std::shared_ptr<const GeneticCode> gc,
      std::unique_ptr<CodonFrequencySetInterface> codonFreqs);

  RELAX* clone() const override { return new RELAX(*this); }

protected:
  void updateMatrices_() override;

public:
  std::string getName() const override { return "RELAX"; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_RELAX_H
