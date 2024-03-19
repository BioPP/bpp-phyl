// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_YNGP_M1_H
#define BPP_PHYL_MODEL_CODON_YNGP_M1_H


#include "YNGP_M.h"

namespace bpp
{
/**
 * @brief The Yang et al (2000) M1 substitution model for codons, with
 * the more realistic modification in Wong & al (2004).
 * @author Laurent GuÃÂ©guen
 *
 * This model is a mixture of models as described in YN98 class, the
 * mixture being defined on the selection parameter to allow it to
 * vary among sites. A site is either negatively selected @f$ 0 <
 * \omega_0 < 1 @f$ (with probability @f$p_0 @f$), or neutral (@f$
 * \omega_1 = 1 @f$) with probability @f$1-p_0 @f$.
 *
 * The synonymous rates must be the same between both models, so the
 * overall rates of the models are modified to respect this constraint
 * and such that the mean rate of the mixed model equals one.
 *
 * This model includes 3 parameters (@f$\kappa@f$, @f$ p0 @f$ and
 * @f$\omega@f$). The codon frequencies @f$\pi_j@f$ are either
 * observed or infered.
 *
 * References:
 *
 * Yang, Z., R. Nielsen, N. Goldman, and A.-M. K. Pedersen (2000)
 * Genetics 155:431-449.
 *
 * Wong, W. S. W., Z. Yang, N. Goldman, and R. Nielsen. (2004)
 * Genetics 168:1041--1051.
 */
class YNGP_M1 :
  public YNGP_M
{
public:
  YNGP_M1(
      std::shared_ptr<const GeneticCode> gc,
      std::unique_ptr<CodonFrequencySetInterface> codonFreqs);

  YNGP_M1* clone() const override { return new YNGP_M1(*this); }

public:
  std::string getName() const override { return "YNGP_M1"; }

protected:
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_YNGP_M1_H
