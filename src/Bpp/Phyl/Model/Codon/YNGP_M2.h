// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_YNGP_M2_H
#define BPP_PHYL_MODEL_CODON_YNGP_M2_H

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

#include "YNGP_M.h"

namespace bpp
{
/**
 * @brief The Yang et al (2000) M2 substitution model for codons, with
 * the more realistic modification in Wong & al (2004).
 * @author Laurent Guéguen
 *
 * This model is a mixture of models as described in YN98 class, the
 * mixture being defined on the selection parameter to allow it to
 * vary among sites. A site is either negatively selected @f$ 0 <
 * \omega_0 < 1 @f$ (with probability @f$p_0 @f$), or neutral (@f$
 * \omega_1 = 1 @f$) with probability @f$p_1 @f$, or positively
 * selected @f$ 1 < \omega_2 @f$ (with probability @f$p_2 @f$). This
 * model includes 5 parameters (@f$kappa@f$, @f$ theta1=p0,
 * theta2=\frac{p1}{p1+p2}, omega0 @f$ and @f$ omega2 @f$). The codon
 * frequencies @f$\pi_j@f$ are either observed or infered.
 *
 * References:
 *
 * Yang, Z., R. Nielsen, N. Goldman, and A.-M. K. Pedersen (2000)
 * Genetics 155:431-449.
 *
 * Wong, W. S. W., Z. Yang, N. Goldman, and R. Nielsen. (2004)
 * Genetics 168:1041--1051.
 */
class YNGP_M2 :
  public YNGP_M
{
public:
  YNGP_M2(
      std::shared_ptr<const GeneticCode> gc,
      std::unique_ptr<CodonFrequencySetInterface> codonFreqs);

  YNGP_M2* clone() const override { return new YNGP_M2(*this); }

public:
  std::string getName() const override { return "YNGP_M2"; }

protected:
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_YNGP_M2_H
