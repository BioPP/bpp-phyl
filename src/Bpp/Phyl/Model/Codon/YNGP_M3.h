// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_YNGP_M3_H
#define BPP_PHYL_MODEL_CODON_YNGP_M3_H

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

#include "YNGP_M.h"

namespace bpp
{
/**
 * @brief The Yang et al (2000) M3 substitution model for codons.
 * @author Laurent Guéguen
 *
 * This model is a mixture of models as described in YN98 class, the
 * mixture being defined on the selection parameter to allow it to
 * vary among sites. There are $K$ selection parameters @f$ \omega_0 <
 * ... \omega_{K-1} @f$, with their respective probabilities @f$ p_0,
 * ..., p_{K-1} @f$ with @f$ p_0+p_1+...+p_{K-1}=1@f$. To garantee
 * that the @f$\omega_i@f$ are in increasing order, we define
 * @f$\delta_i=\omega_i - \omega_{i-1}@f$.
 *
 * This model includes 2*K parameters (@f$\kappa@f$, relative
 * probabilities @f$ theta1, theta2, ..., thetaK-1 @f$ and @f$omega0,
 * delta1, deltaK-1@f$). The codon frequencies @f$\pi_j@f$ are either
 * observed or infered.
 *
 * References:
 *
 * Yang, Z., R. Nielsen, N. Goldman, and A.-M. K. Pedersen (2000)
 * Genetics 155:431-449.
 *
 */
class YNGP_M3 :
  public YNGP_M
{
public:
  YNGP_M3(
      std::shared_ptr<const GeneticCode> gc,
      std::unique_ptr<CodonFrequencySetInterface> codonFreqs,
      unsigned int nclass = 3);

  YNGP_M3* clone() const override { return new YNGP_M3(*this); }

public:
  std::string getName() const override { return "YNGP_M3"; }

protected:
  void updateMatrices_() override;

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_YNGP_M3_H
