// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_YNGP_M7_H
#define BPP_PHYL_MODEL_CODON_YNGP_M7_H

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

#include "YNGP_M.h"

namespace bpp
{
/**
 * @brief The Yang et al (2000) M7 substitution model for codons.
 * @author Laurent Guéguen
 *
 * This model is a mixture of models as described in YN98 class, the
 * mixture being defined on the selection parameter to allow it to
 * vary among sites, following a Beta distribution.
 *
 * This model includes 3 parameters (@f$\kappa@f$, @f$ p @f$ and
 * @f$q@f$) of the Beta distribution. The codon frequencies
 * @f$\pi_j@f$ are either observed or infered.
 *
 * References:
 *
 * Yang, Z., R. Nielsen, N. Goldman, and A.-M. K. Pedersen (2000)
 * Genetics 155:431-449.
 */
class YNGP_M7 :
  public YNGP_M
{
public:
  /**
   * @brief Constructor that requires the number of classes of the
   * BetaDiscreteDistribution.
   */
  YNGP_M7(
      std::shared_ptr<const GeneticCode> gc,
      std::unique_ptr<CodonFrequencySetInterface> codonFreqs,
      unsigned int nclass);

  YNGP_M7* clone() const override { return new YNGP_M7(*this); }

  YNGP_M7(const YNGP_M7& mod2) :
    AbstractParameterAliasable(mod2),
    AbstractWrappedModel(mod2),
    AbstractWrappedTransitionModel(mod2),
    AbstractTotallyWrappedTransitionModel(mod2),
    AbstractBiblioTransitionModel(mod2),
    YNGP_M(mod2)
  {}

  YNGP_M7& operator=(const YNGP_M7& mod2)
  {
    YNGP_M::operator=(mod2);
    return *this;
  }
  
  std::string getName() const override { return "YNGP_M7"; }

protected:
  
  void updateMatrices_() override;

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_YNGP_M7_H
