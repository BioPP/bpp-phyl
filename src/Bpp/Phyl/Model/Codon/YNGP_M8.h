// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_YNGP_M8_H
#define BPP_PHYL_MODEL_CODON_YNGP_M8_H

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

#include "YNGP_M.h"

namespace bpp
{
/**
 * @brief The Yang et al (2000) M8 substitution model for codons.
 * @author Laurent GuÃÂ©guen
 *
 * This model is a mixture of models as described in YN98 class, the
 * mixture being defined on the selection parameter oomega to allow it
 * to vary among sites, following a mixture of a Beta distribution and
 * of another value above 1.
 *
 * This model includes 5 parameters (@f$\kappa@f$, @f$ p @f$ and
 * @f$q@f$ of the @f$ Beta(p,q) @f$ distribution, @f$p0@f$ the weight
 * of the Beta distribution and @f$\omega @f$ the selection parameter
 * above 1 (with weight @f$ 1-p0 @f$)). The codon frequencies @f$
 * \pi_j @f$ are either observed or inferred.
 *
 * In option, the model YNGP_M8a is similar with @f$\omega=1 @f$ fixed
 * (and then only 4 parameters).
 *
 * References:
 *
 * Yang, Z., R. Nielsen, N. Goldman, and A.-M. K. Pedersen (2000)
 * Genetics 155:431-449.
 *
 */
class YNGP_M8 :
  public YNGP_M
{
private:
  /**
   * @brief If parameter omega=1
   */
  bool neutral_;

public:
  /**
   * @brief Constructor that requires the number of classes of the
   * BetaDiscreteDistribution.
   */
  YNGP_M8(
      std::shared_ptr<const GeneticCode> gc,
      std::unique_ptr<CodonFrequencySetInterface> codonFreqs,
      unsigned int nbclass,
      bool neutral = false);

  YNGP_M8* clone() const override { return new YNGP_M8(*this); }

  YNGP_M8(const YNGP_M8& mod2) :
    AbstractParameterAliasable(mod2),
    AbstractWrappedModel(mod2),
    AbstractWrappedTransitionModel(mod2),
    AbstractTotallyWrappedTransitionModel(mod2),
    AbstractBiblioTransitionModel(mod2),
    YNGP_M(mod2),
    neutral_(mod2.neutral_)
  {}

  YNGP_M8& operator=(const YNGP_M8& mod2)
  {
    YNGP_M::operator=(mod2);
    neutral_ = mod2.neutral_;
    return *this;
  }

  std::string getName() const override { return neutral_ ? "YNGP_M8a" : "YNGP_M8"; }

protected:
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_YNGP_M8_H
