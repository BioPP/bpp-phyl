// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_YNGP_M10_H
#define BPP_PHYL_MODEL_CODON_YNGP_M10_H

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

#include "YNGP_M.h"

namespace bpp
{
/**
 * @brief The Yang et al (2000) M10 substitution model for codons.
 * @author Laurent GuÃÂ©guen
 *
 * This model is a mixture of models as described in YN98 class, the
 * mixture being defined on the selection parameter omega to allow it
 * to vary among sites, following a mixture of a Beta distribution and
 * of a switched Gamma distribution.
 *
 * This model includes 5 parameters (@f$\kappa@f$, @f$ p @f$ and
 * @f$q@f$ of the @f$ Beta(p,q) @f$ distribution, @f$ \alpha @f$ and
 * @f$\beta@f$ of the @f$ 1 + Gamma(\alpha,\beta) @f$
 * distribution,@f$p0@f$ the weight of the Beta distribution. The
 * codon frequencies @f$ \pi_j @f$ are either observed or inferred.
 *
 * References:
 *
 * Yang, Z., R. Nielsen, N. Goldman, and A.-M. K. Pedersen (2000)
 * Genetics 155:431-449.
 *
 */
class YNGP_M10 :
  public YNGP_M
{
private:
  unsigned int nBeta_, nGamma_;

public:
  /**
   * @brief Constructor that requires the number of classes of the
   * BetaDiscreteDistribution and the GammaDiscreteDistribution.
   */
  YNGP_M10(
      std::shared_ptr<const GeneticCode> gc,
      std::unique_ptr<CodonFrequencySetInterface> codonFreqs,
      unsigned int nbBeta,
      unsigned int nbGamma);

  YNGP_M10* clone() const override { return new YNGP_M10(*this); }

  YNGP_M10(const YNGP_M10& mod2) :
    AbstractParameterAliasable(mod2),
    AbstractWrappedModel(mod2),
    AbstractWrappedTransitionModel(mod2),
    AbstractTotallyWrappedTransitionModel(mod2),
    AbstractBiblioTransitionModel(mod2),
    YNGP_M(mod2),
    nBeta_(mod2.nBeta_),
    nGamma_(mod2.nGamma_)
  {}

  YNGP_M10& operator=(const YNGP_M10& mod2)
  {
    YNGP_M::operator=(mod2);
    nBeta_ = mod2.nBeta_;
    nGamma_ = mod2.nGamma_;
    return *this;
  }

protected:
  void updateMatrices_() override;

public:
  std::string getName() const override { return "YNGP_M10"; }

  unsigned int getNBeta() const
  {
    return nBeta_;
  }

  unsigned int getNGamma() const
  {
    return nGamma_;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_YNGP_M10_H
