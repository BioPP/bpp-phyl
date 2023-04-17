//
// File: YNGP_M7.h
// Authors:
//   Laurent Gueguen
// Created: 2010-05-08 00:00:00
//

/*
  Copyright or ÃÂ© or Copr. CNRS, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

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
