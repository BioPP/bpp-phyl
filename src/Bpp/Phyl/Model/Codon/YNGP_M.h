// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_YNGP_M_H
#define BPP_PHYL_MODEL_CODON_YNGP_M_H

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

#include "../AbstractBiblioMixedTransitionModel.h"
#include "../FrequencySet/CodonFrequencySet.h"
#include "../MixtureOfASubstitutionModel.h"

namespace bpp
{
/**
 * @brief Abstract generic class for The Yang et al (2000) M
 * substitution models for codons. al (2004).
 * @author Laurent Guéguen
 *
 * References:
 *
 * Yang, Z., R. Nielsen, N. Goldman, and A.-M. K. Pedersen (2000)
 * Genetics 155:431-449.
 *
 * Wong, W. S. W., Z. Yang, N. Goldman, and R. Nielsen. (2004)
 * Genetics 168:1041--1051.
 */
class YNGP_M :
  public AbstractBiblioMixedTransitionModel,
  public virtual TransitionModelInterface
{
protected:
  /**
   * redefined mixed model pointer
   */
  const MixtureOfASubstitutionModel* mixedSubModelPtr_;

  /**
   * @brief indexes of 2 codons states between which the substitution is
   * synonymous, to set a basis to the homogeneization of the rates.
   */
  size_t synfrom_, synto_;

public:
  YNGP_M(const std::string& name) :
    AbstractBiblioMixedTransitionModel(name),
    mixedSubModelPtr_(),
    synfrom_(),
    synto_()
  {}

  YNGP_M(const YNGP_M& mod2) :
    AbstractBiblioMixedTransitionModel(mod2),
    mixedSubModelPtr_(),
    synfrom_(mod2.synfrom_),
    synto_(mod2.synto_)
  {
    mixedSubModelPtr_ = dynamic_cast<const MixtureOfASubstitutionModel*>(&mixedModel());
  }

  virtual YNGP_M* clone() const override = 0;

  YNGP_M& operator=(const YNGP_M& mod2)
  {
    AbstractBiblioMixedTransitionModel::operator=(mod2);

    mixedSubModelPtr_ = dynamic_cast<const MixtureOfASubstitutionModel*>(&mixedModel());

    synfrom_ = mod2.synfrom_;
    synto_ = mod2.synto_;

    return *this;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_YNGP_M_H
