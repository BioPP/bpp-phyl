// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_DFP07_H
#define BPP_PHYL_MODEL_CODON_DFP07_H

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

#include "../AbstractBiblioMixedTransitionModel.h"
#include "../FrequencySet/CodonFrequencySet.h"
#include "../MixtureOfASubstitutionModel.h"
#include "../Protein/ProteinSubstitutionModel.h"
#include "CodonSameAARateSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Class for non-synonymous substitution models on codons with
 * parameterized equilibrium frequencies and nucleotidic models, with
 * allowed multiple substitutions as parameterized in DFP model, with
 * correction to mimic AA substitution rates from a given protein
 * substitution model.
 *
 * Reference: Adi Doron-Faigenboim, Tal Pupko, 2007, A Combined
 * Empirical and Mechanistic Codon Model, Molecular Biology and
 * Evolution, Volume 24, Issue 2, Pages 388Ã¢ÂÂ397,
 * https://doi.org/10.1093/molbev/msl175
 *
 * This class should be used with models which equilibrium
 * distribution is fixed, ans does not depend on the parameters.
 * Otherwise there may be problems of identifiability of the
 * parameters.
 *
 * See description in AbstractDFPSubstitutionModel and
 * AbstractCodonFrequenciesSubstitutionModel and
 * CodonSameAARateSubstitutionModel class.
 *
 *
 * The additional parameters to AbstractDFPSubstitutionModel and
 * AbstractCodonFrequenciesSubstitutionModel are the rates of
 * nonsynonymous over synonymous substitutions.
 *
 * This rates are a mixture of neutral model (omega=1) and negative
 * selection (omega<1).
 *
 * Parameter p0 is the proportion of negative "omega<1" sub-model.
 */
class DFP07 :
  public AbstractBiblioMixedTransitionModel
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
  DFP07(
      std::shared_ptr<const GeneticCode> gCode,
      std::unique_ptr<ProteinSubstitutionModelInterface> pAAmodel,
      std::unique_ptr<CodonFrequencySetInterface> codonFreqs);

  DFP07(const DFP07& mod2) :
    AbstractParameterAliasable(mod2),
    AbstractWrappedModel(mod2),
    AbstractWrappedTransitionModel(mod2),
    AbstractTotallyWrappedTransitionModel(mod2),
    AbstractBiblioTransitionModel(mod2),
    AbstractBiblioMixedTransitionModel(mod2),
    mixedSubModelPtr_(nullptr),
    synfrom_(mod2.synfrom_),
    synto_(mod2.synto_)
  {
    mixedSubModelPtr_ = dynamic_cast<const MixtureOfASubstitutionModel*>(&mixedModel());
  }

  virtual DFP07* clone() const override
  {
    return new DFP07(*this);
  }

  DFP07& operator=(const DFP07& mod2)
  {
    const auto& eq = AbstractBiblioMixedTransitionModel::operator=(mod2);

    synfrom_ = mod2.synfrom_;
    synto_ = mod2.synto_;

    mixedSubModelPtr_ = dynamic_cast<const MixtureOfASubstitutionModel*>(&eq.mixedModel());

    return *this;
  }

  const ProteinSubstitutionModelInterface& proteinModel() const
  {
    return dynamic_cast<const CodonSameAARateSubstitutionModel&>(nModel(0)).proteinModel();
  }

  std::string getName() const override { return "DFP07"; }

protected:
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_DFP07_H
