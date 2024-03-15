// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_PROTEIN_LGL08_CAT_H
#define BPP_PHYL_MODEL_PROTEIN_LGL08_CAT_H


#include "../AbstractBiblioMixedTransitionModel.h"
#include "../AbstractSubstitutionModel.h"
#include "../MixtureOfSubstitutionModels.h"
#include "ProteinSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The Le et al  (2008) CAT substitution model for proteins.
 * @author Mathieu Groussin
 *
 * This model is a mixture of N profiles empirically built with an EM algorithm
 * (see ref). The submodels are called C1, C2, ..., CN. For each model, exchangeabilities are equal (F81 model).
 *
 *
 * This model includes 2N-2 parameters :
 *
 * - relrate1, ..., relrate(N-1) are the relative rates of model C1, ..., C(N-1);
 * - relproba1, ..., relproba(N-1) are the relative proportions of model C1, ..., C(N-1);
 *
 * Important: See the relation between these parameters and the
 * rates and probabilities of the models in the description of
 * MixtureOfSubstitutionModels class.
 *
 * Reference:
 *
 * Le S.Q., Gascuel O. & Lartillot N. (2008) Bioinformatics 24:2317Ã¢ÂÂ2323.
 */

class LGL08_CAT :
  public AbstractBiblioMixedTransitionModel
{
public:
  class EmbeddedModel :
    public AbstractReversibleProteinSubstitutionModel
  {
private:
    double proportion_;
    std::string name_;

public:
    EmbeddedModel(
	std::shared_ptr<const ProteicAlphabet> alpha,
	std::string name,
	unsigned int nbCat = 10);

    virtual ~EmbeddedModel() {}

    EmbeddedModel* clone() const override { return new EmbeddedModel(*this); }

    std::string getName() const override { return name_; }

    double getProportion() const { return proportion_; }
  };

public:
  /**
   * @brief Build a CAT model, with original equilibrium frequencies, probabilities and rates.
   *
   * @param alpha A proteic alphabet.
   * @param nbCat number of profiles
   *
   */
  LGL08_CAT(std::shared_ptr<const ProteicAlphabet> alpha, unsigned int nbCat = 10);

  LGL08_CAT* clone() const override { return new LGL08_CAT(*this); }

  LGL08_CAT(const LGL08_CAT& mod2) :
    AbstractParameterAliasable(mod2),
    AbstractWrappedModel(mod2),
    AbstractWrappedTransitionModel(mod2),
    AbstractTotallyWrappedTransitionModel(mod2),
    AbstractBiblioTransitionModel(mod2),
    AbstractBiblioMixedTransitionModel(mod2)
  {}

  LGL08_CAT& operator=(const LGL08_CAT& mod2)
  {
    AbstractBiblioMixedTransitionModel::operator=(mod2);
    return *this;
  }

  uint getNumberOfCategories() const
  {
    return static_cast<uint>(mixedModelPtr_->getNumberOfModels());
  }

  std::string getName() const override { return "LGL08_CAT";}
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_LGL08_CAT_H
