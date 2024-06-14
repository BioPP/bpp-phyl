// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_PROTEIN_LLG08_UL2_H
#define BPP_PHYL_MODEL_PROTEIN_LLG08_UL2_H


#include "../AbstractBiblioMixedTransitionModel.h"
#include "../AbstractSubstitutionModel.h"
#include "../MixtureOfSubstitutionModels.h"
#include "ProteinSubstitutionModel.h"

using namespace std;

namespace bpp
{
/**
 * @brief The Le et al  (2008) UL2 substitution model for proteins.
 * @author Laurent Guéguen
 *
 * This model is a mixture of two models built by an unsupervised
 * method (see ref). The submodels are called M1 & M2.
 *
 *
 * This model includes 2 parameters :
 *
 * - relrate1 is the relative rate of model of M1;
 * - relproba1 is the proportion of model M1;
 *
 * Important: See the relation between these parameters and the
 * rates and probabilities of the models in the description of
 * MixtureOfSubstitutionModels class.
 *
 * Reference:
 *
 * Le S.Q., Lartillot N., Gascuel O. (2008) Phil. Trans. R. Soc. B 363:3965--3976.
 */

class LLG08_UL2 :
  public AbstractBiblioMixedTransitionModel
{
public:
  class EmbeddedModel :
    public AbstractReversibleProteinSubstitutionModel
  {
private:
    double proportion_;
    string name_;

public:
    EmbeddedModel(std::shared_ptr<const ProteicAlphabet> alpha, string name);
    EmbeddedModel* clone() const override { return new EmbeddedModel(*this); }
    string getName() const override { return name_; }
    double getProportion() const { return proportion_;}
  };

public:
  /**
   * @brief Build a  UL2 model, with original equilibrium frequencies, probabilities and rates.
   *
   * @param alpha A proteic alphabet.
   */
  LLG08_UL2(std::shared_ptr<const ProteicAlphabet> alpha);

  LLG08_UL2* clone() const override { return new LLG08_UL2(*this); }

  std::string getName() const override { return "LLG08_UL2"; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_LLG08_UL2_H
