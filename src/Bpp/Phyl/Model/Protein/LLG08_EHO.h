// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_PROTEIN_LLG08_EHO_H
#define BPP_PHYL_MODEL_PROTEIN_LLG08_EHO_H


#include "../AbstractBiblioMixedTransitionModel.h"
#include "../AbstractSubstitutionModel.h"
#include "../MixtureOfTransitionModels.h"
#include "ProteinSubstitutionModel.h"

using namespace std;

namespace bpp
{
/**
 * @brief The Le et al  (2008) EH0 substitution model for proteins.
 * @author Laurent Guéguen
 *
 * This model is a mixture of three models corresponding to
 * extended/helix/other sites in proteins. The models are considered
 * in this order.
 *
 *
 * This model includes 4 parameters :
 *
 * - relrate1 is the relative rate of model of extended sites;
 * - relrate2 is the relative rate of helix sites;
 * - relproba1 is the proportion  of extended sites;
 * - relproba2 is the ratio of the proportions of helix sites over the sum of
 * the proportion of helix sites plus the proportion of other sites.
 *
 * Important: See the relation between these parameters and the rates
 * and probabilities of the models in the description of
 * MixtureOfSubstitutionModels class.
 *
 * Reference:
 *
 * Le S.Q., Lartillot N., Gascuel O. (2008) Phil. Trans. R. Soc. B 363:3965--3976.
 */
class LLG08_EHO :
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
    string getName() const override { return name_;}
    double getProportion() const { return proportion_;}
  };

public:
  /**
   * @brief Build a EH0 model, with original equilibrium frequencies, probabilities and rates.
   *
   * @param alpha A proteic alphabet.
   */
  LLG08_EHO(std::shared_ptr<const ProteicAlphabet> alpha);

  LLG08_EHO* clone() const override { return new LLG08_EHO(*this); }

  std::string getName() const override { return "LLG08_EHO"; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_LLG08_EHO_H
