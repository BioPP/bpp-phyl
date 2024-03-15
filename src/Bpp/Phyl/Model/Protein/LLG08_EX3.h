// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_PROTEIN_LLG08_EX3_H
#define BPP_PHYL_MODEL_PROTEIN_LLG08_EX3_H


#include "../AbstractBiblioMixedTransitionModel.h"
#include "../AbstractSubstitutionModel.h"
#include "../MixtureOfSubstitutionModels.h"
#include "ProteinSubstitutionModel.h"

using namespace std;

namespace bpp
{
/**
 * @brief The Le et al  (2008) EX3 substitution model for proteins.
 * @author Laurent Guéguen
 *
 * This model is a mixture of three models corresponding to
 * buried/intermediate/highly exposed sites in proteins. The models
 * are considered in this order.
 *
 *
 * This model includes 4 parameters :
 *
 * - relrate1 is the relative rate of model of buried sites;
 * - relrate2 is the relative rate of intermediate sites;
 * - relproba1 is the proportion  of buried sites;
 * - relproba2 is the ratio of the proportions of intermediate sites
 * over the sum of the proportion of intermediate sites plus the
 * proportion of highly exposed sites.
 *
 * Important: See the relation between these parameters and the
 * rates and probabilities of the models in the description of
 * MixtureOfSubstitutionModels class.
 *
 * Reference:
 *
 * Le S.Q., Lartillot N., Gascuel O. (2008) Phil. Trans. R. Soc. B 363:3965--3976.
 */
class LLG08_EX3 :
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
   * @brief Build a  EX3 model, with original equilibrium frequencies, probabilities and rates.
   *
   * @param alpha A proteic alphabet.
   */
  LLG08_EX3(std::shared_ptr<const ProteicAlphabet> alpha);

  LLG08_EX3* clone() const override { return new LLG08_EX3(*this); }

  std::string getName() const override { return "LLG08_EX3"; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_LLG08_EX3_H
