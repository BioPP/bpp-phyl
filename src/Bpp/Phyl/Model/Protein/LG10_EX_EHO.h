// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_PROTEIN_LG10_EX_EHO_H
#define BPP_PHYL_MODEL_PROTEIN_LG10_EX_EHO_H


#include "../AbstractBiblioMixedTransitionModel.h"
#include "../AbstractSubstitutionModel.h"
#include "ProteinSubstitutionModel.h"

using namespace std;

namespace bpp
{
/**
 * @brief The Le and Gascuel (2010) EX_EHO substitution model for proteins.
 * @author Mathieu Groussin
 *
 * This model is a mixture of six models. It combines two previously published models, EX2 and EHO.
 * Sites are classified as: exposed & extended, buried & extended, exposed & alpha-helix
 * buried & alpha-helix, exposed & other, or buried & other.
 * The submodels are called BUR_EXT, BUR_HEL, BUR_OTH, EXP_EXT, EXP_HEL and EXP_OTH.
 *
 * The model includes 10 parameters :
 *
 * - relrate1 to relrate5 are the relative rates of the submodels;
 * - relproba1 to relproba5 are the relative proportions of the submodels;
 *
 * Important: See the relation between these parameters and the
 * rates and probabilities of the models in the description of
 * MixtureOfSubstitutionModels class.
 *
 * Reference:
 *
 * Le S.Q., Gascuel O. (2010) Syst. Biol. 59(3):277Ã¢ÂÂ287
 */
class LG10_EX_EHO :
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
    EmbeddedModel* clone() const { return new EmbeddedModel(*this); }
    string getName() const { return name_;}
    double getProportion() const { return proportion_;}
  };

public:
  /**
   * @brief Build a EX_EHO model, with original equilibrium frequencies, probabilities and rates.
   *
   * @param alpha A proteic alphabet.
   *
   */

  LG10_EX_EHO(std::shared_ptr<const ProteicAlphabet> alpha);

  LG10_EX_EHO* clone() const { return new LG10_EX_EHO(*this); }

  std::string getName() const { return "LG10_EX_EHO"; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_LG10_EX_EHO_H
