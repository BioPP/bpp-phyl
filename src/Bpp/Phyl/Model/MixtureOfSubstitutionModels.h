// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_MIXTUREOFSUBSTITUTIONMODELS_H
#define BPP_PHYL_MODEL_MIXTUREOFSUBSTITUTIONMODELS_H

#include <Bpp/Numeric/VectorTools.h>
#include <cstring> // C lib for string copy
#include <map>
#include <string>
#include <vector>

#include "MixtureOfTransitionModels.h"

namespace bpp
{
class MixtureOfSubstitutionModels :
  public MixtureOfTransitionModels
{
public:
  /**
   * @brief Constructor of a MixtureOfSubstitutionModels, where all
   * the models have rate 1 and equal probability.
   *
   * @param alpha pointer to the Alphabet
   * @param vpModel vector of pointers to SubstitutionModels. All the
   *   SubstitutionModels are owned by the instance.
   * @warning providing a vpModel with size 0 will generate a segmentation fault!
   */
  MixtureOfSubstitutionModels(
      std::shared_ptr<const Alphabet> alpha,
      std::vector< std::unique_ptr<TransitionModelInterface>>& vpModel) :
    AbstractParameterAliasable("Mixture."),
    AbstractTransitionModel(alpha, vpModel.size() ? vpModel[0]->getStateMap() : 0, "Mixture."),
    MixtureOfTransitionModels(alpha, vpModel)
  {
    // Check that all models are substitutionmodels
    for (const auto& model : modelsContainer_)
    {
      if (!dynamic_cast<const SubstitutionModelInterface*>(model.get()))
        throw Exception("MixtureOfSubstitutionModels can only be built with SubstitutionModels, not " + model->getName());
    }
  }

  /**
   * @brief Constructor of a MixtureOfSubstitutionModels.
   *
   * @param alpha pointer to the Alphabet
   * @param vpModel vector of pointers to SubstitutionModels. All the
   *   SubstitutionModels are owned by the instance.
   * @param vproba vector of the probabilities of the models
   * @param vrate vector of the rates of the models
   * @warning providing a vpModel with size 0 will generate a segmentation fault!
   *
   * See above the constraints on the rates and the probabilities of
   * the vectors.
   */
  MixtureOfSubstitutionModels(
      std::shared_ptr<const Alphabet> alpha,
      std::vector< std::unique_ptr<TransitionModelInterface>>& vpModel,
      Vdouble& vproba, Vdouble& vrate) :
    AbstractParameterAliasable("Mixture."),
    AbstractTransitionModel(alpha, vpModel.size() ? vpModel[0]->getStateMap() : 0, "Mixture."),
    MixtureOfTransitionModels(alpha, vpModel, vproba, vrate)
  {
    /*
     * Check that all models are substitutionmodels
     */

    for (const auto& model : modelsContainer_)
    {
      if (!dynamic_cast<const SubstitutionModelInterface*>(model.get()))
        throw Exception("MixtureOfSubstitutionModels can only be built with SubstitutionModels, not " + model->getName());
    }
  }


  MixtureOfSubstitutionModels(const MixtureOfSubstitutionModels& model) :
    AbstractParameterAliasable(model),
    AbstractTransitionModel(model),
    MixtureOfTransitionModels(model)
  {}


  MixtureOfSubstitutionModels& operator=(const MixtureOfSubstitutionModels& model)
  {
    MixtureOfTransitionModels::operator=(model);
    return *this;
  }

  MixtureOfSubstitutionModels* clone() const override { return new MixtureOfSubstitutionModels(*this); }

public:
  /**
   * @brief retrieve a pointer to the substitution model with the given name.
   */
  const SubstitutionModelInterface& subModel(const std::string& name) const
  {
    return dynamic_cast<const SubstitutionModelInterface&>(model(name));
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_MIXTUREOFSUBSTITUTIONMODELS_H
