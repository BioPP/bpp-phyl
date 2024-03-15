// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_MIXTUREOFASUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_MIXTUREOFASUBSTITUTIONMODEL_H

#include <Bpp/Numeric/VectorTools.h>
#include <cstring> // C lib for string copy
#include <map>
#include <string>
#include <vector>

#include "MixtureOfATransitionModel.h"

namespace bpp
{
class MixtureOfASubstitutionModel :
  public MixtureOfATransitionModel
{
public:
  /**
   * @brief Constructor of a MixtureOfASubstitutionModel, where all
   * the models have rate 1 and equal probability.
   *
   * @param alpha pointer to the Alphabet
   * @param model pointer to the SubstitutionModel that will be mixed
   * @param parametersDistributionsList list from parameters names towards discrete distributions to will define the mixtures.
   * @param ffrom   index of the starting codon that will be used to homogeneize the rates of the submodels
   * @param tto     index of the arriving codon that will be used to homogeneize the rates of the submodels
   *
   *   If ffrom and tto are not -1, for all submodels the transition
   *   rate ffrom->tto is the same. Otherwise, all submodels are
   *   normalized to have a substitution/time unit at equilibrium.
   */
  MixtureOfASubstitutionModel(
      std::shared_ptr<const Alphabet> alpha,
      std::unique_ptr<SubstitutionModelInterface> model,
      std::map<std::string, std::unique_ptr<DiscreteDistributionInterface>>& parametersDistributionsList,
      int ffrom = -1,
      int tto = -1) :
    AbstractParameterAliasable(model->getNamespace()),
    AbstractTransitionModel(alpha, model->getStateMap(), model->getNamespace()),
    MixtureOfATransitionModel(alpha, std::move(model), parametersDistributionsList, ffrom, tto)
  {}

  MixtureOfASubstitutionModel(const MixtureOfASubstitutionModel& model) :
    AbstractParameterAliasable(model),
    AbstractTransitionModel(model),
    MixtureOfATransitionModel(model)
  {}


  MixtureOfASubstitutionModel& operator=(const MixtureOfASubstitutionModel& model)
  {
    MixtureOfATransitionModel::operator=(model);
    return *this;
  }

  MixtureOfASubstitutionModel* clone() const override { return new MixtureOfASubstitutionModel(*this); }

protected:

  void updateMatrices_() override
  {
    MixtureOfATransitionModel::updateMatrices_();
    // setting the rates, if to_ & from_ are different from -1

    if (to_ >= 0 && from_ >= 0)
    {
      Vdouble vd;
      for (size_t j = 0; j < modelsContainer_.size(); ++j)
      {
        vd.push_back(1 / subNModel(j).Qij(static_cast<size_t>(from_), static_cast<size_t>(to_)));
      }
      setVRates(vd);
    }
  }

public:
  
  /**
   * @brief retrieve a pointer to the subsitution model with the given name.
   */
  const SubstitutionModelInterface& subModel(const std::string& name) const
  {
    return dynamic_cast<const SubstitutionModelInterface&>(model(name));
  }

  const SubstitutionModelInterface& subNModel(size_t i) const
  {
    return dynamic_cast<const SubstitutionModelInterface&>(nModel(i));
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_MIXTUREOFASUBSTITUTIONMODEL_H
