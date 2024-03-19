// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_FROMMIXTURESUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_FROMMIXTURESUBSTITUTIONMODEL_H


#include "AbstractSubstitutionModel.h"
#include "AbstractWrappedModel.h"
#include "MixedTransitionModel.h"

namespace bpp
{
/**
 * @brief Model taken from a SubModel of a
 * Mixture of SubstitutionModels.
 *
 * It has the same parameters as the SubModel.
 */
class FromMixtureSubstitutionModel :
  public virtual AbstractTotallyWrappedSubstitutionModel
{
private:
  /**
   * @brief The subModel taken from the AbstractTotallyWrappedSubstitutionModel.
   *
   * This subModel is normalized, even if it is not in the mixture.
   */
  std::unique_ptr<SubstitutionModelInterface> subModel_;

  /**
   * @brief The name of the mixture model (for io purpose).
   */
  std::string mixtName_;

public:
  FromMixtureSubstitutionModel(
      const MixedTransitionModelInterface& mixedModel,
      const std::string& subModelName, const std::string& mixtDesc);

  FromMixtureSubstitutionModel(
      const MixedTransitionModelInterface& mixedModel,
      size_t subModelNumber,
      const std::string& mixtDesc);

  FromMixtureSubstitutionModel(const FromMixtureSubstitutionModel& fmsm);

  FromMixtureSubstitutionModel& operator=(const FromMixtureSubstitutionModel& fmsm);

  virtual ~FromMixtureSubstitutionModel() {}

  FromMixtureSubstitutionModel* clone() const override { return new FromMixtureSubstitutionModel(*this); }

public:
  const SubstitutionModelInterface& substitutionModel() const override
  {
    return *subModel_.get();
  }

protected:
  SubstitutionModelInterface& substitutionModel_() override
  {
    return *subModel_;
  }

public:
  /**
   * @
   * brief Methods to supersede AbstractSubstitutionModel methods.
   *
   * @{
   */

  /**
   * @brief Tells the model that a parameter value has changed.
   *
   * This updates the matrices consequently.
   */
  void fireParameterChanged(const ParameterList& parameters) override
  {
    model_().matchParametersValues(parameters);
  }

  virtual void setNamespace(const std::string& name) override
  {
    AbstractParameterAliasable::setNamespace(name);
    model_().setNamespace(name);
  }

  virtual void addRateParameter() override
  {
    model_().addRateParameter();
    addParameter_(new Parameter(getNamespace() + "rate", model().getRate(), Parameter::R_PLUS_STAR));
  }

  /*
   * @}
   */
  std::string getName() const override
  {
    size_t posp = mixtName_.find("(");
    return mixtName_.substr(0, posp) + "_" + model().getName() + mixtName_.substr(posp);
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_FROMMIXTURESUBSTITUTIONMODEL_H
