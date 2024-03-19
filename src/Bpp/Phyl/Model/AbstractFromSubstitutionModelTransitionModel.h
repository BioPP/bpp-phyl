// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_ABSTRACTFROMSUBSTITUTIONMODELTRANSITIONMODEL_H
#define BPP_PHYL_MODEL_ABSTRACTFROMSUBSTITUTIONMODELTRANSITIONMODEL_H


#include "AbstractWrappedModel.h"

namespace bpp
{
/**
 * @brief Virtual class of a Transition Model related to a given
 * SubstitutionModel.
 *
 * It has the same parameters as the SubModel.
 */
class AbstractFromSubstitutionModelTransitionModel :
  public virtual AbstractWrappedTransitionModel
{
protected:
  /**
   * @brief The related model.
   */
  std::unique_ptr<SubstitutionModelInterface> subModel_;

  /**
   * The number of states
   */
  size_t size_;

  /**
   * @brief These ones are for bookkeeping:
   */
  mutable RowMatrix<double> pij_t;
  mutable RowMatrix<double> dpij_t;
  mutable RowMatrix<double> d2pij_t;

  std::string nestedPrefix_;

public:
  AbstractFromSubstitutionModelTransitionModel(
      std::unique_ptr<SubstitutionModelInterface> subModel,
      const std::string& prefix);

  AbstractFromSubstitutionModelTransitionModel(const AbstractFromSubstitutionModelTransitionModel& fmsm);

  AbstractFromSubstitutionModelTransitionModel& operator=(const AbstractFromSubstitutionModelTransitionModel& fmsm);

  virtual ~AbstractFromSubstitutionModelTransitionModel() {}

public:
  const SubstitutionModelInterface& substitutionModel() const
  {
    return *subModel_;
  }

  const TransitionModelInterface& transitionModel() const override
  {
    return *subModel_;
  }

  const BranchModelInterface& model() const override
  {
    return *subModel_;
  }

  bool computeFrequencies() const override
  {
    return subModel_->computeFrequencies();
  }

  /**
   * @return Set if equilibrium frequencies should be computed from
   * the generator
   */
  void computeFrequencies(bool yn) override
  {
    subModel_->computeFrequencies(yn);
  }

  /**
   * @}
   */

protected:
  Vdouble& getFrequencies_() override
  {
    return subModel_->getFrequencies_();
  }

  SubstitutionModelInterface& substitutionModel_()
  {
    return *subModel_;
  }


  TransitionModelInterface& transitionModel_() override
  {
    return *subModel_;
  }

  BranchModelInterface& model_()
  {
    return *subModel_;
  }

public:
  virtual void addRateParameter() override
  {
    model_().addRateParameter();
    addParameter_(new Parameter(getNamespace() + "rate", model().getRate(), Parameter::R_PLUS_STAR));
  }

  virtual void fireParameterChanged(const ParameterList& parameters) override
  {
    AbstractParameterAliasable::fireParameterChanged(parameters);
    model_().matchParametersValues(parameters);
  }

  virtual void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix + nestedPrefix_);
    model_().setNamespace(prefix + nestedPrefix_);
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ABSTRACTFROMSUBSTITUTIONMODELTRANSITIONMODEL_H
