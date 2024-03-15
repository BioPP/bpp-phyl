// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_REGISTERRATESSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_REGISTERRATESSUBSTITUTIONMODEL_H

#include <Bpp/Numeric/AbstractParameterAliasable.h>

#include "../Mapping/SubstitutionRegister.h"
#include "AbstractSubstitutionModel.h"
#include "AbstractWrappedModel.h"
#include "AnonymousSubstitutionModel.h"

namespace bpp
{
/**
 * @brief From a model, substitution rates are set into categories
 * following a given register. Each substitution of a category is then
 * multiplied by a rate parameter specific to this category.
 *
 * It has the same parameters as the SubModel.
 *
 * @see SubstitutionRegister
 */
class RegisterRatesSubstitutionModel :
  public AbstractWrappedSubstitutionModel,
  public AbstractSubstitutionModel
{
private:
  /**
   * @brief The related model.
   */
  std::unique_ptr<SubstitutionModelInterface> originalModel_;

  /**
   * For output
   */
  std::string registerName_;

  /**
   * @brief Vector of register state -> vector of from states ->
   * vector of to states (for acceleration purpose)
   */
  VVVuint vRegStates_;
  size_t nbTypes_;

  /**
   * @brief vector of the rates of the register types
   */
  Vdouble vRates_;

public:
  /**
   * @brief Constructor
   *
   * @param originalModel the substitution model used
   * @param reg the register in which the considered types of event are
   * used.
   * @param isNormalized says if model is normalized (default false)
   */
  RegisterRatesSubstitutionModel(
      std::unique_ptr<SubstitutionModelInterface> originalModel,
      const SubstitutionRegisterInterface& reg,
      bool isNormalized = false);


  RegisterRatesSubstitutionModel(const RegisterRatesSubstitutionModel& fmsm) :
    AbstractParameterAliasable(fmsm),
    AbstractWrappedModel(fmsm),
    AbstractWrappedTransitionModel(fmsm),
    AbstractWrappedSubstitutionModel(fmsm),
    AbstractSubstitutionModel(fmsm),
    originalModel_(fmsm.originalModel_->clone()),
    registerName_(fmsm.registerName_),
    vRegStates_(fmsm.vRegStates_),
    nbTypes_(fmsm.nbTypes_),
    vRates_(fmsm.vRates_)
  {}


  RegisterRatesSubstitutionModel& operator=(const RegisterRatesSubstitutionModel& fmsm)
  {
    AbstractWrappedSubstitutionModel::operator=(fmsm);
    AbstractSubstitutionModel::operator=(fmsm);
    
    originalModel_.reset(fmsm.originalModel_->clone());
    registerName_ = fmsm.registerName_;
    vRegStates_ = fmsm.vRegStates_;
    nbTypes_ = fmsm.nbTypes_;
    vRates_ = fmsm.vRates_;

    return *this;
  }

  virtual ~RegisterRatesSubstitutionModel() {}

  RegisterRatesSubstitutionModel* clone() const override
  {
    return new RegisterRatesSubstitutionModel(*this);
  }

public:
  /**
   * @brief From AbstractWrappedSubstitutionModel
   */
  const SubstitutionModelInterface& substitutionModel() const override
  {
    return *originalModel_;
  }

  const TransitionModelInterface& transitionModel() const override
  {
    return *originalModel_;
  }

protected:
  SubstitutionModelInterface& substitutionModel_() override
  {
    return *originalModel_;
  }


  TransitionModelInterface& transitionModel_() override
  {
    return *originalModel_;
  }

public:
  void fireParameterChanged(const ParameterList& parameters) override
  {
    substitutionModel_().matchParametersValues(parameters);
    AbstractParameterAliasable::fireParameterChanged(parameters);
    updateMatrices_();
  }

  size_t getNumberOfStates() const override
  {
    return stateMap().getNumberOfModelStates();
  }

public:
  std::string getName() const override
  {
    return "FromRegister";
  }

  const std::string& getRegisterName() const
  {
    return registerName_;
  }

  void addRateParameter() override
  {
    throw Exception("RegisterRatesSubstitutionModel::addRateParameter method should not be called, because rates are defined through registers.");
  }

  /**
   * @}
   */


   /**
    * @brief Overrides of AbstractSubstitutionModel and
    * AbstractWrappedSubstitutionModel.
    *
    * @{
    */
   const std::vector<int>& getAlphabetStates() const override
   {
     return AbstractWrappedSubstitutionModel::getAlphabetStates();
   }
  
   std::vector<size_t> getModelStates(int i) const override
   {
     return AbstractWrappedSubstitutionModel::getModelStates(i);
   }
  
   std::vector<size_t> getModelStates(const std::string& s) const override
   {
     return AbstractWrappedSubstitutionModel::getModelStates(s);
   }
  
   int getAlphabetStateAsInt(size_t i) const override
   {
     return AbstractWrappedSubstitutionModel::getAlphabetStateAsInt(i);
   }
  
   std::string getAlphabetStateAsChar(size_t s) const override
   {
     return AbstractWrappedSubstitutionModel::getAlphabetStateAsChar(s);
   }
  
   const Alphabet& alphabet() const override
   {
     return AbstractWrappedSubstitutionModel::alphabet();
   }
  
   std::shared_ptr<const Alphabet> getAlphabet() const override
   {
     return AbstractWrappedSubstitutionModel::getAlphabet();
   }
  
   const StateMapInterface& stateMap() const override
   {
     return AbstractWrappedSubstitutionModel::stateMap();
   }
  
   std::shared_ptr<const StateMapInterface> getStateMap() const override
   {
     return AbstractWrappedSubstitutionModel::getStateMap();
   }

   const FrequencySetInterface& frequencySet() const override
   {
     return AbstractWrappedSubstitutionModel::frequencySet();
   }

   /** @} */

  void setRate(double rate) override { model_().setRate(rate); }

  double getRate() const override { return model().getRate(); }

private:

  void setRegStates_();

protected:

  void updateMatrices_() override;

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_REGISTERRATESSUBSTITUTIONMODEL_H
