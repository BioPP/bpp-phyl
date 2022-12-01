//
// File: RegisterRatesSubstitutionModel.h
// Authors:
//   Laurent Gueguen
// Created: lundi 16 octobre 2017, ÃÂ  16h 38
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

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
  public AbstractTotallyWrappedSubstitutionModel
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
      std::unique_ptr<const SubstitutionModelInterface> originalModel,
      const SubstitutionRegisterInterface& reg,
      bool isNormalized = false);


  RegisterRatesSubstitutionModel(const RegisterRatesSubstitutionModel& fmsm) :
      AbstractWrappedModel(fmsm),
      AbstractWrappedTransitionModel(fmsm),
      AbstractTotallyWrappedTransitionModel(fmsm),
      AbstractWrappedSubstitutionModel(fmsm),
      AbstractTotallyWrappedSubstitutionModel(fmsm),
    originalModel_(fmsm.originalModel_->clone()),
    registerName_(fmsm.registerName_),
    vRegStates_(fmsm.vRegStates_),
    nbTypes_(fmsm.nbTypes_),
    vRates_(fmsm.vRates_)
  {}


  RegisterRatesSubstitutionModel& operator=(const RegisterRatesSubstitutionModel& fmsm)
  {
    AbstractTotallyWrappedSubstitutionModel::operator=(fmsm);
    
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
  SubstitutionModelInterface& substitutionModel() override
  {
    return *originalModel_;
  }


  TransitionModelInterface& transitionModel() override
  {
    return *originalModel_;
  }

public:
  void fireParameterChanged(const ParameterList& parameters) override
  {
    substitutionModel().matchParametersValues(parameters);

    AbstractParameterAliasable::fireParameterChanged(parameters);

    updateMatrices();
  }

  size_t getNumberOfStates() const override
  {
    return stateMap().getNumberOfModelStates();
  }

public:
  void updateMatrices() override;

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
  void setRate(double rate) override { model().setRate(rate); }

  double getRate() const override { return model().getRate(); }

private:
  void setRegStates_();
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_REGISTERRATESSUBSTITUTIONMODEL_H
