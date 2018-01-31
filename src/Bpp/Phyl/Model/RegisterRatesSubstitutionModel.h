//
// File: RegisterRatesSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: lundi 16 octobre 2017, à 16h 38
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _REGISTER_RATES_TRANSITION_MODEL_H_
#define _REGISTER_RATES_TRANSITION_MODEL_H_

#include <Bpp/Numeric/AbstractParameterAliasable.h>

#include "AbstractWrappedModel.h"
#include "AnonymousSubstitutionModel.h"
#include "AbstractSubstitutionModel.h"

#include "../Mapping/SubstitutionRegister.h"

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
    virtual public AbstractWrappedSubstitutionModel,
    virtual public AbstractSubstitutionModel
  {
  private:
    /*
     * @brief The related model.
     *
     */
    
    std::unique_ptr<SubstitutionModel> originalModel_;

    /*
     * For output
     *
     */

    std::string registerName_;

    /*
     * @brief Vector of register state -> vector of from states ->
     * vector of to states (for acceleration purpose)
     *
     */ 

    
    VVVuint vRegStates_;
    size_t nbTypes_;
    
    /*
     * @brief vector of the rates of the register types
     *
     */

    Vdouble vRates_;
    
  public:
    /*
     * @brief Constructor
     *
     * @param originalModel the substitution model used
     * @param reg the register in which the considered types of event are
     * used.
     * @param isNormalized says if model is normalized (default false)
     *
     */
    
    RegisterRatesSubstitutionModel(const SubstitutionModel& originalModel, const SubstitutionRegister& reg, bool isNormalized = false);
    

    RegisterRatesSubstitutionModel(const RegisterRatesSubstitutionModel& fmsm) :
      AbstractParameterAliasable(fmsm),
      AbstractSubstitutionModel(fmsm),
      originalModel_(fmsm.originalModel_->clone()),
      registerName_(fmsm.registerName_),
      vRegStates_(fmsm.vRegStates_),
      nbTypes_(fmsm.nbTypes_),
      vRates_(fmsm.vRates_)
    {}
    

    RegisterRatesSubstitutionModel& operator=(const RegisterRatesSubstitutionModel& fmsm)
    {
      AbstractSubstitutionModel::operator=(fmsm);
      originalModel_.reset(fmsm.originalModel_->clone());
      
      registerName_= fmsm.registerName_;
      vRegStates_ = fmsm.vRegStates_;
      nbTypes_ = fmsm.nbTypes_;
      vRates_ = fmsm.vRates_;
        
      return *this;
    }
    
    ~RegisterRatesSubstitutionModel() {}

    RegisterRatesSubstitutionModel* clone() const { return new RegisterRatesSubstitutionModel(*this); }

  public:

    /*
     * @brief clear overrides of AbstractSubstitutionModel and
     * AbstractWrappedSubstitutionModel.
     *
     */
    
    const std::vector<int>& getAlphabetStates() const
    {
      return AbstractWrappedSubstitutionModel::getAlphabetStates();
    }
    
    std::vector<long unsigned int> getModelStates(int i) const
    {
      return AbstractWrappedSubstitutionModel::getModelStates(i);
    }

    std::vector<long unsigned int> getModelStates(const std::string& s) const
    {
      return AbstractWrappedSubstitutionModel::getModelStates(s);
    }

    int getAlphabetStateAsInt(size_t i) const{
      return AbstractWrappedSubstitutionModel::getAlphabetStateAsInt(i);
    }

    std::string getAlphabetStateAsChar(size_t s) const
    {
      return AbstractWrappedSubstitutionModel::getAlphabetStateAsChar(s);
    }

    const Alphabet* getAlphabet() const
    {
      return AbstractWrappedSubstitutionModel::getAlphabet();
    }

    const StateMap& getStateMap() const
    {
      return AbstractWrappedSubstitutionModel::getStateMap();
    }

    
    /*
     * @brief From AbstractWrappedSubstitutionModel
     *
     */
      
    const SubstitutionModel& getSubstitutionModel() const
    {
      return *originalModel_.get();
    }

    const TransitionModel& getModel() const
    {
      return *originalModel_.get();
    }

  protected:
    SubstitutionModel& getSubstitutionModel()
    {
      return *originalModel_.get();
    }

    
    TransitionModel& getModel()
    {
      return *originalModel_.get();
    }

  public:

    /*
     * @brief From AbstractSubstitutionModel
     *
     */
    
    void fireParameterChanged(const ParameterList& parameters)
    {
      getSubstitutionModel().matchParametersValues(parameters);

      AbstractParameterAliasable::fireParameterChanged(parameters);
      
      updateMatrices();
    }
    
    size_t getNumberOfStates() const
    {
      return getModel().getNumberOfStates();
    }

  public:

    void updateMatrices();
    
    std::string getName() const
    {
      return "FromRegister";
    }

    const std::string& getRegisterName() const
    {
      return registerName_;
    }

    void addRateParameter()
    {
      throw Exception("RegisterRatesSubstitutionModel::addRateParameter method should not be called, because rates are defined through registers.");
    }

    /*
     * @}
     *
     */

    void setRate(double rate) { getModel().setRate(rate); }

    double getRate() const { return getModel().getRate(); }


  private:

    void setRegStates_();
      
  };
} // end of namespace bpp.

#endif  // _REGISTER_RATES_TRANSITION_MODEL_H_
