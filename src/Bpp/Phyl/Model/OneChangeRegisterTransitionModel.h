// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_ONECHANGEREGISTERTRANSITIONMODEL_H
#define BPP_PHYL_MODEL_ONECHANGEREGISTERTRANSITIONMODEL_H


#include "../Mapping/SubstitutionRegister.h"
#include "AbstractFromSubstitutionModelTransitionModel.h"
#include "AnonymousSubstitutionModel.h"

namespace bpp
{
/**
 * @brief From a model, compute transition probabilities given there
 * is at least a change of a category (ie a non null number in a
 * register) in the branch.
 *
 * It has the same parameters as the SubModel.
 *
 * @see SubstitutionRegister
 */
class OneChangeRegisterTransitionModel :
  public AbstractFromSubstitutionModelTransitionModel
{
private:
  /**
   * Boolean matrix of the sustitutions that are NOT considered (ie
   * for which the changes generator equal the ones of the original model).
   */
  RowMatrix<uint> noChangedStates_;

  /**
   * The SubstitutionModel in which generator has registered changes
   * set to 0.
   */
  std::unique_ptr<AnonymousSubstitutionModel> modelChanged_;

  /**
   * For output
   */
  std::string registerName_;

  /**
   * Vector of considered categories numbers in the register
   */
  std::vector<size_t> vNumRegs_;

public:
  /*
   * @brief Constructor
   *
   * @param originalModel the substitution model used
   * @param reg the register in which the considered type of event is
   * defined.
   * @param numReg the number of the considered event in the
   * register.
   */
  OneChangeRegisterTransitionModel(
      std::unique_ptr<SubstitutionModelInterface> originalModel,
      const SubstitutionRegisterInterface& reg, size_t numReg);

  /**
   * @brief Constructor
   *
   * @param originalModel the substitution model used
   * @param reg the register in which the considered type of event is
   * defined.
   * @param vNumRegs the vector of numbers of the considered event
   * in the register.
   */
  OneChangeRegisterTransitionModel(
      std::unique_ptr<SubstitutionModelInterface> originalModel,
      const SubstitutionRegisterInterface& reg,
      std::vector<size_t> vNumRegs);

  OneChangeRegisterTransitionModel(const OneChangeRegisterTransitionModel& fmsm) :
    AbstractParameterAliasable(fmsm),
    AbstractWrappedModel(fmsm),
    AbstractWrappedTransitionModel(fmsm),
    AbstractFromSubstitutionModelTransitionModel(fmsm),
    noChangedStates_(fmsm.noChangedStates_),
    modelChanged_(fmsm.modelChanged_->clone()),
    registerName_(fmsm.registerName_),
    vNumRegs_(fmsm.vNumRegs_)
  {}


  OneChangeRegisterTransitionModel& operator=(const OneChangeRegisterTransitionModel& fmsm)
  {
    AbstractWrappedModel::operator=(fmsm);
    AbstractWrappedTransitionModel::operator=(fmsm);
    AbstractFromSubstitutionModelTransitionModel::operator=(fmsm);
    noChangedStates_ = fmsm.noChangedStates_;
    modelChanged_ = std::unique_ptr<AnonymousSubstitutionModel>(fmsm.modelChanged_->clone());
    registerName_ = fmsm.registerName_;
    vNumRegs_ = fmsm.vNumRegs_;

    return *this;
  }

  virtual ~OneChangeRegisterTransitionModel() {}

  OneChangeRegisterTransitionModel* clone() const override {
    return new OneChangeRegisterTransitionModel(*this);
  }

public:
  void fireParameterChanged(const ParameterList& parameters) override
  {
    AbstractFromSubstitutionModelTransitionModel::fireParameterChanged(parameters);
    updateMatrices_();
  }

  double Pij_t    (size_t i, size_t j, double t) const override;
  double dPij_dt  (size_t i, size_t j, double t) const override;
  double d2Pij_dt2(size_t i, size_t j, double t) const override;

  const Matrix<double>& getPij_t(double t) const override;

  const Matrix<double>& getdPij_dt(double t) const override;

  const Matrix<double>& getd2Pij_dt2(double t) const override;

  double freq(size_t i) const override
  {
    return transitionModel().freq(i);
  }

  const Vdouble& getFrequencies() const override
  {
    return transitionModel().getFrequencies();
  }

  const FrequencySetInterface& frequencySet() const override
  {
    return transitionModel().frequencySet();
  }

  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount) override
  {
    transitionModel_().setFreqFromData(data, pseudoCount);
  }

  virtual void setFreq(std::map<int, double>& m) override
  {
    transitionModel_().setFreq(m);
  }

  double getRate() const override { return transitionModel().getRate(); }

  void setRate(double rate) override { return transitionModel_().setRate(rate); }

  double getInitValue(size_t i, int state) const override
  { 
    return model().getInitValue(i, state);
  }

  std::string getName() const override
  {
    return "OneChange";
  }

  const std::string& getRegisterName() const
  {
    return registerName_;
  }

  const std::vector<size_t>& getRegisterNumbers() const
  {
    return vNumRegs_;
  }

  /*
   * @}
   */

protected:  
  
  void updateMatrices_();

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ONECHANGEREGISTERTRANSITIONMODEL_H
