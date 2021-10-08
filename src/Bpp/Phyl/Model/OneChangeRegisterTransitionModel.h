//
// File: OneChangeRegisterTransitionModel.h
// Created by: Laurent Gueguen
// Created on: samedi 24 octobre 2015, à 18h 28
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

#ifndef _ONE_CHANGE_REGISTER_TRANSITION_MODEL_H_
#define _ONE_CHANGE_REGISTER_TRANSITION_MODEL_H_

#include "AbstractFromSubstitutionModelTransitionModel.h"
#include "AnonymousSubstitutionModel.h"

#include "../Mapping/SubstitutionRegister.h"

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
  /*
   * Boolean matrix of the sustitutions that are NOT considered (ie
   * for which the changes generator equal the ones of the original model).
   *
   */

  RowMatrix<uint> noChangedStates_;

  /*
   * The SubstitutionModel in which generator has registered changes
   * set to 0.
   */

  std::unique_ptr<AnonymousSubstitutionModel> modelChanged_;

  /*
   * For output
   *
   */

  std::string registerName_;

  /*
   *
   * Vector of considered categories numbers in the register
   *
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
   *
   */

  OneChangeRegisterTransitionModel(const SubstitutionModel& originalModel, const SubstitutionRegister& reg, size_t numReg);

  /*
   * @brief Constructor
   *
   * @param originalModel the substitution model used
   * @param reg the register in which the considered type of event is
   * defined.
   * @param vNumRegs the vector of numbers of the considered event
   * in the register.
   *
   */

  OneChangeRegisterTransitionModel(const SubstitutionModel& originalModel, const SubstitutionRegister& reg, std::vector<size_t> vNumRegs);

  OneChangeRegisterTransitionModel(const OneChangeRegisterTransitionModel& fmsm) :
    AbstractParameterAliasable(fmsm),
    AbstractFromSubstitutionModelTransitionModel(fmsm),
    noChangedStates_(fmsm.noChangedStates_),
    modelChanged_(std::unique_ptr<AnonymousSubstitutionModel>(fmsm.modelChanged_->clone())),
    registerName_(fmsm.registerName_),
    vNumRegs_(fmsm.vNumRegs_)
  {}


  OneChangeRegisterTransitionModel& operator=(const OneChangeRegisterTransitionModel& fmsm)
  {
    AbstractFromSubstitutionModelTransitionModel::operator=(fmsm);
    noChangedStates_ = fmsm.noChangedStates_;
    modelChanged_ = std::unique_ptr<AnonymousSubstitutionModel>(fmsm.modelChanged_->clone());
    registerName_ = fmsm.registerName_;
    vNumRegs_ = fmsm.vNumRegs_;

    return *this;
  }

  ~OneChangeRegisterTransitionModel()
  {}

  OneChangeRegisterTransitionModel* clone() const { return new OneChangeRegisterTransitionModel(*this); }

public:
  void fireParameterChanged(const ParameterList& parameters)
  {
    AbstractFromSubstitutionModelTransitionModel::fireParameterChanged(parameters);
    updateMatrices();
  }

  double Pij_t    (size_t i, size_t j, double t) const;
  double dPij_dt  (size_t i, size_t j, double t) const;
  double d2Pij_dt2(size_t i, size_t j, double t) const;

  const Matrix<double>& getPij_t(double t) const;

  const Matrix<double>& getdPij_dt(double t) const;

  const Matrix<double>& getd2Pij_dt2(double t) const;


  double freq(size_t i) const { return getTransitionModel().freq(i); }

  const Vdouble& getFrequencies() const { return getTransitionModel().getFrequencies(); }

  const std::shared_ptr<FrequencySet> getFrequencySet() const {return getTransitionModel().getFrequencySet(); }

  void setFreqFromData(const SequencedValuesContainer& data, double pseudoCount)
  {
    getTransitionModel().setFreqFromData(data, pseudoCount);
  }

  virtual void setFreq(std::map<int, double>& m)
  {
    getTransitionModel().setFreq(m);
  }

  double getRate() const { return getTransitionModel().getRate(); }

  void setRate(double rate) { return getTransitionModel().setRate(rate); }

  double getInitValue(size_t i, int state) const { return getModel().getInitValue(i, state); }

  void updateMatrices();

  std::string getName() const
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
   *
   */
};
} // end of namespace bpp.

#endif// _ONE_CHANGE_REGISTER_TRANSITION_MODEL_H_
