// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_ONECHANGETRANSITIONMODEL_H
#define BPP_PHYL_MODEL_ONECHANGETRANSITIONMODEL_H


#include "AbstractFromSubstitutionModelTransitionModel.h"

namespace bpp
{
/**
 * @brief From a model, compute transition probabilities given there
 * is at least a change in the branch.
 *
 * It has the same parameters as the SubModel.
 */

class OneChangeTransitionModel :
  public AbstractFromSubstitutionModelTransitionModel
{
public:
  OneChangeTransitionModel(std::unique_ptr<SubstitutionModelInterface> originalModel) :
    AbstractParameterAliasable("Onechange."),
    AbstractWrappedModel("Onechange."),
    AbstractWrappedTransitionModel("Onechange."),
    AbstractFromSubstitutionModelTransitionModel(std::move(originalModel), "OneChange.")
  {}

  OneChangeTransitionModel(const OneChangeTransitionModel& fmsm) :
    AbstractParameterAliasable(fmsm),
    AbstractWrappedModel(fmsm),
    AbstractWrappedTransitionModel(fmsm),
    AbstractFromSubstitutionModelTransitionModel(fmsm)
  {}

  OneChangeTransitionModel& operator=(const OneChangeTransitionModel& fmsm)
  {
    AbstractWrappedModel::operator=(fmsm);
    AbstractWrappedTransitionModel::operator=(fmsm);
    AbstractFromSubstitutionModelTransitionModel::operator=(fmsm);
    return *this;
  }

  virtual ~OneChangeTransitionModel() {}

  OneChangeTransitionModel* clone() const override { return new OneChangeTransitionModel(*this); }

public:
  double Pij_t    (size_t i, size_t j, double t) const override;
  double dPij_dt  (size_t i, size_t j, double t) const override;
  double d2Pij_dt2(size_t i, size_t j, double t) const override;

  const Matrix<double>& getPij_t(double t) const override;

  const Matrix<double>& getdPij_dt(double t) const override;

  const Matrix<double>& getd2Pij_dt2(double t) const override;

  double freq(size_t i) const override { return transitionModel().freq(i); }

  const Vdouble& getFrequencies() const override { return transitionModel().getFrequencies(); }

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

  /**
   * @}
   */
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ONECHANGETRANSITIONMODEL_H
