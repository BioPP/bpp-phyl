// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_TRANSITIONFROMTRANSITIONMODEL_H
#define BPP_PHYL_MODEL_TRANSITIONFROMTRANSITIONMODEL_H


#include "AbstractWrappedModel.h"

namespace bpp
{
/**
 * @brief From a transition model, compute the transition function
 * probabilities.
 *
 * This class is (up to now) mostly for test purpose.
 *
 * It has the same parameters as the SubModel.
 */
class TransitionFromTransitionModel :
  public AbstractWrappedModel
{
private:
  /**
   * @brief The related model.
   */
  std::shared_ptr<TransitionModelInterface> subModel_;

  /**
   * The number of states
   */
  size_t size_;

  /**
   * @brief Reference time to avoid recomuputation of transition
   * matrix when time has not changed. If <0, it means that
   * transition matrix should be recomputed (for ex if parameters
   * have changed).
   */
  mutable double tref_;

  /**
   * @brief Transition Matrices owned by the submodel.
   */

  /**
   * @brief These ones are for bookkeeping:
   */
  mutable const Matrix<double>* Pij_t, * dPij_dt, * d2Pij_dt2;

  /**
   * @brief Used return vectors
   */
  mutable Eigen::VectorXd Pi_, dPi_, d2Pi_;

public:
  TransitionFromTransitionModel(std::shared_ptr<TransitionModelInterface> originalModel) :
    AbstractParameterAliasable("TransitionFrom." + originalModel->getNamespace()),
    AbstractWrappedModel("TransitionFrom." + originalModel->getNamespace()),
    subModel_(originalModel),
    size_(originalModel->getNumberOfStates()),
    tref_(-1), Pij_t(0), dPij_dt(0), d2Pij_dt2(0), Pi_(size_), dPi_(size_), d2Pi_(size_)
  {
    subModel_->setNamespace(getNamespace());
    addParameters_(subModel_->getParameters());
  }

  TransitionFromTransitionModel(const TransitionFromTransitionModel& fmsm) :
    AbstractParameterAliasable(fmsm),
    AbstractWrappedModel(fmsm),
    subModel_(fmsm.subModel_->clone()),
    size_(fmsm.size_),
    tref_(-1), Pij_t(0), dPij_dt(0), d2Pij_dt2(0), Pi_(size_), dPi_(size_), d2Pi_(size_)
  {}

  TransitionFromTransitionModel& operator=(const TransitionFromTransitionModel& fmsm)
  {
    AbstractWrappedModel::operator=(fmsm);

    subModel_ = std::shared_ptr<TransitionModelInterface>(fmsm.subModel_->clone());
    size_ = fmsm.size_;
    Pi_.resize(Eigen::Index(size_));
    dPi_.resize(Eigen::Index(size_));
    d2Pi_.resize(Eigen::Index(size_));

    tref_ = -1;

    return *this;
  }

  virtual ~TransitionFromTransitionModel() {}

  TransitionFromTransitionModel* clone() const override { return new TransitionFromTransitionModel(*this); }

public:
  void fireParameterChanged(const ParameterList& parameters) override
  {
    AbstractParameterAliasable::fireParameterChanged(parameters);
    if (model_().matchParametersValues(parameters))
      tref_ = -1;
  }

  const BranchModelInterface& model() const override
  {
    return *subModel_;
  }

  const TransitionModelInterface& transitionModel() const
  {
    return *subModel_;
  }

  const Eigen::VectorXd& Lik_t    (const Eigen::VectorXd& from, double t) const override;
  const Eigen::VectorXd& dLik_dt  (const Eigen::VectorXd& from, double t) const override;
  const Eigen::VectorXd& d2Lik_dt2(const Eigen::VectorXd& from, double t) const override;

  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount)
  {
    transitionModel_().setFreqFromData(data, pseudoCount);
  }

  virtual void setFreq(std::map<int, double>& m)
  {
    transitionModel_().setFreq(m);
  }

  double getRate() const override { return transitionModel().getRate(); }

  void setRate(double rate) override { return transitionModel_().setRate(rate); }

  double getInitValue(size_t i, int state) const override
  {
    return transitionModel().getInitValue(i, state);
  }

  std::string getName() const override
  {
    return "TransitionFrom";
  }

  void addRateParameter() override
  {
    model_().addRateParameter();
    addParameter_(new Parameter(getNamespace() + "rate", model().getRate(), Parameter::R_PLUS_STAR));
  }

protected:
  BranchModelInterface& model_()
  {
    return *subModel_;
  }

  TransitionModelInterface& transitionModel_()
  {
    return *subModel_;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_TRANSITIONFROMTRANSITIONMODEL_H
