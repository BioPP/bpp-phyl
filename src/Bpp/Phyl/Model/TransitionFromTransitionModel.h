//
// File: TransitionFromTransitionModel.h
// Created by: Laurent Gueguen
// Created on: dimanche 26 janvier 2020, à 07h 52
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

#ifndef _MULTINOMIAL_FROM_TRANSITION_MODEL_H_
#define _MULTINOMIAL_FROM_TRANSITION_MODEL_H_

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
  virtual public AbstractParameterAliasable,
  virtual public AbstractWrappedModel
{
private:
  /*
   * @brief The related model.
   *
   */

  std::unique_ptr<TransitionModel> subModel_;

  /*
   * The number of states
   *
   */

  size_t size_;

  /*
   * @brief Reference time to avoid recomuputation of transition
   * matrix when time has not changed. If <0, it means that
   * transition matrix should be recomputed (for ex if parameters
   * have changed).
   *
   */

  mutable double tref_;

  /*
   * @brief Transition Matrices owned by the submodel.
   *
   */

  /**
   * @brief These ones are for bookkeeping:
   */
  mutable const Matrix<double>* Pij_t, * dPij_dt, * d2Pij_dt2;

  /*
   * @brief Used return vectors
   *
   */

  mutable Eigen::VectorXd Pi_, dPi_, d2Pi_;

protected:
  BranchModel& getModel()
  {
    return *subModel_.get();
  }

  TransitionModel& getTransitionModel()
  {
    return *subModel_.get();
  }

public:
  TransitionFromTransitionModel(const TransitionModel& originalModel) :
    AbstractParameterAliasable("TransitionFrom." + originalModel.getNamespace()),
    AbstractWrappedModel(),
    subModel_(std::unique_ptr<TransitionModel>(originalModel.clone())),
    size_(originalModel.getNumberOfStates()),
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
    AbstractParameterAliasable::operator=(fmsm);

    subModel_ = std::unique_ptr<TransitionModel>(fmsm.subModel_->clone());
    size_ = fmsm.size_;
    Pi_.resize(Eigen::Index(size_));
    dPi_.resize(Eigen::Index(size_));
    d2Pi_.resize(Eigen::Index(size_));

    tref_ = -1;

    return *this;
  }

  ~TransitionFromTransitionModel() {}

  TransitionFromTransitionModel* clone() const { return new TransitionFromTransitionModel(*this); }

public:
  void fireParameterChanged(const ParameterList& parameters)
  {
    AbstractParameterAliasable::fireParameterChanged(parameters);
    if (getModel().matchParametersValues(parameters))
      tref_ = -1;
  }

  const BranchModel& getModel() const
  {
    return *subModel_.get();
  }

  const TransitionModel& getTransitionModel() const
  {
    return *subModel_.get();
  }

  const Eigen::VectorXd& Lik_t    (const Eigen::VectorXd& from, double t) const;
  const Eigen::VectorXd& dLik_dt  (const Eigen::VectorXd& from, double t) const;
  const Eigen::VectorXd& d2Lik_dt2(const Eigen::VectorXd& from, double t) const;

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

  double getInitValue(size_t i, int state) const { return getTransitionModel().getInitValue(i, state); }

  std::string getName() const
  {
    return "TransitionFrom";
  }

  void addRateParameter()
  {
    getModel().addRateParameter();
    addParameter_(new Parameter(getNamespace() + "rate", getModel().getRate(), Parameter::R_PLUS_STAR));
  }
};
} // end of namespace bpp.

#endif// _MULTINOMIAL_FROM_TRANSITION_MODEL_H_
