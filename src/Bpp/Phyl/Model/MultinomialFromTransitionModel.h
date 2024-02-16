//
// File: MultinomialFromTransitionModel.h
// Authors:
//   Laurent Gueguen
// Created: dimanche 26 janvier 2020, Ã  07h 52
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_MODEL_MULTINOMIALFROMTRANSITIONMODEL_H
#define BPP_PHYL_MODEL_MULTINOMIALFROMTRANSITIONMODEL_H

#include <unsupported/Eigen/SpecialFunctions>

#include "AbstractWrappedModel.h"

// namespace std
// {
// //  template<class U = Eigen::VectorXd>
//   template<>


// }

namespace bpp
{
/**
 * @brief From a model, compute the likelihood of counts given an
 * ancestral state.
 *
 * It has the same parameters as the SubModel.
 *
 * If counts are not integers, rounded values are used.
 */
class MultinomialFromTransitionModel :
  public AbstractWrappedModel
{
  typedef bool lessEigenType(Eigen::VectorXd const&, Eigen::VectorXd const&);

private:

  /**
   * @brief The related model.
   */
  std::shared_ptr<TransitionModelInterface> subModel_;  // -> shared? 

  /*
   * The number of states
   *
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

  /**
   * @brief map to store constant values of multinomial
   */
  mutable std::map<Eigen::VectorXd, double, bool (*)(const Eigen::VectorXd&, const Eigen::VectorXd&)> mapFact_;

  static bool lessEigen(Eigen::VectorXd const& a, Eigen::VectorXd const& b)
  {
    assert(a.size() == b.size());
    for (auto i = 0; i < a.size(); ++i)
    {
      if (a[i] < b[i]) return true;
      if (a[i] > b[i]) return false;
    }
    return false;
  }

protected:
  BranchModelInterface& model()
  {
    return *subModel_;
  }

  TransitionModelInterface& transitionModel()
  {
    return *subModel_;
  }

public:
  MultinomialFromTransitionModel(std::shared_ptr<TransitionModelInterface> originalModel) :
    AbstractParameterAliasable("MultinomialFrom." + originalModel->getNamespace()),
    AbstractWrappedModel("MultinomialFrom." + originalModel->getNamespace()),
    subModel_(originalModel),
    size_(originalModel->getNumberOfStates()),
    tref_(NumConstants::MINF()), Pij_t(0), dPij_dt(0), d2Pij_dt2(0), Pi_(size_), dPi_(size_), d2Pi_(size_), mapFact_(&lessEigen)
  {
    subModel_->setNamespace(getNamespace());
    addParameters_(subModel_->getParameters());
  }

  MultinomialFromTransitionModel(const MultinomialFromTransitionModel& fmsm) :
    AbstractParameterAliasable(fmsm),
    AbstractWrappedModel(fmsm),
    subModel_(fmsm.subModel_->clone()),
    size_(fmsm.size_),
    tref_(NumConstants::MINF()), Pij_t(0), dPij_dt(0), d2Pij_dt2(0), Pi_(size_), dPi_(size_), d2Pi_(size_), mapFact_(&lessEigen)
  {}

  MultinomialFromTransitionModel& operator=(const MultinomialFromTransitionModel& fmsm)
  {
    AbstractWrappedModel::operator=(fmsm);

    subModel_ = std::shared_ptr<TransitionModelInterface>(fmsm.subModel_->clone());
    size_ = fmsm.size_;
    Pi_.resize(Eigen::Index(size_));
    dPi_.resize(Eigen::Index(size_));
    d2Pi_.resize(Eigen::Index(size_));

    tref_ = NumConstants::MINF();
    mapFact_.clear();

    return *this;
  }

  virtual ~MultinomialFromTransitionModel() {}

  MultinomialFromTransitionModel* clone() const override { return new MultinomialFromTransitionModel(*this); }

public:
  void fireParameterChanged(const ParameterList& parameters) override
  {
    AbstractParameterAliasable::fireParameterChanged(parameters);
    if (model().matchParametersValues(parameters))
      tref_ = NumConstants::MINF();
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

  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0)
  {
    transitionModel().setFreqFromData(data, pseudoCount);
  }

  virtual void setFreq(std::map<int, double>& m)
  {
    transitionModel().setFreq(m);
  }

  double getRate() const override { return transitionModel().getRate(); }

  void setRate(double rate) override { return transitionModel().setRate(rate); }

  double getInitValue(size_t i, int state) const override { return transitionModel().getInitValue(i, state); }

  std::string getName() const override
  {
    return "MultinomialFrom";
  }

  void addRateParameter() override
  {
    model().addRateParameter();
    addParameter_(new Parameter(getNamespace() + "rate", model().getRate(), Parameter::R_PLUS_STAR));
  }

private:
  /**
   * @brief Fills the res vector of the likelihoods of the counts
   * given ancestral states & transition matrix.
   */
  void compute_Multinomial_(const Eigen::VectorXd& counts) const;
  void compute_dMultinomial_dt_(const Eigen::VectorXd& counts) const;
  void compute_d2Multinomial_dt2_(const Eigen::VectorXd& counts) const;

  static uint factorial(uint n)
  {
    return (n == 0 || n == 1) ? 1 : factorial(n - 1) * n;
  }

  /**
   * @brief Compute the log of the normalization term for the
   * multinomial for any count, and keeps it in a map to avoid
   * recalcultation.
   *
   */
  
  double getFact_(const Eigen::VectorXd& counts) const
  {
    auto it = mapFact_.find(counts);
    if (it == mapFact_.end())
    {
      auto lsto(std::lgamma(counts.sum() + 1));
      auto lr((counts.array() + 1.).lgamma().sum());
      auto fact = lsto - lr;

      mapFact_[counts] = fact;
      return fact;
    }
    else
      return it->second;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_MULTINOMIALFROMTRANSITIONMODEL_H
