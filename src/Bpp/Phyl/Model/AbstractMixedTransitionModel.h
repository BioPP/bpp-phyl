// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_ABSTRACTMIXEDTRANSITIONMODEL_H
#define BPP_PHYL_MODEL_ABSTRACTMIXEDTRANSITIONMODEL_H

#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <cstring> // C lib for string copy
#include <map>
#include <string>
#include <vector>

#include "AbstractSubstitutionModel.h"
#include "MixedTransitionModel.h"

namespace bpp
{
/**
 * @brief Partial implementation for Mixed Transition models,
 * defined as a mixture of "simple" substitution models. Each model
 * has a specific probability and rate, with the constraint that the
 * expectation (on the distribution of the models) of the rate of
 * all the models equals one.
 *
 * In this kind of model, there is no generator.
 *
 * @author Laurent Guéguen
 */
class AbstractMixedTransitionModel :
  public virtual MixedTransitionModelInterface,
  public virtual AbstractTransitionModel
{
protected:
  /**
   * @brief vector of pointers to TransitionModels.
   *
   * Beware: these TransitionModels are owned by the object, so
   * will be deleted at destruction
   */
  std::vector< std::shared_ptr<TransitionModelInterface>> modelsContainer_;

  /**
   * @brief vector of the probabilities of the models
   */
  std::vector<double> vProbas_;

  /**
   * @brief vector of the rates of the models.
   *
   * For the computation of the transition probabilities, the rates
   * are included in the submodels while updating the mixture, so
   * there is no need to multiply here the transition times with the
   * rates.
   *
   * The mean (on the distribution of the models) of the elements of
   * this vector equals the overall rate of the mixture model, that
   * is rate_;
   */
  std::vector<double> vRates_;

public:
  AbstractMixedTransitionModel(
      std::shared_ptr<const Alphabet>,
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& prefix);

  AbstractMixedTransitionModel(const AbstractMixedTransitionModel&);

  AbstractMixedTransitionModel& operator=(const AbstractMixedTransitionModel&);

  virtual AbstractMixedTransitionModel* clone() const override = 0;

public:
  /**
   * @brief returns the number of models in the mixture
   */
  virtual size_t getNumberOfModels() const override
  {
    return modelsContainer_.size();
  }

  /**
   * @brief Returns a specific model from the mixture
   */
  const TransitionModelInterface& nModel(size_t i) const override
  {
    return *modelsContainer_[i];
  }

  std::shared_ptr<const TransitionModelInterface> getNModel(size_t i) const override
  {
    return modelsContainer_[i];
  }

  /**
   * @brief Returns the rate of a specific model from the mixture
   */
  double getNRate(size_t i) const override
  {
    return vRates_[i];
  }

  /**
   * @brief Set the rate of the model and the submodels.
   * @param rate must be positive.
   */
  virtual void setRate(double rate) override;

  /**
   * @brief Sets the rates of the submodels to be proportional to a
   * given vector, with the constraint that the mean rate of the
   * mixture equals rate_.
   * @param vd a vector of positive values such that the rates of
   * the respective submodels are in the same proportions (ie this
   * vector does not need to be normalized).
   */
  virtual void setVRates(const Vdouble& vd) override;

  /**
   * @brief Normalizes the rates of the submodels so that the mean
   * rate of the mixture equals rate_.
   */
  virtual void normalizeVRates() override;

  /**
   * @brief Returns the vector of all the rates of the mixture
   */
  const std::vector<double>& getVRates() const override
  {
    return vRates_;
  }

  /**
   * @brief Returns the probability of a specific model from the
   * mixture
   */
  virtual double getNProbability(size_t i) const override
  {
    return vProbas_[i];
  }

  /**
   * @brief Returns the vector of probabilities
   *
   */
  virtual const std::vector<double>& getProbabilities() const override
  {
    return vProbas_;
  }

  /**
   * @brief Sets the  probability of a specific model from the mixture
   */
  virtual void setNProbability(size_t i, double prob) override
  {
    if (prob < 0)
      prob = 0;
    if (prob > 1)
      prob = 1;

    vProbas_[i] = prob;
  }

  /**
   * @brief From TransitionModel interface
   */
  virtual const Matrix<double>& getPij_t(double t) const override;
  virtual const Matrix<double>& getdPij_dt(double t) const override;
  virtual const Matrix<double>& getd2Pij_dt2(double t) const override;

  /**
   * @return Says if equilibrium frequencies should be computed (all
   * models are likewise, may be refined)
   */
  bool computeFrequencies() const override
  {
    return modelsContainer_[0]->computeFrequencies();
  }

  /**
   * @return Set if equilibrium frequencies should be computed
   */
  void computeFrequencies(bool yn) override
  {
    for (auto& sm : modelsContainer_)
    {
      sm->computeFrequencies(yn);
    }
  }

  void setFreq(std::map<int, double>& frequ) override
  {
    for (auto& sm : modelsContainer_)
    {
      sm->setFreq(frequ);
    }
  }

  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount) override
  {
    std::map<int, double> freqs;
    SequenceContainerTools::getFrequencies(data, freqs, pseudoCount);
    setFreq(freqs);
  }

protected:
  TransitionModelInterface& nModel_(size_t i) override
  {
    return *modelsContainer_[i];
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ABSTRACTMIXEDTRANSITIONMODEL_H
