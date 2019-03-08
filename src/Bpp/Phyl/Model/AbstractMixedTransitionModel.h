//
// File: AbstractMixedTransitionModel.h
// Created by: Laurent Gueguen
// On: vendredi 19 novembre 2010, à 15h 48
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

#ifndef _ABSTRACT_MIXED_TRANSITION_MODEL_H_
#define _ABSTRACT_MIXED_TRANSITION_MODEL_H_

#include "MixedTransitionModel.h"
#include "AbstractSubstitutionModel.h"

#include <vector>
#include <string>
#include <map>
#include <cstring> // C lib for string copy

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
 *
 */

  class AbstractMixedTransitionModel :
    public virtual MixedTransitionModel,
    public virtual AbstractTransitionModel
  {
  protected:
    /**
     * @brief vector of pointers to TransitionModels.
     *
     * Beware: these TransitionModels are owned by the object, so
     * will be deleted at destruction
     */
    std::vector<TransitionModel*> modelsContainer_;

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
    AbstractMixedTransitionModel(const Alphabet*, std::shared_ptr<const StateMap> stateMap, const std::string& prefix);

    AbstractMixedTransitionModel(const AbstractMixedTransitionModel&);

    AbstractMixedTransitionModel& operator=(const AbstractMixedTransitionModel&);

    virtual ~AbstractMixedTransitionModel();

    virtual AbstractMixedTransitionModel* clone() const = 0;

  public:
    /**
     * @brief returns the number of models in the mixture
     */
    virtual size_t getNumberOfModels() const
    {
      return modelsContainer_.size();
    }

    /**
     * @brief Returns a specific model from the mixture
     */
    const TransitionModel* getNModel(size_t i) const
    {
      return modelsContainer_[i];
    }

    TransitionModel* getNModel(size_t i)
    {
      return modelsContainer_[i];
    }

    /**
     * @brief Returns the rate of a specific model from the mixture
     */
    double getNRate(size_t i) const
    {
      return vRates_[i];
    }

    /**
     * @brief Set the rate of the model and the submodels.
     * @param rate must be positive.
     */

    virtual void setRate(double rate);

    /**
     * @brief Sets the rates of the submodels to be proportional to a
     * given vector, with the constraint that the mean rate of the
     * mixture equals rate_.

     * @param vd a vector of positive values such that the rates of
     * the respective submodels are in the same proportions (ie this
     * vector does not need to be normalized).
     */

    virtual void setVRates(const Vdouble& vd);

    /**
     * @brief Normalizes the rates of the submodels so that the mean
     * rate of the mixture equals rate_.
     */

    virtual void normalizeVRates();

    /**
     * @brief Returns the vector of all the rates of the mixture
     */

    const std::vector<double>& getVRates() const
    {
      return vRates_;
    }

    /**
     * @brief Returns the probability of a specific model from the
     * mixture
     */
    virtual double getNProbability(size_t i) const
    {
      return vProbas_[i];
    }

    /**
     * @brief Returns the vector of probabilities
     *
     */

    virtual const std::vector<double>& getProbabilities() const
    {
      return vProbas_;
    }

    /**
     * @brief Sets the  probability of a specific model from the mixture
     */
    virtual void setNProbability(size_t i, double prob)
    {
      if ((prob >= 0) && (prob <= 1))
        vProbas_[i] = prob;
    }

    /**
     * @brief From TransitionModel interface
     *
     */

    virtual size_t getNumberOfStates() const;

    virtual const Matrix<double>& getPij_t(double t) const;
    virtual const Matrix<double>& getdPij_dt(double t) const;
    virtual const Matrix<double>& getd2Pij_dt2(double t) const;

    /**
     * @return Says if equilibrium frequencies should be computed (all
     * models are likewise, may be refined)
     */
    
    bool computeFrequencies() const
    {
      return modelsContainer_[0]->computeFrequencies();
    }

    /**
     * @return Set if equilibrium frequencies should be computed
     */
    
    void computeFrequencies(bool yn)
    {
      for (auto& sm : modelsContainer_)
        sm->computeFrequencies(yn);
    }

  };
} // end of namespace bpp.

#endif  // _ABSTRACT_MIXED_TRANSITION_MODEL_H_
