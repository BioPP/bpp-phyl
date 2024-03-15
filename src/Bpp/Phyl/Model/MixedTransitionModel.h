// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_MIXEDTRANSITIONMODEL_H
#define BPP_PHYL_MODEL_MIXEDTRANSITIONMODEL_H

#include <cstring> // C lib for string copy
#include <map>
#include <string>
#include <vector>

#include "SubstitutionModel.h"

namespace bpp
{
/**
 * @brief Interface for Transition models, defined as a mixture
 * of "simple" transition models.
 * @author Laurent Guéguen
 */
class MixedTransitionModelInterface :
  public virtual TransitionModelInterface
{
public:
  MixedTransitionModelInterface() {}

  virtual ~MixedTransitionModelInterface(){}

  virtual MixedTransitionModelInterface* clone() const override = 0;

public:
  /**
   * @brief Returns a specific model from the mixture
   */
  virtual const TransitionModelInterface& nModel(size_t i) const = 0;
  
  virtual std::shared_ptr<const TransitionModelInterface> getNModel(size_t i) const = 0;
  
  /**
   * @brief Returns the  probability of a specific model from the mixture
   */
  virtual double getNProbability(size_t i) const = 0;

  virtual const std::vector<double>& getProbabilities() const = 0;

  /**
   * @brief Sets the  probability of a specific model from the mixture
   */
  virtual void setNProbability(size_t i, double prob) = 0;

  virtual size_t getNumberOfModels() const = 0;

  /**
   * @brief Returns the rates of the submodels.
   */
  virtual const std::vector<double>& getVRates() const = 0;

  /**
   * @brief Returns the rate of a specific submodel.
   */
  virtual double getNRate(size_t i) const = 0;

  /**
   * @brief Sets the rates of the submodels to be proportional to a
   * given vector, and normalizes them so that the mean rate of the
   * mixture equals rate_.
   * @param vd a vector of positive values such that the rates of
   * the respective submodels are in the same proportions (ie this
   * vector does not need to be normalized).
   */
  virtual void setVRates(const Vdouble& vd) = 0;

  /**
   * @brief Normalizes the rates of the submodels so that the mean
   * rate of the mixture equals rate_.
   */
  virtual void normalizeVRates() = 0;

  /**
   * @brief Access the submodel with the given name.
   *
   * @throw NullPointerException if no model with the given name is found.
   */
  virtual const TransitionModelInterface& model(const std::string& name) const = 0;

  /**
   * @brief Returns the vector of numbers of the submodels in the
   * mixture that match a description.
   */
  virtual Vuint getSubmodelNumbers(const std::string& desc) const = 0;

protected:

  virtual TransitionModelInterface& nModel_(size_t i) = 0;

  friend class AbstractBiblioMixedTransitionModel;
  friend class InMixedSubstitutionModel;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_MIXEDTRANSITIONMODEL_H
