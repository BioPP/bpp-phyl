//
// File: AbstractBiblioMixedTransitionModel.h
// Authors:
//   Laurent Guéguen
// Created: lundi 18 juillet 2011, ÃÂ  15h 17
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_MODEL_ABSTRACTBIBLIOMIXEDTRANSITIONMODEL_H
#define BPP_PHYL_MODEL_ABSTRACTBIBLIOMIXEDTRANSITIONMODEL_H


#include "AbstractBiblioSubstitutionModel.h"
#include "MixedTransitionModel.h"

namespace bpp
{
/**
 * @brief Abstract class for mixture models based on the bibliography.
 * @author Laurent Guéguen
 */
class AbstractBiblioMixedTransitionModel :
  public virtual MixedTransitionModelInterface,
  public virtual AbstractBiblioTransitionModel
{
protected:
  std::unique_ptr<MixedTransitionModelInterface> mixedModelPtr_;

public:
  AbstractBiblioMixedTransitionModel(const std::string& prefix);

  AbstractBiblioMixedTransitionModel(const AbstractBiblioMixedTransitionModel& model);

  AbstractBiblioMixedTransitionModel& operator=(const AbstractBiblioMixedTransitionModel& model);

  virtual ~AbstractBiblioMixedTransitionModel();

public:
  /**
   * @brief Returns the submodel from the mixture.
   */
  const TransitionModelInterface& nModel(size_t i) const override
  {
    return mixedModel().nModel(i);
  }

  std::shared_ptr<const TransitionModelInterface> getNModel(size_t i) const override
  {
    return mixedModel().getNModel(i);
  }

  /**
   * @brief Returns the  probability of a specific model from the mixture
   */
  double getNProbability(size_t i) const override
  {
    return mixedModel().getNProbability(i);
  }

  /**
   * @brief Returns the vector of the probabilities of the
   * submodels of the mixture.
   *
   */
  const std::vector<double>& getProbabilities() const override
  {
    return mixedModel().getProbabilities();
  }

  /**
   * @brief Sets the probabilities of the submodels of the mixture.
   *
   */
  void setNProbability(size_t i, double prob) override
  {
    mixedModel_().setNProbability(i, prob);
  }

  /**
   * @brief Returns the number of submodels
   */
  size_t getNumberOfModels() const override
  {
    return mixedModel().getNumberOfModels();
  }

  /**
   * @brief sets the rates of the submodels.
   *
   **/
  void setVRates(const Vdouble& vd) override
  {
    mixedModel_().setVRates(vd);
  }

  /**
   * @brief normalizes the rates of the submodels.
   *
   **/
  void normalizeVRates() override
  {
    mixedModel_().normalizeVRates();
  }

  /**
   * @brief Returns the vector of all the rates of the mixture
   */
  const std::vector<double>& getVRates() const override
  {
    return mixedModel().getVRates();
  }

  /**
   * @brief Returns the rate of a specific model from the mixture
   */
  double getNRate(size_t i) const override
  {
    return mixedModel().getNRate(i);
  }

  /**
   * @brief retrieve a pointer to the submodel with the given name.
   *
   * Return Null if not found.
   */
  using AbstractWrappedModel::model;
  const TransitionModelInterface& model(const std::string& name) const override
  {
    return mixedModel().model(name);
  }

  /**
   * @brief Returns the vector of numbers of the submodels in the
   * mixture that match a description.
   */
  Vuint getSubmodelNumbers(const std::string& desc) const override;

  const TransitionModelInterface& transitionModel() const override { return *mixedModelPtr_; }

  const MixedTransitionModelInterface& mixedModel() const { return *mixedModelPtr_; }

protected:
  
  TransitionModelInterface& transitionModel_() override
  {
    return *mixedModelPtr_;
  }

  MixedTransitionModelInterface& mixedModel_() { return *mixedModelPtr_; }

  const FrequencySetInterface& frequencySet_() const
  {
    return mixedModelPtr_->nModel(0).frequencySet();
  }

  TransitionModelInterface& nModel_(size_t i) override
  {
    return mixedModel_().nModel_(i);
  }

  
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ABSTRACTBIBLIOMIXEDTRANSITIONMODEL_H
