//
// File: AbstractBiblioMixedTransitionModel.h
// Authors:
//   Laurent Gueguen
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
  std::shared_ptr<MixedTransitionModelInterface> pmixmodel_;

public:
  AbstractBiblioMixedTransitionModel(const std::string& prefix);

  AbstractBiblioMixedTransitionModel(const AbstractBiblioMixedTransitionModel& model);

  AbstractBiblioMixedTransitionModel& operator=(const AbstractBiblioMixedTransitionModel& model);

  virtual ~AbstractBiblioMixedTransitionModel();

public:
  /**
   * @brief Returns the submodel from the mixture.
   */
  const TransitionModelInterface& nModel(size_t i) const
  {
    return mixedModel().nModel(i);
  }

  std::shared_ptr<const TransitionModelInterface> getNModel(size_t i) const
  {
    return mixedModel().getNModel(i);
  }

  TransitionModelInterface& nModel(size_t i)
  {
    return mixedModel().nModel(i);
  }

  std::shared_ptr<TransitionModelInterface> getNModel(size_t i)
  {
    return mixedModel().getNModel(i);
  }

  /**
   * @brief Returns the  probability of a specific model from the mixture
   */
  double getNProbability(size_t i) const
  {
    return mixedModel().getNProbability(i);
  }

  /**
   * @brief Returns the vector of the probabilities of the
   * submodels of the mixture.
   *
   */
  const std::vector<double>& getProbabilities() const
  {
    return mixedModel().getProbabilities();
  }

  /**
   * @brief Sets the probabilities of the submodels of the mixture.
   *
   */
  void setNProbability(size_t i, double prob)
  {
    mixedModel().setNProbability(i, prob);
  }

  /**
   * @brief Returns the number of submodels
   *
   */
  size_t getNumberOfModels() const
  {
    return mixedModel().getNumberOfModels();
  }

  /**
   * @brief sets the rates of the submodels.
   *
   **/
  void setVRates(const Vdouble& vd)
  {
    mixedModel().setVRates(vd);
  }

  /**
   * @brief normalizes the rates of the submodels.
   *
   **/
  void normalizeVRates()
  {
    mixedModel().normalizeVRates();
  }

  /**
   * @brief Returns the vector of all the rates of the mixture
   */
  const std::vector<double>& getVRates() const
  {
    return mixedModel().getVRates();
  }

  /**
   * @brief Returns the rate of a specific model from the mixture
   */
  double getNRate(size_t i) const
  {
    return mixedModel().getNRate(i);
  }

  /**
   * @brief retrieve a pointer to the submodel with the given name.
   *
   * Return Null if not found.
   */
  using AbstractWrappedModel::model;
  const TransitionModelInterface& model(const std::string& name) const
  {
    return mixedModel().model(name);
  }
  std::shared_ptr<const TransitionModelInterface> getModel(const std::string& name) const
  {
    return mixedModel().getModel(name);
  }

  /**
   * @brief Returns the vector of numbers of the submodels in the
   * mixture that match a description.
   */
  Vuint getSubmodelNumbers(const std::string& desc) const;

  std::shared_ptr<const TransitionModelInterface> getTransitionModel() const { return pmixmodel_; }
  
  const TransitionModelInterface& transitionModel() const { return *pmixmodel_; }

  std::shared_ptr<const MixedTransitionModelInterface> getMixedModel() const { return pmixmodel_; }
  
  const MixedTransitionModelInterface& mixedModel() const { return *pmixmodel_; }

protected:
  TransitionModelInterface& transitionModel()
  {
    return *pmixmodel_;
  }

  MixedTransitionModelInterface& mixedModel() { return *pmixmodel_; }

  std::shared_ptr<const FrequencySetInterface> getFrequencySet() const
  {
    return pmixmodel_->getNModel(0)->getFrequencySet();
  }
  
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ABSTRACTBIBLIOMIXEDTRANSITIONMODEL_H
