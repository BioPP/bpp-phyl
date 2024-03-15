// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_ABSTRACTBIBLIOSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_ABSTRACTBIBLIOSUBSTITUTIONMODEL_H

#include <Bpp/Numeric/AbstractParameterAliasable.h>

#include "AbstractWrappedModel.h"
#include "SubstitutionModel.h"

namespace bpp
{
/**
 * @brief Partial implementation of the SubstitutionModel interface
 *   for models that are set for matching the bibliography, and are
 *   only defined through a link to a "real" model.
 */
class AbstractBiblioTransitionModel :
  public virtual AbstractTotallyWrappedTransitionModel
{
protected:
  /**
   * @brief Tools to make the link between the Parameters of the
   * object and those of pmixmodel_.
   */
  std::map<std::string, std::string> mapParNamesFromPmodel_;

  ParameterList lParPmodel_;

public:
  AbstractBiblioTransitionModel(const std::string& prefix);

  AbstractBiblioTransitionModel(const AbstractBiblioTransitionModel& model);

  AbstractBiblioTransitionModel& operator=(const AbstractBiblioTransitionModel& model);

  virtual ~AbstractBiblioTransitionModel() {}

protected:

  virtual void updateMatrices_();

public:

  /**
   * @brief get the name of a parameter from its name in a submodel
   *
   * @param name the name of the parameter in the submodel
   */
  std::string getParNameFromPmodel(const std::string& name) const;

  /**
   * @brief get the name of a parameter in the submodel from its apparent name
   *
   * @param name the name of the parameter
   */
  std::string getPmodelParName(const std::string& name) const;

  /**
   * @brief Methods to supersede TransitionModel methods.
   *
   * @{
   */
  void addRateParameter() override;

  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override;

  void setFreq(std::map<int, double>& frequ) override;

  /*
   * @}
   *
   */

  /**
   * @brief Methods to supersede AbstractTransitionModel methods.
   *
   * @{
   */

  /**
   * @brief Tells the model that a parameter value has changed.
   *
   * This updates the matrices consequently.
   */
  virtual void fireParameterChanged(const ParameterList& parameters) override
  {
    if (parameters.hasParameter(getNamespace() + "rate"))
    {
      model_().setRate(parameters.getParameterValue(getNamespace() + "rate"));
      if (parameters.size() != 1)
        updateMatrices_();
    }
    else
      updateMatrices_();
  }

  void setNamespace(const std::string& name) override;

  /*
   * @}
   */
};

class AbstractBiblioSubstitutionModel :
  public virtual AbstractBiblioTransitionModel,
  public virtual AbstractTotallyWrappedSubstitutionModel
{
	
public:
  AbstractBiblioSubstitutionModel(const std::string& prefix) :
    AbstractBiblioTransitionModel(prefix),
    AbstractTotallyWrappedSubstitutionModel(prefix)
  {}

  AbstractBiblioSubstitutionModel(const AbstractBiblioSubstitutionModel& model) :
    AbstractBiblioTransitionModel(model)
  {}


  AbstractBiblioSubstitutionModel& operator=(const AbstractBiblioSubstitutionModel& model)
  {
    AbstractBiblioTransitionModel::operator=(model);
    return *this;
  }

  virtual ~AbstractBiblioSubstitutionModel() {}

protected:

  void updateMatrices_() override
  {
    AbstractBiblioTransitionModel::updateMatrices_();
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ABSTRACTBIBLIOSUBSTITUTIONMODEL_H
