//
// File: FromMixtureSubstitutionModel.h
// Authors:
//   Laurent Gueguen
// Created: samedi 24 octobre 2015, Ã  18h 28
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

#ifndef BPP_PHYL_MODEL_FROMMIXTURESUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_FROMMIXTURESUBSTITUTIONMODEL_H


#include "AbstractSubstitutionModel.h"
#include "AbstractWrappedModel.h"
#include "MixedTransitionModel.h"

namespace bpp
{
/**
 * @brief Model taken from a SubModel of a
 * Mixture of SubstitutionModels.
 *
 * It has the same parameters as the SubModel.
 */
class FromMixtureSubstitutionModel :
  public virtual AbstractTotallyWrappedSubstitutionModel
{
private:
  /**
   * @brief The subModel taken from the AbstractTotallyWrappedSubstitutionModel.
   *
   * This subModel is normalized, even if it is not in the mixture.
   */
  std::unique_ptr<SubstitutionModelInterface> subModel_;

  /**
   * @brief The name of the mixture model (for io purpose).
   */
  std::string mixtName_;

public:
  FromMixtureSubstitutionModel(
      const MixedTransitionModelInterface& mixedModel,
      const std::string& subModelName,const std::string& mixtDesc);

  FromMixtureSubstitutionModel(
      const MixedTransitionModelInterface& mixedModel,
      size_t subModelNumber,
      const std::string& mixtDesc);

  FromMixtureSubstitutionModel(const FromMixtureSubstitutionModel& fmsm);

  FromMixtureSubstitutionModel& operator=(const FromMixtureSubstitutionModel& fmsm);

  virtual ~FromMixtureSubstitutionModel() {}

  FromMixtureSubstitutionModel* clone() const override { return new FromMixtureSubstitutionModel(*this); }

public:
  const SubstitutionModelInterface& substitutionModel() const override
  {
    return *subModel_.get();
  }

protected:
  SubstitutionModelInterface& substitutionModel() override
  {
    return *subModel_;
  }

public:
  /**
   * @ 
   *brief Methods to supersede AbstractSubstitutionModel methods.
   *
   * @{
   */

  /**
   * @brief Tells the model that a parameter value has changed.
   *
   * This updates the matrices consequently.
   */
  void fireParameterChanged(const ParameterList& parameters) override
  {
    model().matchParametersValues(parameters);
  }

  virtual void setNamespace(const std::string& name) override
  {
    AbstractParameterAliasable::setNamespace(name);
    model().setNamespace(name);
  }

  virtual void addRateParameter() override
  {
    model().addRateParameter();
    addParameter_(new Parameter(getNamespace() + "rate", model().getRate(), Parameter::R_PLUS_STAR));
  }

  /*
   * @}
   */
  std::string getName() const override
  {
    size_t posp = mixtName_.find("(");
    return mixtName_.substr(0, posp) + "_" + model().getName() + mixtName_.substr(posp);
  }

protected:
  void updateMatrices() override {}
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_FROMMIXTURESUBSTITUTIONMODEL_H
