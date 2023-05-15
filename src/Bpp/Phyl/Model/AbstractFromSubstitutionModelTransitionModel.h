//
// File: AbstractFromSubstitutionModelTransitionModel.h
// Authors:
//   Laurent Gueguen
// Created: lundi 24 avril 2017, Ã  22h 55
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

#ifndef BPP_PHYL_MODEL_ABSTRACTFROMSUBSTITUTIONMODELTRANSITIONMODEL_H
#define BPP_PHYL_MODEL_ABSTRACTFROMSUBSTITUTIONMODELTRANSITIONMODEL_H


#include "AbstractWrappedModel.h"

namespace bpp
{
/**
 * @brief Virtual class of a Transition Model related to a given
 * SubstitutionModel.
 *
 * It has the same parameters as the SubModel.
 */
class AbstractFromSubstitutionModelTransitionModel :
  public virtual AbstractWrappedTransitionModel
{
protected:
  /**
   * @brief The related model.
   */
  std::unique_ptr<SubstitutionModelInterface> subModel_;

  /**
   * The number of states
   */
  size_t size_;

  /**
   * @brief These ones are for bookkeeping:
   */
  mutable RowMatrix<double> pijt_;
  mutable RowMatrix<double> dpijt_;
  mutable RowMatrix<double> d2pijt_;

  std::string nestedPrefix_;

public:
  AbstractFromSubstitutionModelTransitionModel(
    std::unique_ptr<SubstitutionModelInterface> subModel,
    const std::string& prefix);

  AbstractFromSubstitutionModelTransitionModel(const AbstractFromSubstitutionModelTransitionModel& fmsm);

  AbstractFromSubstitutionModelTransitionModel& operator=(const AbstractFromSubstitutionModelTransitionModel& fmsm);

  virtual ~AbstractFromSubstitutionModelTransitionModel() {}

public:
  const SubstitutionModelInterface& substitutionModel() const
  {
    return *subModel_;
  }

  const TransitionModelInterface& transitionModel() const override
  {
    return *subModel_;
  }

  const BranchModelInterface& model() const override
  {
    return *subModel_;
  }

  bool computeFrequencies() const override
  {
    return subModel_->computeFrequencies();
  }

  /**
   * @return Set if equilibrium frequencies should be computed from
   * the generator
   */
  void computeFrequencies(bool yn) override
  {
    subModel_->computeFrequencies(yn);
  }

  /**
   * @}
   */

protected:
  Vdouble& getFrequencies_() override
  {
    return subModel_->getFrequencies_();
  }

  SubstitutionModelInterface& substitutionModel_()
  {
    return *subModel_;
  }


  TransitionModelInterface& transitionModel_() override
  {
    return *subModel_;
  }

  BranchModelInterface& model_()
  {
    return *subModel_;
  }

public:
  virtual void addRateParameter() override
  {
    model_().addRateParameter();
    addParameter_(new Parameter(getNamespace() + "rate", model().getRate(), Parameter::R_PLUS_STAR));
  }

  virtual void fireParameterChanged(const ParameterList& parameters) override
  {
    AbstractParameterAliasable::fireParameterChanged(parameters);
    model_().matchParametersValues(parameters);
  }

  virtual void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix + nestedPrefix_);
    model_().setNamespace(prefix + nestedPrefix_);
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ABSTRACTFROMSUBSTITUTIONMODELTRANSITIONMODEL_H
