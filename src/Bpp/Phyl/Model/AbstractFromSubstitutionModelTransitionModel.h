//
// File: AbstractFromModelTransitionModel.h
// Created by: Laurent Gueguen
// Created on: lundi 24 avril 2017, à 22h 55
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

#ifndef _ABSTRACT_FROM_SUBSTITUTION_MODEL_TRANSITION_MODEL_H_
#define _ABSTRACT_FROM_SUBSTITUTION_MODEL_TRANSITION_MODEL_H_

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
    public virtual AbstractWrappedModel,
    virtual public AbstractParameterAliasable
  {
  protected:
    /*
     * @brief The related model.
     *
     */
    
    std::unique_ptr<SubstitutionModel> subModel_;

    /*
     * The number of states
     *
     */
    
     size_t size_;

    /**
     * @brief These ones are for bookkeeping:
     */
    mutable RowMatrix<double> pij_t;
    mutable RowMatrix<double> dpij_t;
    mutable RowMatrix<double> d2pij_t;

  public:
    AbstractFromSubstitutionModelTransitionModel(const SubstitutionModel& subModel, const std::string& prefix);

    AbstractFromSubstitutionModelTransitionModel(const AbstractFromSubstitutionModelTransitionModel& fmsm);

    AbstractFromSubstitutionModelTransitionModel& operator=(const AbstractFromSubstitutionModelTransitionModel& fmsm);
    
    virtual ~AbstractFromSubstitutionModelTransitionModel() {};
    
    virtual AbstractFromSubstitutionModelTransitionModel* clone() const = 0;
    
  public:
    const SubstitutionModel& getSubstitutionModel() const
    {
      return *subModel_.get();
    }

    const TransitionModel& getModel() const
    {
      return *subModel_.get();
    }

    bool computeFrequencies() const
    {
      return subModel_->computeFrequencies();
    }

    /**
     * @return Set if equilibrium frequencies should be computed from
     * the generator
     */
    
    void computeFrequencies(bool yn)
    {
      subModel_->computeFrequencies(yn);
    }

    /*
     * @}
     *
     */

  protected:

    Vdouble& getFrequencies_()
    {
      return subModel_->getFrequencies_();
    }

    SubstitutionModel& getSubstitutionModel()
    {
      return *subModel_.get();
    }

    
    TransitionModel& getModel()
    {
      return *subModel_.get();
    }

  public:
    virtual void addRateParameter()
    {
      getModel().addRateParameter();
      addParameter_(new Parameter(getNamespace() + "rate", getModel().getRate(), &Parameter::R_PLUS_STAR));
    }

    virtual void fireParameterChanged(const ParameterList& parameters)
    {
      AbstractParameterAliasable::fireParameterChanged(parameters);
      getModel().matchParametersValues(parameters);
    }

    virtual void setNamespace(const std::string& name)
    {
      AbstractParameterAliasable::setNamespace(name);
      getModel().setNamespace(name);
    }

  };
} // end of namespace bpp.

#endif  // _ABSTRACT_FROM_SUBSTITUTION_MODEL_TRANSITION_MODEL_H_
