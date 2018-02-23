//
// File: FromMixtureSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: samedi 24 octobre 2015, à 18h 28
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

#ifndef _FROM_MIXTURE_SUBSTITUTIONMODEL_H_
#define _FROM_MIXTURE_SUBSTITUTIONMODEL_H_

#include "AbstractSubstitutionModel.h"
#include "MixedSubstitutionModel.h"
#include "AbstractWrappedModel.h"

namespace bpp
{
/**
 * @brief Model taken from a SubModel of a
 * Mixture of SubstitutionModels.
 *
 * It has the same parameters as the SubModel.
 */
  class FromMixtureSubstitutionModel :
    public virtual AbstractTotallyWrappedSubstitutionModel,
    public AbstractParameterAliasable
  {
  private:
    /*
     * @brief The subModel taken from the MixtureOfSubstitutionModels.
     *
     * This subModel is normalized, even if it is not in the mixture.
     *
     */

    std::unique_ptr<SubstitutionModel> subModel_;

    /*
     * @brief The name of the mixture model (for io purpose).
     */

    std::string mixtName_;

  public:
    FromMixtureSubstitutionModel(const MixedSubstitutionModel& mixedModel, const std::string& subModelName, const std::string& mixtDesc);

    FromMixtureSubstitutionModel(const MixedSubstitutionModel& mixedModel, size_t subModelNumber, const std::string& mixtDesc);

    FromMixtureSubstitutionModel(const FromMixtureSubstitutionModel& fmsm);

    FromMixtureSubstitutionModel& operator=(const FromMixtureSubstitutionModel& fmsm);

    ~FromMixtureSubstitutionModel() {}

    FromMixtureSubstitutionModel* clone() const { return new FromMixtureSubstitutionModel(*this); }

  public:
    const SubstitutionModel& getSubstitutionModel() const
    {
      return *subModel_.get();
    }

  protected:
    SubstitutionModel& getSubstitutionModel()
    {
      return *subModel_.get();
    }

  public:
    /*
     *@ brief Methods to supersede AbstractSubstitutionModel methods.
     *
     * @{
     */

    /**
     * @brief Tells the model that a parameter value has changed.
     *
     * This updates the matrices consequently.
     */
    void fireParameterChanged(const ParameterList& parameters)
    {
      getModel().matchParametersValues(parameters);
    }

    virtual void setNamespace(const std::string& name)
    {
      AbstractParameterAliasable::setNamespace(name);
      getModel().setNamespace(name);
    }

    virtual void addRateParameter()
    {
      getModel().addRateParameter();
      addParameter_(new Parameter(getNamespace() + "rate", getModel().getRate(), &Parameter::R_PLUS_STAR)); 
    }
    
    /*
     * @}
     */
    
    std::string getName() const
    {
      size_t posp = mixtName_.find("(");
      
      return mixtName_.substr(0, posp) + "_" + getModel().getName() + mixtName_.substr(posp);
    }

  protected:

    void updateMatrices() 
    {
    }
    
    
  };
} // end of namespace bpp.

#endif  // _FROM_MIXTURE_SUBSTITUTIONMODEL_H_
