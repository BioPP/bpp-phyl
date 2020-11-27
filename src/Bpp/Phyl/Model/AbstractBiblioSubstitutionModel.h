//
// File: AbstractBiblioSubstitutionModel.h
// Created by: Laurent Guéguen
// Created on: vendredi 8 juillet 2011, à 20h 17
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

#ifndef _ABSTRACT_BIBLIO_SUBSTITUTIONMODEL_H_
#define _ABSTRACT_BIBLIO_SUBSTITUTIONMODEL_H_

#include "SubstitutionModel.h"
#include "AbstractWrappedModel.h"

#include <Bpp/Numeric/AbstractParameterAliasable.h>

namespace bpp
{
/**
 * @brief Partial implementation of the SubstitutionModel interface
 *   for models that are set for matching the bibliography, and are
 *   only defined through a link to a "real" model.
 *
 */

  class AbstractBiblioTransitionModel :
    public virtual AbstractTotallyWrappedTransitionModel,
    public AbstractParameterAliasable
  {
  protected:
    /**
     * @brief Tools to make the link between the Parameters of the
     * object and those of pmixmodel_.
     *
     */

    std::map<std::string, std::string> mapParNamesFromPmodel_;

    ParameterList lParPmodel_;

  public:
    AbstractBiblioTransitionModel(const std::string& prefix);

    AbstractBiblioTransitionModel(const AbstractBiblioTransitionModel& model);

    AbstractBiblioTransitionModel& operator=(const AbstractBiblioTransitionModel& model);

    virtual ~AbstractBiblioTransitionModel() {}

    virtual AbstractBiblioTransitionModel* clone() const = 0;

  protected:
    virtual void updateMatrices();

  public:

    
    /*
     * @brief get the name of a parameter from its name in a submodel
     *
     * @var name the name of the parameter in the submodel
     *
     */
    
    std::string getParNameFromPmodel(const std::string& name) const;
    
    /*
     * @brief get the name of a parameter in the submodel from its apparent name
     *
     * @var name the name of the parameter
     *
     */
    
    std::string getPmodelParName(const std::string& name) const;
    
    /*
     *@ brief Methods to supersede TransitionModel methods.
     *
     * @{
     */

    void addRateParameter();

    void setFreqFromData(const SequencedValuesContainer& data, double pseudoCount = 0);

    void setFreq(std::map<int, double>& frequ);

    /*
     * @}
     *
     */

    /*
     *@ brief Methods to supersede AbstractTransitionModel methods.
     *
     * @{
     */
    
    /**
     * @brief Tells the model that a parameter value has changed.
     *
     * This updates the matrices consequently.
     */
    virtual void fireParameterChanged(const ParameterList& parameters)
    {
      if (parameters.hasParameter(getNamespace()+"rate"))
      {
        getModel().setRate(parameters.getParameterValue(getNamespace()+"rate"));
        if (parameters.size()!=1)
          updateMatrices();
      }
      else
        updateMatrices();      
    }

    void setNamespace(const std::string& name);
  
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
      AbstractTotallyWrappedSubstitutionModel()
    {}

      AbstractBiblioSubstitutionModel(const AbstractBiblioSubstitutionModel& model):
      AbstractBiblioTransitionModel(model)
    {}
    

    AbstractBiblioSubstitutionModel& operator=(const AbstractBiblioSubstitutionModel& model)
    {
      AbstractBiblioTransitionModel::operator=(model);
      return *this;
    }      

    virtual ~AbstractBiblioSubstitutionModel() {}

    virtual AbstractBiblioSubstitutionModel* clone() const = 0;

    void updateMatrices()
    {
      AbstractBiblioTransitionModel::updateMatrices();
    }

  };


} // end of namespace bpp.


#endif  // _ABSTRACTBIBLIOSUBSTITUTIONMODEL_H_

