//
// File: InMixtureSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: vendredi 22 septembre 2017, à 09h 57
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

#ifndef _IN_MIXED_SUBSTITUTIONMODEL_H_
#define _IN_MIXED_SUBSTITUTIONMODEL_H_

#include "AbstractSubstitutionModel.h"
#include "MixedSubstitutionModel.h"
#include "WrappedModel.h"

namespace bpp
{
/**
 * @brief SubModel taken from a MixedSubstitutionModel, kept in the
 * context of the MixedSubstitutionModel (see
 * FromMixtureSubstitutionModel for an out of context subModel). So
 * "rate" and "scale" are set for the MixedSubstitutionModel.
 *
 * But method getRate returns the specific rate of the subModel.
 *
 * It owns the MixedSubstitutionModel.
 *
 * It has the same parameters as the MixedSubstitutionModel.
 */
  
  class InMixedSubstitutionModel :
    public virtual WrappedSubstitutionModel,
    public AbstractParameterAliasable
  {
  private:
    /*
     * @brief The MixedOfSubstitutionModels.
     *
     */

    std::unique_ptr<MixedSubstitutionModel> mixedModel_;

    /*
     * @brief the number of the submodel
     *
     */

    size_t subModelNumber_;
  
    /*
     * @brief The name of the mixture model (for io purpose).
     */

    std::string mixtName_;

  public:
    InMixedSubstitutionModel(const MixedSubstitutionModel& mixedModel, const std::string& subModelName, const std::string& mixtDesc);

    InMixedSubstitutionModel(const MixedSubstitutionModel& mixedModel, size_t subModelNumber, const std::string& mixtDesc);

    InMixedSubstitutionModel(const InMixedSubstitutionModel& fmsm);

    InMixedSubstitutionModel& operator=(const InMixedSubstitutionModel& fmsm);

    InMixedSubstitutionModel* clone() const { return new InMixedSubstitutionModel(*this); }

  public:
    const MixedSubstitutionModel& getMixedModel() const
    {
      return *mixedModel_.get();
    }

    const SubstitutionModel& getSubstitutionModel() const
    {
      return *mixedModel_->getNModel(subModelNumber_);
    }

    size_t getSubModelNumber() const
    {
      return subModelNumber_;
    }
    
  protected:
    MixedSubstitutionModel& getMixedModel()
    {
      return *mixedModel_.get();
    }

    SubstitutionModel& getSubstitutionModel()
    {
      return *mixedModel_->getNModel(subModelNumber_);
    }

  public:
    /*
     *@ brief Methods to supersede WrappedSubstitutionModel methods.
     *
     * @{
     */

    void enableEigenDecomposition(bool yn) { getMixedModel().enableEigenDecomposition(yn); }

    bool enableEigenDecomposition() { return getMixedModel().enableEigenDecomposition(); }

    void setRate(double rate) { return getMixedModel().setRate(rate); }

    void addRateParameter()
    {
      getMixedModel().addRateParameter();
      addParameter_(new Parameter(getNamespace() + "rate", getMixedModel().getRate(), &Parameter::R_PLUS_STAR));
    }

    void setFreqFromData(const SequencedValuesContainer& data, double pseudoCount = 0){getMixedModel().setFreqFromData(data, pseudoCount); }

    void setFreq(std::map<int, double>& frequ) {getMixedModel().setFreq(frequ); }

    /*
     * @}
     *
     */

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
      AbstractParameterAliasable::fireParameterChanged(parameters);
      getMixedModel().matchParametersValues(parameters);
    }

    void setNamespace(const std::string& name)
    {
      AbstractParameterAliasable::setNamespace(name);
      getMixedModel().setNamespace(name);
    }

  public:
    bool isScalable() const 
    {
      return getMixedModel().isScalable();
    }

    void setScalable(bool scalable)
    {
      getMixedModel().setScalable(scalable);
    }

    void normalize()
    {
      getMixedModel().normalize();
    }
    
    double getScale() const { return getMixedModel().getScale(); }

    void setScale(double scale) { getMixedModel().setScale(scale); }

    /*
     * @}
     */

    std::string getName() const
    {
      return mixedModel_->getName();
    }


  };
} // end of namespace bpp.

#endif  // _IN_MIXED_SUBSTITUTIONMODEL_H_
