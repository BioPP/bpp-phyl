//
// File: InMixtureSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: vendredi 22 septembre 2017, � 09h 57
//

/*
  Copyright or � or Copr. Bio++ Development Team, (November 16, 2004)

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
#include "MixedTransitionModel.h"
#include "AbstractWrappedModel.h"

namespace bpp
{
/**
 * @brief SubModel taken from a MixedTransitionModel, kept in the
 * context of the MixedTransitionModel (see
 * FromMixtureSubstitutionModel for an out of context subModel). So
 * "rate" and "scale" are set for the MixedTransitionModel.
 *
 * But method getRate returns the specific rate of the subModel.
 *
 * It owns the MixedTransitionModel.
 *
 * It has the same parameters as the MixedTransitionModel.
 */
  
  class InMixedSubstitutionModel :
    public virtual AbstractWrappedSubstitutionModel,
    public AbstractParameterAliasable
  {
  private:
    /*
     * @brief The MixedOfTransitionModels.
     *
     */

    std::unique_ptr<MixedTransitionModel> mixedModel_;

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
    InMixedSubstitutionModel(const MixedTransitionModel& mixedModel, const std::string& subModelName, const std::string& mixtDesc);

    InMixedSubstitutionModel(const MixedTransitionModel& mixedModel, size_t subModelNumber, const std::string& mixtDesc);

    InMixedSubstitutionModel(const InMixedSubstitutionModel& fmsm);

    InMixedSubstitutionModel& operator=(const InMixedSubstitutionModel& fmsm);

    InMixedSubstitutionModel* clone() const { return new InMixedSubstitutionModel(*this); }

  public:
    
    const MixedTransitionModel& getMixedModel() const
    {
      return *mixedModel_.get();
    }

    const SubstitutionModel& getSubstitutionModel() const
    {
      return dynamic_cast<const SubstitutionModel&>(*mixedModel_->getNModel(subModelNumber_));
    }

    size_t getSubModelNumber() const
    {
      return subModelNumber_;
    }

    bool computeFrequencies() const
    {
      return getMixedModel().computeFrequencies();
    }

    /**
     * @return Set if equilibrium frequencies should be computed from
     * the generator
     */
    
    void computeFrequencies(bool yn)
    {
      getMixedModel().computeFrequencies(yn);
    }

  protected:

    Vdouble& getFrequencies_()
    {
      return getMixedModel().getFrequencies_();
    }

    /*
     * @}
     *
     */

    MixedTransitionModel& getMixedModel()
    {
      return *mixedModel_.get();
    }

    SubstitutionModel& getSubstitutionModel()
    {
      return dynamic_cast<SubstitutionModel&>(*mixedModel_->getNModel(subModelNumber_));
    }

  public:

    
    /*
     *@ brief Methods to supersede WrappedSubstitutionModel methods.
     *
     * @{
     */

    double freq(size_t i) const { return getModel().freq(i); }

    double Pij_t    (size_t i, size_t j, double t) const { return getModel().Pij_t(i, j, t); }
    double dPij_dt  (size_t i, size_t j, double t) const { return getModel().dPij_dt (i, j, t); }
    double d2Pij_dt2(size_t i, size_t j, double t) const { return getModel().d2Pij_dt2(i, j, t); }

    const Vdouble& getFrequencies() const { return getModel().getFrequencies(); }

    const Matrix<double>& getPij_t(double t) const { return getModel().getPij_t(t); }

    const Matrix<double>& getdPij_dt(double t) const { return getModel().getdPij_dt(t); }

    const Matrix<double>& getd2Pij_dt2(double t) const { return getModel().getd2Pij_dt2(t); }

    double getInitValue(size_t i, int state) const
    {
      return getModel().getInitValue(i,state);
    }
    
    void setFreqFromData(const SequencedValuesContainer& data, double pseudoCount = 0)
    {
      getMixedModel().setFreqFromData(data, pseudoCount);
    }
    
    void setFreq(std::map<int, double>& frequencies)
    {
      getMixedModel().setFreq(frequencies);
    }

    /*
     *@ brief Methods to supersede SubstitutionModel methods.
     *
     * @{
     */

    double Qij(size_t i, size_t j) const { return getSubstitutionModel().Qij(i, j); }

    const Matrix<double>& getGenerator() const { return getSubstitutionModel().getGenerator(); }

    const Matrix<double>& getExchangeabilityMatrix() const { return getSubstitutionModel().getExchangeabilityMatrix(); }

    double Sij(size_t i, size_t j) const { return getSubstitutionModel().Sij(i, j); }

    void enableEigenDecomposition(bool yn) { getSubstitutionModel().enableEigenDecomposition(yn); }

    bool enableEigenDecomposition() { return getSubstitutionModel().enableEigenDecomposition(); }

    bool isDiagonalizable() const { return getSubstitutionModel().isDiagonalizable(); }

    bool isNonSingular() const { return getSubstitutionModel().isNonSingular(); }

    const Vdouble& getEigenValues() const { return getSubstitutionModel().getEigenValues(); }

    const Vdouble& getIEigenValues() const { return getSubstitutionModel().getIEigenValues(); }

    const Matrix<double>& getRowLeftEigenVectors() const { return getSubstitutionModel().getRowLeftEigenVectors(); }

    const Matrix<double>& getColumnRightEigenVectors() const { return getSubstitutionModel().getColumnRightEigenVectors(); }


    /*
     * @}
     *
     */

    bool isScalable() const 
    {
      return getSubstitutionModel().isScalable();
    }

    void setScalable(bool scalable)
    {
      getSubstitutionModel().setScalable(scalable);
    }

 
    double getScale() const { return getSubstitutionModel().getScale(); }

    void setScale(double scale) { getSubstitutionModel().setScale(scale); }


    void normalize()
    {
      getSubstitutionModel().normalize();
    }

    void setDiagonal()
    {
      getSubstitutionModel().setDiagonal();
    }

    double getRate() const
    {
      return getModel().getRate();
    }

    void setRate(double rate)
    {
      return getMixedModel().setRate(rate);
    }
    
    void addRateParameter()
    {
      getMixedModel().addRateParameter();
      addParameter_(new Parameter(getNamespace() + "rate", getMixedModel().getRate(), Parameter::R_PLUS_STAR));
    }

    /*
     * @}
     *
     */

    /*
     *@ brief Methods to supersede AbstractSubstitutionnModel methods.
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
      getMixedModel().matchParametersValues(parameters);
    }

    void setNamespace(const std::string& name)
    {
      AbstractParameterAliasable::setNamespace(name);
      getMixedModel().setNamespace(name);
    }


    /*
     * @}
     */

    std::string getName() const
    {
      return mixedModel_->getName();
    }

  protected:

    void updateMatrices() 
    {
    }
  };
  
} // end of namespace bpp.

#endif  // _IN_MIXED_SUBSTITUTIONMODEL_H_
