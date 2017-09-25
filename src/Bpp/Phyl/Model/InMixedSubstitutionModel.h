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
    public virtual SubstitutionModel,
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

    ~InMixedSubstitutionModel() {}

    InMixedSubstitutionModel* clone() const { return new InMixedSubstitutionModel(*this); }

  public:
    const MixedSubstitutionModel& getMixedModel() const
    {
      return *mixedModel_.get();
    }

    const SubstitutionModel& getSubModel() const
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

  public:
    /*
     *@ brief Methods to supersede SubstitutionModel methods.
     *
     * @{
     */

    const std::vector<int>& getAlphabetStates() const { return getSubModel().getAlphabetStates(); }

    const StateMap& getStateMap() const { return getSubModel().getStateMap(); }

    int getAlphabetStateAsInt(size_t i) const { return getSubModel().getAlphabetStateAsInt(i); }

    std::string getAlphabetStateAsChar(size_t i) const { return getSubModel().getAlphabetStateAsChar(i); }

    std::vector<size_t> getModelStates(int code) const { return getSubModel().getModelStates(code); }

    std::vector<size_t> getModelStates(const std::string& code) const { return getSubModel().getModelStates(code); }

    double freq(size_t i) const { return getSubModel().freq(i); }

    double Qij(size_t i, size_t j) const { return getSubModel().Qij(i, j); }

    double Pij_t    (size_t i, size_t j, double t) const { return getSubModel().Pij_t(i, j, t); }
    double dPij_dt  (size_t i, size_t j, double t) const { return getSubModel().dPij_dt (i, j, t); }
    double d2Pij_dt2(size_t i, size_t j, double t) const { return getSubModel().d2Pij_dt2(i, j, t); }

    const Vdouble& getFrequencies() const { return getSubModel().getFrequencies(); }

    const Matrix<double>& getGenerator() const { return getSubModel().getGenerator(); }

    const Matrix<double>& getExchangeabilityMatrix() const { return getSubModel().getExchangeabilityMatrix(); }

    double Sij(size_t i, size_t j) const { return getSubModel().Sij(i, j); }

    const Matrix<double>& getPij_t(double t) const { return getSubModel().getPij_t(t); }

    const Matrix<double>& getdPij_dt(double t) const { return getSubModel().getdPij_dt(t); }

    const Matrix<double>& getd2Pij_dt2(double t) const { return getSubModel().getd2Pij_dt2(t); }

    void enableEigenDecomposition(bool yn) { getMixedModel().enableEigenDecomposition(yn); }

    bool enableEigenDecomposition() { return getMixedModel().enableEigenDecomposition(); }

    bool isDiagonalizable() const { return getSubModel().isDiagonalizable(); }

    bool isNonSingular() const { return getSubModel().isNonSingular(); }

    const Vdouble& getEigenValues() const { return getSubModel().getEigenValues(); }

    const Vdouble& getIEigenValues() const { return getSubModel().getIEigenValues(); }

    const Matrix<double>& getRowLeftEigenVectors() const { return getSubModel().getRowLeftEigenVectors(); }

    const Matrix<double>& getColumnRightEigenVectors() const { return getSubModel().getColumnRightEigenVectors(); }

    double getRate() const { return getSubModel().getRate(); }

    void setRate(double rate) { return getMixedModel().setRate(rate); }

    void addRateParameter()
    {
      getMixedModel().addRateParameter();
      addParameter_(new Parameter(getNamespace() + "rate", getMixedModel().getRate(), &Parameter::R_PLUS_STAR));
    }

    void setFreqFromData(const SequenceContainer& data, double pseudoCount = 0){getMixedModel().setFreqFromData(data, pseudoCount); }

    void setFreq(std::map<int, double>& frequ) {getMixedModel().setFreq(frequ); }

    const Alphabet* getAlphabet() const { return getSubModel().getAlphabet(); }

    size_t getNumberOfStates() const { return getSubModel().getNumberOfStates(); }

    double getInitValue(size_t i, int state) const throw (BadIntException) { return getSubModel().getInitValue(i, state); }

    const FrequenciesSet* getFrequenciesSet() const {return getSubModel().getFrequenciesSet(); }

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

    std::string getName() const;

  };
} // end of namespace bpp.

#endif  // _IN_MIXED_SUBSTITUTIONMODEL_H_
