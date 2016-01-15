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
#include "MixtureOfSubstitutionModels.h"

namespace bpp
{
/**
 * @brief Model taken from a SubModel of a
 * MixtureOfSubstitutionModels.
 *
 * It has the same parameters as the SubModel.
 *
 */

  class FromMixtureSubstitutionModel :
    public virtual SubstitutionModel,
    public AbstractParameterAliasable
  {
  private:
    /*
     * @brief The subModel taken from the MixtureOfSubstitutionModels.
     *
     * This subModel is normalized, even if it is not in the mixture.
     *
     */
  
    std::auto_ptr<SubstitutionModel> subModel_;
  
    /*
     * @brief The name of the mixture model (for io purpose).
     */

    std::string mixtName_;

  public:
    FromMixtureSubstitutionModel(const MixedSubstitutionModel& mixedModel, const std::string& subModelName, const std::string& mixtDesc);

    FromMixtureSubstitutionModel(const FromMixtureSubstitutionModel& fmsm);

    FromMixtureSubstitutionModel& operator=(const FromMixtureSubstitutionModel& fmsm);
    
    ~FromMixtureSubstitutionModel() {}

    FromMixtureSubstitutionModel* clone() const { return new FromMixtureSubstitutionModel(*this); }
  
  public:
    const SubstitutionModel& getModel() const
    {
      return *subModel_.get();
    }

  protected:
    void updateMatrices();

    SubstitutionModel& getModel()
    {
      return *subModel_.get();
    }


  public:
    /*
     *@ brief Methods to supersede SubstitutionModel methods.
     *
     * @{
     */

    const std::vector<int>& getAlphabetStates() const { return getModel().getAlphabetStates(); }

    const StateMap& getStateMap() const { return getModel().getStateMap(); }

    int getAlphabetStateAsInt(size_t i) const { return getModel().getAlphabetStateAsInt(i); }
  
    std::string getAlphabetStateAsChar(size_t i) const { return getModel().getAlphabetStateAsChar(i); }

    std::vector<size_t> getModelStates(int code) const { return getModel().getModelStates(code); }
  
    std::vector<size_t> getModelStates(const std::string& code) const { return getModel().getModelStates(code); }

    double freq(size_t i) const { return getModel().freq(i); }

    double Qij(size_t i, size_t j) const { return getModel().Qij(i, j); }

    double Pij_t    (size_t i, size_t j, double t) const { return getModel().Pij_t(i, j, t); }
    double dPij_dt  (size_t i, size_t j, double t) const { return getModel().dPij_dt (i, j, t); }
    double d2Pij_dt2(size_t i, size_t j, double t) const { return getModel().d2Pij_dt2(i, j, t); }

    const Vdouble& getFrequencies() const { return getModel().getFrequencies(); }

    const Matrix<double>& getGenerator() const { return getModel().getGenerator(); }

    const Matrix<double>& getExchangeabilityMatrix() const { return getModel().getExchangeabilityMatrix(); }

    double Sij(size_t i, size_t j) const { return getModel().Sij(i, j); }

    const Matrix<double>& getPij_t(double t) const { return getModel().getPij_t(t); }

    const Matrix<double>& getdPij_dt(double t) const { return getModel().getdPij_dt(t); }

    const Matrix<double>& getd2Pij_dt2(double t) const { return getModel().getd2Pij_dt2(t); }

    void enableEigenDecomposition(bool yn) { getModel().enableEigenDecomposition(yn); }

    bool enableEigenDecomposition() { return getModel().enableEigenDecomposition(); }

    bool isDiagonalizable() const { return getModel().isDiagonalizable(); }

    bool isNonSingular() const { return getModel().isNonSingular(); }

    const Vdouble& getEigenValues() const { return getModel().getEigenValues(); }

    const Vdouble& getIEigenValues() const { return getModel().getIEigenValues(); }

    const Matrix<double>& getRowLeftEigenVectors() const { return getModel().getRowLeftEigenVectors(); }

    const Matrix<double>& getColumnRightEigenVectors() const { return getModel().getColumnRightEigenVectors(); }

    double getRate() const { return getModel().getRate(); }

    void setRate(double rate) { return getModel().setRate(rate); }

    void addRateParameter()
    {
      getModel().addRateParameter();
      addParameter_(new Parameter(getNamespace() + "rate", getModel().getRate(), &Parameter::R_PLUS_STAR));
    }

    void setFreqFromData(const SequenceContainer& data, double pseudoCount = 0){getModel().setFreqFromData(data, pseudoCount);}

    void setFreq(std::map<int, double>& frequ) {getModel().setFreq(frequ);}
    
    const Alphabet* getAlphabet() const { return getModel().getAlphabet(); }

    size_t getNumberOfStates() const { return getModel().getNumberOfStates(); }

    double getInitValue(size_t i, int state) const throw (BadIntException) { return getModel().getInitValue(i, state); }

    const FrequenciesSet* getFrequenciesSet() const {return getModel().getFrequenciesSet(); }

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
      getModel().matchParametersValues(parameters);
    }

    void setNamespace(const std::string& name)
    {
      AbstractParameterAliasable::setNamespace(name);
      getModel().setNamespace(name);
    }
    
  public:
    double getScale() const { return getModel().getScale(); }

    void setScale(double scale) { getModel().setScale(scale); }

    /*
     * @}
     */

    std::string getName() const;

  };
} // end of namespace bpp.

#endif  // _FROM_MIXTURE_SUBSTITUTIONMODEL_H_

