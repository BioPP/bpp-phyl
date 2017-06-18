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

#include "SubstitutionModel.h"

namespace bpp
{
/**
 * @brief Virtual class of a Transition Model related to a given
 * SubstitutionModel.
 *
 * It has the same parameters as the SubModel.
 */

  class AbstractFromSubstitutionModelTransitionModel :
    public virtual TransitionModel,
    public AbstractParameterAliasable
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
    const SubstitutionModel& getModel() const
    {
      return *subModel_.get();
    }

  protected:
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

    const Alphabet* getAlphabet() const { return getModel().getAlphabet(); }

    size_t getNumberOfStates() const { return getModel().getNumberOfStates(); }

    const std::vector<int>& getAlphabetStates() const { return getModel().getAlphabetStates(); }

    const StateMap& getStateMap() const { return getModel().getStateMap(); }

    int getAlphabetStateAsInt(size_t i) const { return getModel().getAlphabetStateAsInt(i); }

    std::string getAlphabetStateAsChar(size_t i) const { return getModel().getAlphabetStateAsChar(i); }

    std::vector<size_t> getModelStates(int code) const { return getModel().getModelStates(code); }

    std::vector<size_t> getModelStates(const std::string& code) const { return getModel().getModelStates(code); }

    virtual double freq(size_t i) const = 0;

    virtual double Pij_t (size_t i, size_t j, double t) const = 0;
    
    virtual double dPij_dt (size_t i, size_t j, double t) const = 0;
    
    virtual double d2Pij_dt2(size_t i, size_t j, double t) const = 0;

    virtual const Vdouble& getFrequencies() const = 0;

    virtual const Matrix<double>& getPij_t(double t) const = 0;

    virtual const Matrix<double>& getdPij_dt(double t) const = 0;

    virtual const Matrix<double>& getd2Pij_dt2(double t) const = 0;

    double getRate() const { return getModel().getRate(); }

    void setRate(double rate) { return getModel().setRate(rate); }

    void addRateParameter()
    {
      getModel().addRateParameter();
      addParameter_(new Parameter(getNamespace() + "rate", getModel().getRate(), &Parameter::R_PLUS_STAR));
    }

    void setFreqFromData(const SequenceContainer& data, double pseudoCount = 0){getModel().setFreqFromData(data, pseudoCount); }

    void setFreq(std::map<int, double>& frequ) {getModel().setFreq(frequ); }

    double getInitValue(size_t i, int state) const throw (BadIntException) { return getModel().getInitValue(i, state); }

    virtual const FrequenciesSet* getFrequenciesSet() const = 0;

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

    /*
     * @}
     */

    virtual std::string getName() const
    {
      return getModel().getName();
    }
    
  };
} // end of namespace bpp.

#endif  // _ABSTRACT_FROM_SUBSTITUTION_MODEL_TRANSITION_MODEL_H_
