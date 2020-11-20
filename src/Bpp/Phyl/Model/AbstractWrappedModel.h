//
// File: AbtractWrappedModel.h
// Created by: Laurent Guéguen
// Created on: mardi 26 septembre 2017, à 16h 18
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

#ifndef _ABSTRACT_WRAPPED_MODEL_H_
#define _ABSTRACT_WRAPPED_MODEL_H_

#include "WrappedModel.h"
#include <Bpp/Seq/Container/SequencedValuesContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

namespace bpp
{
/**
 * @brief Abstract class of Wrapping model class, where all methods
 * are redirected from getModel().
 * 
 *
 *
 */

  class AbstractWrappedModel :
    public virtual WrappedModel
  {
  public:
    AbstractWrappedModel() {}
    virtual ~AbstractWrappedModel() {}
    
  public:
    /*
     *@ brief Methods to supersede TransitionModel methods.
     *
     * @{
     */

    const std::vector<int>& getAlphabetStates() const { return getModel().getAlphabetStates(); }

    const StateMap& getStateMap() const { return getModel().getStateMap(); }

    std::shared_ptr<const StateMap> shareStateMap() const { return getModel().shareStateMap(); }

    int getAlphabetStateAsInt(size_t i) const { return getModel().getAlphabetStateAsInt(i); }
  
    std::string getAlphabetStateAsChar(size_t i) const { return getModel().getAlphabetStateAsChar(i); }

    std::vector<size_t> getModelStates(int code) const { return getModel().getModelStates(code); }
  
    std::vector<size_t> getModelStates(const std::string& code) const { return getModel().getModelStates(code); }


    const Alphabet* getAlphabet() const { return getModel().getAlphabet(); }

    size_t getNumberOfStates() const { return getModel().getNumberOfStates(); }

    const std::shared_ptr<FrequencySet> getFrequencySet() const { return getModel().getFrequencySet();}

    /*
     * @}
     */

    virtual std::string getName() const
    {
      return getModel().getName();
    }

    /*
     * @}
     *
     */
  };
  
  class AbstractWrappedTransitionModel :
    public virtual AbstractWrappedModel,
    public virtual WrappedTransitionModel
  {
  protected:
    BranchModel& getModel()
    {
      return getTransitionModel();
    }

  public:
    const std::shared_ptr<FrequencySet> getFrequencySet() const { return getTransitionModel().getFrequencySet();}

    const BranchModel& getModel() const
    {
      return getTransitionModel();
    }
  };
  

  class AbstractTotallyWrappedTransitionModel :
    public virtual AbstractWrappedTransitionModel
  {
  public:
    AbstractTotallyWrappedTransitionModel() {}
    virtual ~AbstractTotallyWrappedTransitionModel() {}
    
  public:
    /*
     *@ brief Methods to supersede TransitionModel methods.
     *
     * @{
     */

    double freq(size_t i) const { return getTransitionModel().freq(i); }

    double Pij_t    (size_t i, size_t j, double t) const { return getTransitionModel().Pij_t(i, j, t); }
    double dPij_dt  (size_t i, size_t j, double t) const { return getTransitionModel().dPij_dt (i, j, t); }
    double d2Pij_dt2(size_t i, size_t j, double t) const { return getTransitionModel().d2Pij_dt2(i, j, t); }

    const Vdouble& getFrequencies() const { return getTransitionModel().getFrequencies(); }

    const Matrix<double>& getPij_t(double t) const { return getTransitionModel().getPij_t(t); }

    const Matrix<double>& getdPij_dt(double t) const { return getTransitionModel().getdPij_dt(t); }

    const Matrix<double>& getd2Pij_dt2(double t) const { return getTransitionModel().getd2Pij_dt2(t); }

    double getInitValue(size_t i, int state) const
    {
      return getTransitionModel().getInitValue(i,state);
    }
    
    double getRate() const
    {
      return getTransitionModel().getRate();
    }

    void setRate(double rate)
    {
      return getTransitionModel().setRate(rate);
    }

    void setFreqFromData(const SequencedValuesContainer& data, double pseudoCount = 0)
    {
      std::map<int, double> freqs;
      SequenceContainerTools::getFrequencies(data, freqs, pseudoCount);
      // Re-compute generator and eigen values:
      getTransitionModel().setFreq(freqs);
    }
    
    void setFreq(std::map<int, double>& frequencies)
    {
      getTransitionModel().setFreq(frequencies);
    }

    bool computeFrequencies() const
    {
      return getTransitionModel().computeFrequencies();
    }

    /**
     * @return Set if equilibrium frequencies should be computed from
     * the generator
     */
    
    void computeFrequencies(bool yn)
    {
      getTransitionModel().computeFrequencies(yn);
    }

    /*
     * @}
     *
     */

  protected:

    Vdouble& getFrequencies_()
    {
      return getTransitionModel().getFrequencies_();
    }

  };
  
    
  class AbstractWrappedSubstitutionModel :
    public virtual AbstractWrappedTransitionModel,
    public virtual WrappedSubstitutionModel
  {
  public:
    AbstractWrappedSubstitutionModel() {}
    
    virtual ~AbstractWrappedSubstitutionModel() {}
    
    const TransitionModel& getTransitionModel() const
    {
      return getSubstitutionModel();
    }

  protected:
    TransitionModel& getTransitionModel()
    {
      return getSubstitutionModel();
    }

  };

   class AbstractTotallyWrappedSubstitutionModel :
    public virtual AbstractTotallyWrappedTransitionModel,
    public virtual AbstractWrappedSubstitutionModel
  {
  public:
    AbstractTotallyWrappedSubstitutionModel() {}
    
    virtual ~AbstractTotallyWrappedSubstitutionModel() {}
    
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

    void normalize()
    {
      getSubstitutionModel().normalize();
    }

    void setDiagonal()
    {
      getSubstitutionModel().setDiagonal();
    }

    double getScale() const { return getSubstitutionModel().getScale(); }

    void setScale(double scale) { getSubstitutionModel().setScale(scale); }

    /*
     * @}
     */
  };


  
} // end of namespace bpp.


#endif  // _ABSTRACT_WRAPPED_MODEL_SUBSTITUTIONMODEL_H_

