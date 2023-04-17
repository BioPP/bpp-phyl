//
// File: AbstractWrappedModel.h
// Authors:
//   Laurent Guéguen
// Created: mardi 26 septembre 2017, ÃÂ  16h 18
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_MODEL_ABSTRACTWRAPPEDMODEL_H
#define BPP_PHYL_MODEL_ABSTRACTWRAPPEDMODEL_H

#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Container/SequenceData.h>

#include "WrappedModel.h"
#include "AbstractSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Abstract class of Wrapping model class, where all methods
 * are redirected from getModel().
 */
class AbstractWrappedModel :
  public virtual AbstractParameterAliasable,
  public virtual WrappedModelInterface
{
public:
  AbstractWrappedModel(const std::string& prefix) :
    AbstractParameterAliasable(prefix)
  {}
  
  virtual ~AbstractWrappedModel() {}

public:
  /**
   * @ brief Methods to supersede TransitionModel methods.
   *
   * @{
   */
  const std::vector<int>& getAlphabetStates() const override { return model().getAlphabetStates(); }

  const StateMapInterface& stateMap() const override { return model().stateMap(); }

  std::shared_ptr<const StateMapInterface> getStateMap() const override { return model().getStateMap(); }

  int getAlphabetStateAsInt(size_t i) const override { return model().getAlphabetStateAsInt(i); }

  std::string getAlphabetStateAsChar(size_t i) const override { return model().getAlphabetStateAsChar(i); }

  std::vector<size_t> getModelStates(int code) const override { return model().getModelStates(code); }

  std::vector<size_t> getModelStates(const std::string& code) const override { return model().getModelStates(code); }


  const Alphabet& alphabet() const override { return model().alphabet(); }
  
  std::shared_ptr<const Alphabet> getAlphabet() const override { return model().getAlphabet(); }

  size_t getNumberOfStates() const override { return model().getNumberOfStates(); }

  const FrequencySetInterface& frequencySet() const override { return model().frequencySet(); }
  
  /**
   * @}
   */
  virtual std::string getName() const override
  {
    return model().getName();
  }

  /**
   * @}
   */
};

class AbstractWrappedTransitionModel :
  public virtual AbstractWrappedModel,
  public virtual AbstractLkTransitionModel,
  public virtual WrappedTransitionModelInterface
{
public:
  AbstractWrappedTransitionModel(const std::string& prefix):
      AbstractWrappedModel(prefix)
  {}
  
protected:
  BranchModelInterface& model_()
  {
    return transitionModel_();
  }
  
  virtual TransitionModelInterface& transitionModel_() = 0;

public:
  const FrequencySetInterface& frequencySet() const override
  {
    return transitionModel().frequencySet();
  }

  const BranchModelInterface& model() const override
  {
    return transitionModel();
  }
  
};


class AbstractTotallyWrappedTransitionModel :
  public virtual AbstractWrappedTransitionModel
{
public:
  AbstractTotallyWrappedTransitionModel(const std::string& prefix):
      AbstractWrappedTransitionModel(prefix) {}
  
  virtual ~AbstractTotallyWrappedTransitionModel() {}

public:
  /**
   * @brief Methods to supersede TransitionModel methods.
   *
   * @{
   */
  double freq(size_t i) const override { return transitionModel().freq(i); }

  double Pij_t    (size_t i, size_t j, double t) const override { return transitionModel().Pij_t(i, j, t); }
  double dPij_dt  (size_t i, size_t j, double t) const override { return transitionModel().dPij_dt (i, j, t); }
  double d2Pij_dt2(size_t i, size_t j, double t) const override { return transitionModel().d2Pij_dt2(i, j, t); }

  const Vdouble& getFrequencies() const override { return transitionModel().getFrequencies(); }

  const Matrix<double>& getPij_t(double t) const override { return transitionModel().getPij_t(t); }

  const Matrix<double>& getdPij_dt(double t) const override { return transitionModel().getdPij_dt(t); }

  const Matrix<double>& getd2Pij_dt2(double t) const override { return transitionModel().getd2Pij_dt2(t); }

  double getInitValue(size_t i, int state) const override
  {
    return transitionModel().getInitValue(i, state);
  }

  double getRate() const override
  {
    return transitionModel().getRate();
  }

  void setRate(double rate) override
  {
    return transitionModel_().setRate(rate);
  }

  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override
  {
    std::map<int, double> freqs;
    SequenceContainerTools::getFrequencies(data, freqs, pseudoCount);
    // Re-compute generator and eigen values:
    transitionModel_().setFreq(freqs);
  }

  void setFreq(std::map<int, double>& frequencies) override
  {
    transitionModel_().setFreq(frequencies);
  }

  bool computeFrequencies() const override
  {
    return transitionModel().computeFrequencies();
  }

  /**
   * @return Set if equilibrium frequencies should be computed from
   * the generator
   */
  void computeFrequencies(bool yn) override
  {
    transitionModel_().computeFrequencies(yn);
  }

  /**
   * @}
   */

protected:
  Vdouble& getFrequencies_() override
  {
    return transitionModel_().getFrequencies_();
  }
};


class AbstractWrappedSubstitutionModel :
  public virtual AbstractWrappedTransitionModel,
  public virtual WrappedSubstitutionModelInterface
{
public:
  AbstractWrappedSubstitutionModel(const std::string& prefix) :
    AbstractWrappedTransitionModel(prefix)
  {}

  virtual ~AbstractWrappedSubstitutionModel() {}

  const TransitionModelInterface& transitionModel() const
  {
    return substitutionModel();
  }

protected:
  TransitionModelInterface& transitionModel_()
  {
    return substitutionModel_();
  }

  virtual SubstitutionModelInterface& substitutionModel_() = 0;
};

class AbstractTotallyWrappedSubstitutionModel :
  public virtual AbstractTotallyWrappedTransitionModel,
  public virtual AbstractWrappedSubstitutionModel
{
public:
  AbstractTotallyWrappedSubstitutionModel(const std::string& prefix):
       AbstractTotallyWrappedTransitionModel(prefix),
       AbstractWrappedSubstitutionModel(prefix)
  {}

  virtual ~AbstractTotallyWrappedSubstitutionModel() {}

  /**
   * @brief Methods to supersede SubstitutionModel methods.
   *
   * @{
   */
  double Qij(size_t i, size_t j) const { return substitutionModel().Qij(i, j); }

  const Matrix<double>& getGenerator() const { return substitutionModel().getGenerator(); }

  const Matrix<double>& getExchangeabilityMatrix() const { return substitutionModel().getExchangeabilityMatrix(); }

  double Sij(size_t i, size_t j) const { return substitutionModel().Sij(i, j); }

  void enableEigenDecomposition(bool yn) { substitutionModel_().enableEigenDecomposition(yn); }

  bool enableEigenDecomposition() { return substitutionModel_().enableEigenDecomposition(); }

  bool isDiagonalizable() const { return substitutionModel().isDiagonalizable(); }

  bool isNonSingular() const { return substitutionModel().isNonSingular(); }

  const Vdouble& getEigenValues() const { return substitutionModel().getEigenValues(); }

  const Vdouble& getIEigenValues() const { return substitutionModel().getIEigenValues(); }

  const Matrix<double>& getRowLeftEigenVectors() const { return substitutionModel().getRowLeftEigenVectors(); }

  const Matrix<double>& getColumnRightEigenVectors() const { return substitutionModel().getColumnRightEigenVectors(); }


  /**
   * @}
   */
  bool isScalable() const
  {
    return substitutionModel().isScalable();
  }

  void setScalable(bool scalable)
  {
    substitutionModel_().setScalable(scalable);
  }

  void normalize()
  {
    substitutionModel_().normalize();
  }

  void setDiagonal()
  {
    substitutionModel_().setDiagonal();
  }

  double getScale() const
  {
    return substitutionModel().getScale();
  }

  void setScale(double scale)
  {
    substitutionModel_().setScale(scale);
  }

  /**
   * @}
   */
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ABSTRACTWRAPPEDMODEL_H
