// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_MARKOVMODULATEDSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_MARKOVMODULATEDSUBSTITUTIONMODEL_H

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "SubstitutionModel.h"
#include "AbstractSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Partial implementation of the Markov-modulated class of
 * substitution models.
 *
 * This class wraps a substitution model and provide a matrix
 * describing rate changes. The rate matrix must be initialized by
 * derived classes of this class. Using these matrices, the
 * diagonalization procedure of Galtier and Jean-Marie is used.
 *
 * Such models can be described using two matrices:
 * a substitution matrix, @f$M@f$, with size @f$m@f$, which is a
 * "standard" substitution model of any alphabet type,
 * and a rate matrix @f$G@f$ of size @f$g@f$.
 * The generator of the markov-modulated model, @f$Q@f$ can be
 * written using Kronecker matrix operands.
 * @f[
 * Q=D_R \otimes M + G \otimes I_m,
 * @f]
 * where @f$D_R@f$ is the diagonal matrix of all rates, and
 * @f$I_m@f$ is the identity matrix of size @f$m@f$.
 *
 * This generator is normalized so that branch lengths are measured
 * in unit of mean number of substitutions per site, where
 * susbstitution here means "change of alphabet
 * Rate changes are not counted.
 *
 * Galtier N. and Jean-Marie A., Markov-modulated Markov chains and
 * the covarion process of molecular evolution (2004). _Journal of
 * Computational Biology_, 11:727-33.
 */
class MarkovModulatedSubstitutionModel :
  public virtual ReversibleSubstitutionModelInterface,
  public AbstractParameterAliasable,
  public AbstractLkTransitionModel
{
protected:
  std::unique_ptr<ReversibleSubstitutionModelInterface> model_;
  std::shared_ptr<const MarkovModulatedStateMap> stateMap_;
  size_t nbStates_; // Number of states in model
  size_t nbRates_; // Number of rate classes

  /**
   * @name Rate generator.
   *
   * These variables must be initialized in the constructor of the derived class.
   * @{
   */
  RowMatrix<double> rates_;                // All rates values
  RowMatrix<double> ratesExchangeability_; // All rates transitions
  Vdouble ratesFreq_;           // All rates equilibrium frequencies
  /**@}*/
  RowMatrix<double> ratesGenerator_;       // All rates transitions

  /**
   * @brief The generator matrix \f$Q\f$ of the model.
   */
  RowMatrix<double> generator_;

  /**
   * @brief The exchangeability matrix \f$S\f$ of the model.
   */
  RowMatrix<double> exchangeability_;

  /**
   * @brief The \f$U\f$ matrix made of left eigen vectors (by row).
   */
  RowMatrix<double> leftEigenVectors_;

  /**
   * @brief The \f$U^-1\f$ matrix made of right eigen vectors (by column).
   */
  RowMatrix<double> rightEigenVectors_;

  /**
   * @brief The vector of real parts of eigen values.
   */
  Vdouble eigenValues_;

  /**
   * @brief The vector of imaginary parts of the eigen values (zero
   * in case of reversible pmodel).
   */
  Vdouble iEigenValues_;

  /**
   * @brief Tell if the eigen decomposition should be performed.
   */
  bool eigenDecompose_;

  /**
   * @brief Tell if the equilibrium frequencies  should be computed
   * from the generator
   */

  bool compFreq_;

  /**
   * @brief These ones are for bookkeeping:
   */
  mutable RowMatrix<double> pijt_;
  mutable RowMatrix<double> dpijt_;
  mutable RowMatrix<double> d2pijt_;

  /**
   * @brief The vector of equilibrium frequencies.
   */
  Vdouble freq_;

  bool normalizeRateChanges_;

  std::string nestedPrefix_;

public:
  /**
   * @brief Build a new MarkovModulatedSubstitutionModel object.
   *
   * @param model The substitution model to use. Can be of any alphabet type, and will be owned by this instance.
   * @param nbRates The number of rate classes
   * @param normalizeRateChanges Tells if the branch lengths must be computed in terms of rate and state
   * NB: In most cases, this parameter should be set to false.
   * @param prefix The parameter namespace to be forwarded to the AbstractParametrizable constructor.
   * changes instead of state change only.
   */
  MarkovModulatedSubstitutionModel(
      std::unique_ptr<ReversibleSubstitutionModelInterface> model,
      unsigned int nbRates,
      bool normalizeRateChanges,
      const std::string& prefix) :
    AbstractParameterAliasable(prefix),
    model_(std::move(model)),
    stateMap_(std::make_shared<MarkovModulatedStateMap>(model_->stateMap(), nbRates)),
    nbStates_(model_->getNumberOfStates()),
    nbRates_(nbRates),
    rates_(nbRates, nbRates),
    ratesExchangeability_(nbRates, nbRates),
    ratesFreq_(nbRates),
    ratesGenerator_(nbRates, nbRates),
    generator_(),
    exchangeability_(),
    leftEigenVectors_(),
    rightEigenVectors_(),
    eigenValues_(),
    iEigenValues_(),
    eigenDecompose_(true),
    compFreq_(false),
    pijt_(), dpijt_(), d2pijt_(), freq_(),
    normalizeRateChanges_(normalizeRateChanges),
    nestedPrefix_("model_" + model_->getNamespace())
  {
    model_->setNamespace(prefix + nestedPrefix_);
    addParameters_(model_->getIndependentParameters());
  }

  MarkovModulatedSubstitutionModel(const MarkovModulatedSubstitutionModel& model);
  MarkovModulatedSubstitutionModel& operator=(const MarkovModulatedSubstitutionModel& model);

  virtual ~MarkovModulatedSubstitutionModel() {}

  MarkovModulatedSubstitutionModel* clone() const override = 0;

public:
  const Alphabet& alphabet() const override { return model_->alphabet(); }
  
  std::shared_ptr<const Alphabet> getAlphabet() const override { return model_->getAlphabet(); }

  size_t getNumberOfStates() const override { return stateMap_->getNumberOfModelStates(); }

  const StateMapInterface& stateMap() const override { return *stateMap_; }

  std::shared_ptr<const StateMapInterface> getStateMap() const override { return stateMap_; }

  const std::vector<int>& getAlphabetStates() const override { return stateMap_->getAlphabetStates(); }

  std::string getAlphabetStateAsChar(size_t index) const override { return stateMap_->getAlphabetStateAsChar(index); }

  int getAlphabetStateAsInt(size_t index) const override { return stateMap_->getAlphabetStateAsInt(index); }

  std::vector<size_t> getModelStates(int code) const override { return stateMap_->getModelStates(code); }

  std::vector<size_t> getModelStates(const std::string& code) const override { return stateMap_->getModelStates(code); }

  const Vdouble& getFrequencies() const override { return freq_; }

  const Matrix<double>& exchangeabilityMatrix() const override { return exchangeability_; }

  const Matrix<double>& generator() const override { return generator_; }

  const Matrix<double>& getPij_t(double t) const override;
  const Matrix<double>& getdPij_dt(double t) const override;
  const Matrix<double>& getd2Pij_dt2(double t) const override;

  const Vdouble& getEigenValues() const override { return eigenValues_; }
  const Vdouble& getIEigenValues() const override { return iEigenValues_; }

  bool isDiagonalizable() const override { return true; }
  bool isNonSingular() const override { return true; }

  const Matrix<double>& getRowLeftEigenVectors() const override { return leftEigenVectors_; }
  const Matrix<double>& getColumnRightEigenVectors() const override { return rightEigenVectors_; }

  double freq(size_t i) const override { return freq_[i]; }
  double Sij(size_t i, size_t j) const override { return exchangeability_(i, j); }
  double Qij(size_t i, size_t j) const override { return generator_(i, j); }

  double Pij_t    (size_t i, size_t j, double t) const override { return getPij_t(t)(i, j); }
  double dPij_dt  (size_t i, size_t j, double t) const override { return getdPij_dt(t)(i, j); }
  double d2Pij_dt2(size_t i, size_t j, double t) const override { return getd2Pij_dt2(t)(i, j); }

  double getInitValue(size_t i, int state) const override;

  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override
  {
    model_->setFreqFromData(data, pseudoCount);
    updateMatrices_();
  }

  void setFreq(std::map<int, double>& frequencies) override
  {
    model_->setFreq(frequencies);
    updateMatrices_();
  }

  const FrequencySetInterface& frequencySet() const override {
    throw NullPointerException("MarkovModulatedSubstitutionModel::frequencySet. No FrequencySet associated to this model. Frequencies are computed from the FrequencySet of the modulated model.");
  }
  
  const ReversibleSubstitutionModelInterface& nestedModel() const {
    return *model_;
  }

  /**
   * @brief Get the rate category corresponding to a particular state in the compound model.
   *
   * @param i The state.
   * @return The corresponding rate category.
   * @see getState;
   */
  size_t getRate(size_t i) const
  {
    return i / nbStates_;
  }

  double getRate() const override { return model_->getRate(); }

  void setRate(double rate) override { model_->setRate(rate); }

  bool isScalable() const override
  {
    return model_->isScalable();
  }

  void setScalable(bool scalable) override
  {
    model_->setScalable(scalable);
  }

  void normalize() override
  {
    model_->normalize();
    updateMatrices_();
  }

  void setDiagonal() override;

  double getScale() const override 
  {
    std::vector<double> v;
    MatrixTools::diag(generator_, v);
    return -VectorTools::scalar<double, double>(v, freq_);
  }

  void setScale(double scale) override
  {
    model_->setScale(scale);
    updateMatrices_();
  }

  void enableEigenDecomposition(bool yn) override { eigenDecompose_ = yn; }

  bool enableEigenDecomposition() override { return eigenDecompose_; }

  bool computeFrequencies() const override { return compFreq_; }

  void computeFrequencies(bool yn) override { compFreq_ = yn; }


  /**
   * @brief Tells the model that a parameter value has changed.
   *
   * This updates the matrices consequently.
   */
  virtual void fireParameterChanged(const ParameterList& parameters) override
  {
    model_->matchParametersValues(parameters);
    updateRatesModel_();
    updateMatrices_();
  }

  void setNamespace(const std::string& prefix) override;

protected:
  virtual void updateMatrices_();

  /**
   * @brief Update the rates vector, generator and equilibrium frequencies.
   *
   * This method must be implemented by the derived class.
   * It is called by the fireParameterChanged() method.
   */
  virtual void updateRatesModel_() = 0;

  Vdouble& getFrequencies_() override
  {
    return freq_;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_MARKOVMODULATEDSUBSTITUTIONMODEL_H
