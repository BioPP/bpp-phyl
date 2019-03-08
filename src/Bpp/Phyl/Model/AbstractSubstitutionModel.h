//
// File: AbstractSubstitutionModel.h
// Created by: Julien Dutheil
// Created on: Tue May 27 10:31:49 2003
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

#ifndef _ABSTRACTSUBSTITUTIONMODEL_H_
#define _ABSTRACTSUBSTITUTIONMODEL_H_

#include "SubstitutionModel.h"

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/VectorTools.h>

#include <memory>

namespace bpp
{
/**
 * @brief Partial implementation of the SubstitutionModel interface.
 *
 * This abstract class provides some fields, namely:
 * - alphabet_: a pointer toward the alphabet,
 * - size_: the size of the alphabet, a parameter frequently called during various computations,
 * - rate_: the rate of the model
 * - generator_, leftEigenVectors_, rightEigenVectors_: useful matrices,
 * - eigenValues_, iEigenValues_, freq_: useful vectors.
 * - isDiagonalizable_ : boolean value useful for computation of the exponential
 *
 * Access methods for these fields are implemented.
 *
 * This class also provides the updateMatrices() method, which computes eigen values and vectors and fills the corresponding vector (eigenValues_)
 * and matrices (leftEigenVectors_ and rightEigenVectors_) from the generator.
 *
 *
 * The freq_ vector and generator_ matrices are hence the only things to provide to
 * create a substitution model.
 * It is also possible to redefine one of these methods for better efficiency.
 * The Pij_t, dPij_dt and d2Pij_dt2 are particularly inefficient since the matrix formula
 * is used to compute all probabilities, and then the result for the initial and final state
 * of interest is retrieved.
 *
 * @note This class is dedicated to "simple" substitution models, for which the number of states is equivalent to the number of characters in the alphabet.
 * Consider using the MarkovModulatedSubstitutionModel for more complexe cases.
 */

  class AbstractTransitionModel :
    public virtual TransitionModel,
    public virtual AbstractParameterAliasable
  {
  protected:
    /**
     * @brief The alphabet relevant to this model.
     */
    const Alphabet* alphabet_;

    /**
     * @brief The map of model states with alphabet states.
     */
    std::shared_ptr<const StateMap> stateMap_;

    /**
     * @brief The number of states.
     */
    size_t size_;

    /**
     * @brief The rate of the model (default: 1). The generator (and all
     * its vectorial components) is independent of the rate, since it
     * should be normalized.
     */ 
    double rate_;

    /**
     * @brief The vector \f$\pi_e\f$ of equilibrium frequencies.
     */
    Vdouble freq_;

    /**
     * @brief These ones are for bookkeeping:
     */
    mutable RowMatrix<double> pijt_;
    mutable RowMatrix<double> dpijt_;
    mutable RowMatrix<double> d2pijt_;

  public:
    AbstractTransitionModel(const Alphabet* alpha, std::shared_ptr<const StateMap> stateMap, const std::string& prefix);

    AbstractTransitionModel(const AbstractTransitionModel& model) :
      AbstractParameterAliasable(model),
      alphabet_(model.alphabet_),
      stateMap_(model.stateMap_),
      size_(model.size_),
      rate_(model.rate_),
      freq_(model.freq_),
      pijt_(model.pijt_),
      dpijt_(model.dpijt_),
      d2pijt_(model.d2pijt_)
    {}

    AbstractTransitionModel& operator=(const AbstractTransitionModel& model)
    {
      AbstractParameterAliasable::operator=(model);
      alphabet_          = model.alphabet_;
      stateMap_          = model.stateMap_;
      size_              = model.size_;
      rate_              = model.rate_;
      freq_              = model.freq_;
      pijt_              = model.pijt_;
      dpijt_             = model.dpijt_;
      d2pijt_            = model.d2pijt_;
      return *this;
    }
  
    virtual ~AbstractTransitionModel() {}

    virtual AbstractTransitionModel* clone() const = 0;

  public:
    const Alphabet* getAlphabet() const { return alphabet_; }

    const StateMap& getStateMap() const { return *stateMap_; }

    std::shared_ptr<const StateMap> shareStateMap() const { return stateMap_; }

    const std::vector<int>& getAlphabetStates() const { return stateMap_->getAlphabetStates(); }

    std::string getAlphabetStateAsChar(size_t index) const { return stateMap_->getAlphabetStateAsChar(index); }

    int getAlphabetStateAsInt(size_t index) const { return stateMap_->getAlphabetStateAsInt(index); }

    std::vector<size_t> getModelStates(int code) const { return stateMap_->getModelStates(code); }
  
    std::vector<size_t> getModelStates(const std::string& code) const { return stateMap_->getModelStates(code); }

    const Vdouble& getFrequencies() const { return freq_; }

    bool computeFrequencies() const { return false; }
    
    virtual const Matrix<double>& getPij_t(double t) const = 0;
    virtual const Matrix<double>& getdPij_dt(double t) const = 0;
    virtual const Matrix<double>& getd2Pij_dt2(double t) const = 0;

    virtual double freq(size_t i) const { return freq_[i]; }

    virtual double Pij_t    (size_t i, size_t j, double t) const { return getPij_t(t) (i, j); }
    virtual double dPij_dt  (size_t i, size_t j, double t) const { return getdPij_dt(t) (i, j); }
    virtual double d2Pij_dt2(size_t i, size_t j, double t) const { return getd2Pij_dt2(t) (i, j); }

    double getInitValue(size_t i, int state) const;

    void setFreqFromData(const SequencedValuesContainer& data, double pseudoCount = 0);

    virtual void setFreq(std::map<int, double>&);

    /**
     * @brief Tells the model that a parameter value has changed.
     *
     * This updates the matrices consequently.
     */
    virtual void fireParameterChanged(const ParameterList& parameters)
    {
      AbstractParameterAliasable::fireParameterChanged(parameters);
    
      if (parameters.hasParameter(getNamespace()+"rate"))
      {
        rate_=parameters.getParameterValue(getNamespace()+"rate");
      
        if (parameters.size()!=1)
          updateMatrices();
      }
      else
        updateMatrices();      
    }

    /**
     * @brief add a "rate" parameter to the model, that handles the
     * overall rate of the process.
     *
     */
    void addRateParameter();

  protected:
    /**
     * @brief Diagonalize the \f$Q\f$ matrix, and fill the eigenValues_, iEigenValues_, 
     * leftEigenVectors_ and rightEigenVectors_ matrices.
     *
     * The generator_ matrix and freq_ vector must be initialized.
     *
     * Eigen values and vectors are computed from the generator and
     * assigned to the eigenValues_ for the real part, iEigenValues_ for
     * the imaginary part, rightEigenVectors_ and leftEigenVectors_
     * variables. isDiagonalizable_ checks if the generator_ is
     * diagonalizable in R.
     *
     * The optional rate parameter is not taken into account in this
     * method to prevent unnecessary computation.
     *
     * !! Here there is no normalization of the generator.
     * 
     */
    virtual void updateMatrices() = 0;

    /*
     * @brief : To update the eq freq
     *
     */
    
    Vdouble& getFrequencies_() {
      return freq_;
    }

  public:

    /**
     * @brief The rate of the substitution process.
     *
     */
    virtual double getRate() const;

    virtual void setRate(double rate);

  };


  class AbstractSubstitutionModel :
    public AbstractTransitionModel,
    virtual public SubstitutionModel
  {
  protected:
    /**
     * @brief If the model is scalable (ie generator can be normalized
     * automatically).
     */
  
    bool isScalable_;

    /**
     * @brief The generator matrix \f$Q\f$ of the model.
     */
    RowMatrix<double> generator_;

    /**
     * @brief if the Frequencies must be computed from the generator
     */
  
    bool computeFreq_;

    /**
     * @brief The exchangeability matrix \f$S\f$ of the model, defined
     * as \f$ S_{ij}=\frac{Q_{ij}}{\pi_j}\f$. When the model is
     * reversible, this matrix is symetric.
     *
     */
    RowMatrix<double> exchangeability_;

    /**
     * @brief Tell if the eigen decomposition should be performed.
     */
    bool eigenDecompose_;

    /**
     * @brief The vector of eigen values.
     */
    Vdouble eigenValues_;

    /**
     * @brief The vector of the imaginary part of the eigen values.
     */
    Vdouble iEigenValues_;

    /**
     * @brief boolean value for diagonalizability in R of the generator_
     */
    bool isDiagonalizable_;

    /**
     * @brief The \f$U^-1\f$ matrix made of right eigen vectors (by column).
     */
    RowMatrix<double> rightEigenVectors_;

    /**
     * @brief boolean value for non-singularity of rightEigenVectors_
     */
    bool isNonSingular_;

    /**
     * @brief The \f$U\f$ matrix made of left eigen vectors (by row) if
     * rightEigenVectors_ is non-singular.
     */
    RowMatrix<double> leftEigenVectors_;

    /**
     * @brief vector of the powers of generator_ for Taylor development (if
     * rightEigenVectors_ is singular).
     */
    std::vector< RowMatrix<double> > vPowGen_;

    /**
     * @brief For computational issues
     */
    mutable RowMatrix<double> tmpMat_;
  
  public:
    AbstractSubstitutionModel(const Alphabet* alpha, std::shared_ptr<const StateMap> stateMap, const std::string& prefix);

    AbstractSubstitutionModel(const AbstractSubstitutionModel& model) :
      AbstractParameterAliasable(model),
      AbstractTransitionModel(model),
      isScalable_(model.isScalable_),
      generator_(model.generator_),
      computeFreq_(model.computeFreq_),
      exchangeability_(model.exchangeability_),
      eigenDecompose_(model.eigenDecompose_),
      eigenValues_(model.eigenValues_),
      iEigenValues_(model.iEigenValues_),
      isDiagonalizable_(model.isDiagonalizable_),
      rightEigenVectors_(model.rightEigenVectors_),
      isNonSingular_(model.isNonSingular_),
      leftEigenVectors_(model.leftEigenVectors_),
      vPowGen_(model.vPowGen_),
      tmpMat_(model.tmpMat_)
    {}

    AbstractSubstitutionModel& operator=(const AbstractSubstitutionModel& model)
    {
      AbstractTransitionModel::operator=(model);
      isScalable_        = model.isScalable_;
      generator_         = model.generator_;
      computeFreq_       = model.computeFreq_;
      exchangeability_   = model.exchangeability_;
      eigenDecompose_    = model.eigenDecompose_;
      eigenValues_       = model.eigenValues_;
      iEigenValues_      = model.iEigenValues_;
      isDiagonalizable_  = model.isDiagonalizable_;
      rightEigenVectors_ = model.rightEigenVectors_;
      isNonSingular_     = model.isNonSingular_;
      leftEigenVectors_  = model.leftEigenVectors_;
      vPowGen_           = model.vPowGen_;
      tmpMat_            = model.tmpMat_;
      return *this;
    }
  
    virtual ~AbstractSubstitutionModel() {}

    virtual AbstractSubstitutionModel* clone() const = 0;

  public:
    bool computeFrequencies() const { return computeFreq_; }

    void computeFrequencies(bool yn) { computeFreq_=yn; }

    const Matrix<double>& getGenerator() const { return generator_; }

    const Matrix<double>& getExchangeabilityMatrix() const { return exchangeability_; }

    const Matrix<double>& getPij_t(double t) const;
    const Matrix<double>& getdPij_dt(double t) const;
    const Matrix<double>& getd2Pij_dt2(double t) const;

    double Sij(size_t i, size_t j) const { return exchangeability_(i, j); }

    const Vdouble& getEigenValues() const { return eigenValues_; }

    const Vdouble& getIEigenValues() const { return iEigenValues_; }

    bool isDiagonalizable() const { return isDiagonalizable_; }
  
    bool isNonSingular() const { return isNonSingular_; }

    const Matrix<double>& getRowLeftEigenVectors() const { return leftEigenVectors_; }

    const Matrix<double>& getColumnRightEigenVectors() const { return rightEigenVectors_; }

    virtual double Qij(size_t i, size_t j) const { return generator_(i, j); }

    void enableEigenDecomposition(bool yn) { eigenDecompose_ = yn; }

    bool enableEigenDecomposition() { return eigenDecompose_; }

  protected:
    /**
     * @brief Diagonalize the \f$Q\f$ matrix, and fill the eigenValues_, iEigenValues_, 
     * leftEigenVectors_ and rightEigenVectors_ matrices.
     *
     * The generator_ matrix and freq_ vector must be initialized.
     *
     * Eigen values and vectors are computed from the generator and
     * assigned to the eigenValues_ for the real part, iEigenValues_ for
     * the imaginary part, rightEigenVectors_ and leftEigenVectors_
     * variables. isDiagonalizable_ checks if the generator_ is
     * diagonalizable in R.
     *
     * The optional rate parameter is not taken into account in this
     * method to prevent unnecessary computation.
     *
     * !! Here there is no normalization of the generator.
     * 
     */
    virtual void updateMatrices();

  public:

    /**
     * @brief sets if model is scalable, ie scale can be changed.
     * Default : true, set to false to avoid normalization for example.
     *
     */
  
    void setScalable(bool scalable)
    {
      isScalable_=scalable;
    }
  
    /**
     * @brief returns  if model is scalable
     *
     */

    virtual bool isScalable() const
    {
      return isScalable_;
    }

    /**
     * @brief return scale
     *
     */
  
    double getScale() const;

    /**
     * @brief Multiplies the current generator by the given scale.
     *
     * @param scale the scale by which the generator is multiplied.
     *
     */
    void setScale(double scale);

    /**
     * @brief normalize the generator
     *
     */

    void normalize();

    /**
     * @brief set the diagonal of the generator such that sum on each
     * line equals 0.
     *
     */
  
    void setDiagonal();

  };


/**
 * @brief Partial implementation of the ReversibleSubstitutionModel interface.
 *
 * This class overrides the updateMatrices() method, which updates the
 * generator_ matrix from the exchangeability_ matrix and freq_
 * vector. It then computes eigen values and vectors and fills the
 * corresponding vector (eigenValues_) and matrices (leftEigenVectors_
 * and rightEigenVectors_). Because of reversibility,
 * isDiagonalizable_ is set to true.
 *
 * The freq_ vector and exchangeability_ matrices are hence the only
 * things to provide to create a substitution model. It is also
 * possible to redefine one of these methods for better efficiency.
 * The Pij_t, dPij_dt and d2Pij_dt2 are particularly inefficient since
 * the matrix formula is used to compute all probabilities, and then
 * the result for the initial and final state of interest is
 * retrieved.
 *
 * @note This class is dedicated to "simple" substitution models, for
 * which the number of states is equivalent to the number of
 * characters in the alphabet. Consider using the
 * MarkovModulatedSubstitutionModel for more complexe cases.
 */
  class AbstractReversibleSubstitutionModel :
    public AbstractSubstitutionModel,
    public virtual ReversibleSubstitutionModel
  {
  public:
    AbstractReversibleSubstitutionModel(const Alphabet* alpha, std::shared_ptr<const StateMap> stateMap, const std::string& prefix) :
      AbstractParameterAliasable(prefix),
      AbstractSubstitutionModel(alpha, stateMap, prefix)
    {
      isDiagonalizable_ = true;
      isNonSingular_    = true;
      computeFreq_ = false;
    }

    virtual ~AbstractReversibleSubstitutionModel() {}

    virtual AbstractReversibleSubstitutionModel* clone() const = 0;

  protected:

    /**
     * @brief Compute and diagonalize the \f$Q\f$ matrix, and fill the eigenValues_,
     * leftEigenVectors_ and rightEigenVectors_ matrices.
     *
     * The exchangeability_ matrix and freq_ vector must be initialized.
     * This function computes the generator_ matrix with the formula
     * \f[
     * Q = S \times \pi
     * \f]
     * where \f$Q\f$ is the generator matrix, \f$S\f$ is the exchangeability matrix and
     * \f$Pi\f$ the diagonal matrix with frequencies.
     *
     * The generator is then scaled so that
     * \f[
     * \sum_i Q_{i,i} \times \pi_i = -1
     * \f]
     * (\f$\pi_i\f$ are the equilibrium frequencies).
     *
     * Eigen values and vectors are computed from the scaled generator and assigned to the
     * eigenValues_, rightEigenVectors_ and leftEigenVectors_ variables.
     */
    virtual void updateMatrices();

  };

} //end of namespace bpp.

#endif  //_ABSTRACTSUBSTITUTIONMODEL_H_

