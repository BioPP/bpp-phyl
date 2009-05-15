//
// File: MarkovModulatedSubstitutionModel.h
// Created by: Julien Dutheil
// Created on: Sat Aug 05 08:21 2006
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _MARKOVMODULATEDSUBSTITUTIONMODEL_H_
#define _MARKOVMODULATEDSUBSTITUTIONMODEL_H_

#include "SubstitutionModel.h"

//From NumCalc:
#include <NumCalc/AbstractParameterAliasable.h>

namespace bpp
{

/**
 * @brief Partial implementation of the Markov-modulated class of substitution
 * models.
 *
 * This class wraps a substitution model and provide a matrix describing rate changes.
 * The rate matrix must be initialized by derived classes of this class.
 * Using these matrices, the diagonalization procedure of Galtier and Jean-Marie is used.
 *
 * Such models can be described using two matrices:
 * a substitution matrix, @f$M@f$, with size @f$m@f$, which is a "standard" substitution model of any alphabet type,
 * and a rate matrix @f$G@f$ of size @f$g@f$.
 * The generator of the markov-modulated model, @f$Q@f$ can be written using Kronecker matrix operands:
 * @f[
 * Q=D_R \otimes M + G \otimes I_m,
 * @f]
 * where @f$D_R@f$ is the diagonal matrix of all rates, and @f$I_m@f$ is the identity matrix of size @f$m@f$.
 *
 * This generator is normalized so that branch lengths are measured in unit of mean number of substitutions per site,
 * where susbstitution here means "change of alphabet state".
 * Rate changes are not counted.
 *
 * Galtier N. and Jean-Marie A., Markov-modulated Markov chains and the covarion process of molecular evolution (2004).
 * _Journal of Computational Biology_, 11:727-33.
 */
class MarkovModulatedSubstitutionModel:
  public ReversibleSubstitutionModel,
  public AbstractParameterAliasable
{

  protected:
    ReversibleSubstitutionModel* _model;
    unsigned int _nbStates; //Number of states in model
    unsigned int _nbRates; //Number of rate classes
    /**
     * @name Rate generator.
     *
     * These variables must be initialized in the constructor of the derived class.
     * @{
     */
    RowMatrix<double> _rates;                //All rates values
    RowMatrix<double> _ratesExchangeability; //All rates transitions
    Vdouble           _ratesFreq;            //All rates equilibrium frequencies
    /**@}*/
    RowMatrix<double> _ratesGenerator;       //All rates transitions
    
		/**
		 * @brief The generator matrix \f$Q\f$ of the model.
		 */
		RowMatrix<double> _generator;

		/**
		 * @brief The exchangeability matrix \f$S\f$ of the model.
		 */
		RowMatrix<double> _exchangeability;

		/**
		 * @brief The \f$U\f$ matrix made of left eigen vectors (by row).
		 */
		RowMatrix<double> _leftEigenVectors;

		/**
		 * @brief The \f$U^-1\f$ matrix made of right eigen vectors (by column).
		 */
		RowMatrix<double> _rightEigenVectors;

    /**
     * @brief These ones are for bookkeeping:
     */
    mutable RowMatrix<double> _pijt;
    mutable RowMatrix<double> _dpijt;
    mutable RowMatrix<double> _d2pijt;

		/**
		 * @brief The vector of eigen values.
		 */
		Vdouble _eigenValues;

		/**
		 * @brief The vector of equilibrium frequencies.
		 */
		Vdouble _freq;

    bool _normalizeRateChanges;

    string nestedPrefix_;

  public:
    /**
     * @brief Build a new MarkovModulatedSubstitutionModel object.
     *
     * @param model The substitution model to use. Can be of any alphabet type, and will be owned by this instance.
     * @param normalizeRateChanges Tells if the branch lengths must be computed in terms of rate and state
     * @param prefix The parameter namespace to be forwarded to the AbstractParametrizable constructor.
     * changes instead of state change only.
     * NB: In most cases, this parameter should be set to false.
     */
    MarkovModulatedSubstitutionModel(ReversibleSubstitutionModel* model, bool normalizeRateChanges, const string& prefix) :
      AbstractParameterAliasable(prefix),
      _model(model), _nbStates(0), _nbRates(0), _rates(), _ratesExchangeability(),
      _ratesFreq(), _ratesGenerator(), _generator(), _exchangeability(),
      _leftEigenVectors(), _rightEigenVectors(), _eigenValues(), _freq(), _normalizeRateChanges(normalizeRateChanges),
      nestedPrefix_(model->getNamespace())
    {
      _model->setNamespace(prefix + nestedPrefix_);
      addParameters_(_model->getIndependentParameters());
    }
    
    MarkovModulatedSubstitutionModel(const MarkovModulatedSubstitutionModel & model);
    MarkovModulatedSubstitutionModel & operator=(const MarkovModulatedSubstitutionModel & model);

    virtual ~MarkovModulatedSubstitutionModel() { delete _model; }

#ifndef NO_VIRTUAL_COV
    MarkovModulatedSubstitutionModel*
#else
    Clonable*
#endif
    clone() const = 0;

  public:
	  
		const Alphabet * getAlphabet() const { return _model->getAlphabet(); }

    unsigned int getNumberOfStates() const { return _nbStates*_nbRates; }

		const Vdouble & getFrequencies() const { return _freq; }
    
		const Matrix<double> & getExchangeabilityMatrix() const { return _exchangeability; }
    
		const Matrix<double> & getGenerator() const { return _generator; }
    
		const Matrix<double> & getPij_t(double t) const;
		const Matrix<double> & getdPij_dt(double t) const;
		const Matrix<double> & getd2Pij_dt2(double t) const;
    
		const Vdouble & getEigenValues() const { return _eigenValues; }
    
		const Matrix<double> & getRowLeftEigenVectors() const { return _leftEigenVectors; }
		const Matrix<double> & getColumnRightEigenVectors() const { return _rightEigenVectors; }
    
		double freq(int i) const { return _freq[i]; }
    double Sij(int i, int j) const { return _exchangeability(i, j); }
		double Qij(int i, int j) const { return _generator(i, j); }
    
		double Pij_t    (int i, int j, double t) const { return getPij_t(t)(i, j); }
		double dPij_dt  (int i, int j, double t) const { return getdPij_dt(t)(i, j); }
		double d2Pij_dt2(int i, int j, double t) const { return getd2Pij_dt2(t)(i, j); }
    
		double getInitValue(int i, int state) const throw (BadIntException);
    
		void setFreqFromData(const SequenceContainer & data, unsigned int pseudoCount = 0)
    {
      _model->setFreqFromData(data, pseudoCount);
      updateMatrices();
    }

    int getState(int i) const
    {
      return i % _nbStates; 
    }
    
    /**
     * @brief Get the rate category corresponding to a particular state in the compound model.
     *
     * @param i The state.
     * @return The corresponding rate category.
     * @see getState;
     */
    int getRate(int i) const
    {
      return i / _nbStates; 
    }

    /**
		 * @brief Tells the model that a parameter value has changed.
		 *
		 * This updates the matrices consequently.
		 */
		virtual void fireParameterChanged(const ParameterList & parameters)
    {
      _model->matchParametersValues(parameters);
      updateRatesModel();
      updateMatrices();
    }
   
    void setNamespace(const string& prefix);

  protected:
    
    virtual void updateMatrices();

    /**
     * @brief Update the rates vector, generator and equilibrium frequencies.
     *
     * This method must be implemented by the derived class.
     * It is called by the fireParameterChanged() method.
     */
    virtual void updateRatesModel() = 0;

};

} //end of namespace bpp.

#endif //_MARKOVMODULATEDSUBSTITUTIONMODEL_H_

