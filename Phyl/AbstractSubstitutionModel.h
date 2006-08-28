//
// File: AbstractSubstitutionModel.h
// Created by: Julien Dutheil
// Created on: Tue May 27 10:31:49 2003
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

#ifndef _ABSTRACTSUBSTITUTIONMODEL_H_
#define _ABSTRACTSUBSTITUTIONMODEL_H_

#include "SubstitutionModel.h"

/**
 * @brief Partial implementation of the SubstitutionModel interface.
 *
 * This abstract class provides some fields, namely:
 * - _alphabet: a pointer toward the alphabet,
 * - _size: the size of the alphabet, a parameter frequenctly called during various computations,
 * - _generator, _exchangeability, _leftEigenVectors, _rightEigenVectors: usefull matrices,
 * - _eigenValues, _freq: usefull vectors.
 *
 * Access methods for these fields are implemented.
 *
 * This class also provides the updateMatrices() method, which updates
 * the _generator matrix from the _exchangeability matrix and _freq vector.
 * It then computes eigen values and vectors and fills the corresponding vector (_eigenValues)
 * and matrices (_leftEigenVectors and _rightEigenVectors).
 *
 * The _freq vector and _exchangeability matrices are hence the only things to provide to
 * create a substitution model.
 * It is also possible to redefine one of these methods for better efficiency.
 * The Pij_t, dPij_dt and d2Pij_dt2 are particularly inefficient since the matrix formula
 * is used to compute all probabilities, and then the result for the initial and final state
 * of interest is retrieved.
 *
 * @note This class is dedicated to "simple" substitution models, for which the number of states is equivalent to the number of characters in the alphabet.
 * Consider using the MarkovModulatedSubstitutionModel for more complexe cases.
 */
class AbstractSubstitutionModel :
	public virtual SubstitutionModel,
	public virtual AbstractParametrizable
{
	protected:

		/**
		 * @brief The alphabet relevant to this model.
		 */
		const Alphabet * _alphabet;

		/**
		 * @brief The size of the generator, i.e. the number of states.
		 */
		unsigned int _size;

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
		 * @brief The vector of eigen values.
		 */
		Vdouble _eigenValues;

		/**
		 * @brief The vector of equilibrium frequencies.
		 */
		Vdouble _freq;

	public:
		AbstractSubstitutionModel(const Alphabet * alpha);
	
		virtual ~AbstractSubstitutionModel() {}
	
	public:
		const Alphabet * getAlphabet() const { return _alphabet; }
    
		Vdouble getFrequencies() const { return _freq; }
    
		RowMatrix<double> getExchangeabilityMatrix() const { return _exchangeability; }
    
		RowMatrix<double> getGenerator() const { return _generator; }
    
		RowMatrix<double> getPij_t(double t) const;
		RowMatrix<double> getdPij_dt(double t) const;
		RowMatrix<double> getd2Pij_dt2(double t) const;
    
		Vdouble getEigenValues() const { return _eigenValues; }
    
		RowMatrix<double> getRowLeftEigenVectors() const { return _leftEigenVectors; }
		RowMatrix<double> getColumnRightEigenVectors() const { return _rightEigenVectors; }
    
		double freq(int i) const { return _freq[i]; }
    
		double Qij(int i, int j) const { return _generator(i, j); }
    
		double Pij_t    (int i, int j, double t) const { return getPij_t(t)(i, j); }
		double dPij_dt  (int i, int j, double t) const { return getdPij_dt(t)(i, j); }
		double d2Pij_dt2(int i, int j, double t) const { return getd2Pij_dt2(t)(i, j); }

		double getInitValue(int i, int state) const throw (BadIntException);

		void setFreqFromData(const SequenceContainer & data);

		/**
		 * @brief Tells the model that a parameter value has changed.
		 *
		 * This updates the matrices consequently.
		 */
		void fireParameterChanged(const ParameterList & parameters)
		{
			updateMatrices();
		}

		
	protected:

		/**
		 * @brief Compute and diagonalize the \f$Q\f$ matrix, and fill the _eigenValues,
		 * _leftEigenVectors and _rightEigenVectors matrices.
		 *
		 * The _exchangeability matrix and _freq vector must be initialized.
		 * This function computes the _generator matrix with the formula
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
		 * WARNING!!! The exchangeability matrix is not scaled, only the generator matrix.
		 * 
		 * Eigen values and vectors are computed from the scaled generator and assigned to the
		 * _eigenValues, _rightEigenVectors and _leftEigenVectors variables.
		 */
		virtual void updateMatrices();

		/**
		 * @brief Get the scalar product of diagonal elements of the generator
		 * and the frequencies vector.
		 * If the generator is normalized, then scale=1. Otherwise each element
		 * must be multiplied by 1/scale.
		 *
		 * @return Minus the scalar product of diagonal elements and the frequencies vector.
		 */
		double getScale() const;
};

#endif	//_ABSTRACTSUBSTITUTIONMODEL_H_

