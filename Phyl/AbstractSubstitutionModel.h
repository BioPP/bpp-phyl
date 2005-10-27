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
 * @brief Low level implementation of the SubstitutionModel interface.
 */
class AbstractSubstitutionModel :
	public virtual SubstitutionModel,
	public virtual AbstractParametrizable
{
	protected:

		/**
		 * @brief The alphabet this model is to use with.
		 */
		const Alphabet * alphabet;

		/**
		 * @brief The size of the alphabet.
		 */
		unsigned int _size;

		/**
		 * @brief The generator matrix \f$Q\f$ of the model.
		 */
		Mat _generator;

		/**
		 * @brief The exchangeability matrix \f$S\f$ of the model.
		 */
		Mat _exchangeability;

		/**
		 * @brief The \f$U^{-1}\f$ matrix made of horizontal left eigen vectors.
		 */
		Mat _leftEigenVectors;

		/**
		 * @brief The \f$U\f$ matrix made of vertical right eigen vectors.
		 */
		Mat _rightEigenVectors;

		/**
		 * @brief The vector of eigen values.
		 */
		Vec _eigenValues;

		/**
		 * @brief The vector of equilibrium frequencies.
		 */
		Vec _freq;

	public:
		AbstractSubstitutionModel(const Alphabet * alpha);
	
		virtual ~AbstractSubstitutionModel() {}
	
	public:
		const Alphabet * getAlphabet() const;

		Vec getFrequencies() const;
		Mat getExchangeabilityMatrix() const;
		Mat getGenerator() const;
		Mat getPij_t(double t) const;
		Mat getdPij_dt(double t) const;
		Mat getd2Pij_dt2(double t) const;
		Vec eigenValues() const;
		Mat verticalLeftEigenVectors() const;
		Mat horizontalRightEigenVectors() const;
		double freq(int i) const;
		double Qij(int i, int j) const;
		double Pij_t    (int i, int j, double t) const;
		double dPij_dt  (int i, int j, double t) const;
		double d2Pij_dt2(int i, int j, double t) const;

		double getInitValue(int i, int state) const throw (BadIntException);
		void setFreqFromData(const SequenceContainer & data);

		void fireParameterChanged(const ParameterList & parameters)
		{
			updateMatrices();
		}
		
	protected:

		/**
		 * @brief Compute and diagonalize the \f$Q\f$ matrix and fill the _eigenValues,
		 * _leftEigenVectors and _rightEigenVectors fields.
		 *
		 * The _exchangeability matrix and _freq vector must be initialized.
		 * This function computes the _generator matrix with the formula
		 * \f[
		 * Q = S \times \Pi
		 * \f]
		 * where \f$Q\f$ is the genrator matrix, \f$S\f$ is the exchangeability matrix and
		 * \f$Pi\f$ the diagonal matrix with frequencies.
		 *
		 * The generator is then scaled so that
		 * \f[
		 * \sum_i Q_{i,i} \times f_i = -1
		 * \f]
		 * (\f$f_i\f$ are the equilibrium frequencies).
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

