//
// File: AbstractSubstitutionModel.h
// Created by:  <@bogdanof>
// Created on: Tue May 27 10:31:49 2003
//

#ifndef _ABSTRACTSUBSTITUTIONMODEL_H_
#define _ABSTRACTSUBSTITUTIONMODEL_H_

#include "SubstitutionModel.h"

/**
 * @brief Low level implementation of the SubstitutionModel interface.
 */
class AbstractSubstitutionModel : public SubstitutionModel
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
		Matrix _generator;

		/**
		 * @brief The exchangeability matrix \f$S\f$ of the model.
		 */
		Matrix _exchangeability;

		/**
		 * @brief The \f$U^{-1}\f$ matrix made of horizontal left eigen vectors.
		 */
		Matrix _leftEigenVectors;

		/**
		 * @brief The \f$U\f$ matrix made of vertical right eigen vectors.
		 */
		Matrix _rightEigenVectors;

		/**
		 * @brief The vector of eigen values.
		 */
		Vector _eigenValues;

		/**
		 * @brief The vector of equilibrium frequencies.
		 */
		Vector _freq;

		/**
		 * @brief The parameters of this model (needed for the Parametrizable interface).
		 */
		mutable ParameterList _parameters;

	public:
		AbstractSubstitutionModel(const Alphabet * alpha);
	
		virtual ~AbstractSubstitutionModel() {}
	
	public:
		const Alphabet * getAlphabet() const;

		Vector getFrequencies() const;
		Matrix getExchangeabilityMatrix() const;
		Matrix getGenerator() const;
		Matrix getPij_t(double t) const;
		Matrix getdPij_dt(double t) const;
		Matrix getd2Pij_dt2(double t) const;
		Vector eigenValues() const;
		Matrix horizontalLeftEigenVectors() const;
		Matrix verticalRightEigenVectors() const;
		double freq(int i) const;
		double Qij(int i, int j) const;
		double Pij_t    (int i, int j, double t) const;
		double dPij_dt  (int i, int j, double t) const;
		double d2Pij_dt2(int i, int j, double t) const;

		double getInitValue(int i, int state) const throw (BadIntException);
		void setFreqFromData(const SequenceContainer & data);


		/**
		 * @name The Parametrizable interface.
		 *
		 * @{
		 */
		ParameterList getParameters() const ;
	
		double getParameter (const string & name) const
			 throw (ParameterNotFoundException);
			 
		void setAllParametersValues(const ParameterList & params)
			 throw (ParameterNotFoundException, ConstraintException);
			 
		void setParameterValue(const string & name, double value)
			 throw (ParameterNotFoundException, ConstraintException);
			 
		void setParametersValues(const ParameterList & params)
			 throw (ParameterNotFoundException, ConstraintException);
			 
		void matchParametersValues(const ParameterList & params)
			 throw (ConstraintException);
			 
		/** @} */

	protected:
		/**
		 * @brief Compute and diagonalize the $\$Q\f$ matrix and fill the _eigenValues,
		 * _leftEigenVectors and _rightEigenVectors fields.
		 *
		 * This routine uses the MTL interface to Lapack to compute eigen values
		 * and vectors. However, the inverse of the right eigen vectors is used
		 * as left eigen vectors, since the left eigen vectors sent by the lapack
		 * geev routine is not exactly the inverse of the right eigen vectors.
		 * (Dunno why though...)
		 */
		virtual void updateMatrices() = 0;

};


#endif	//_ABSTRACTSUBSTITUTIONMODEL_H_
