//
// File: AbstractSubstitutionModel.h
// Created by:  <@bogdanof>
// Created on: Tue May 27 10:31:49 2003
//

#ifndef _ABSTRACTSUBSTITUTIONMODEL_H_
#define _ABSTRACTSUBSTITUTIONMODEL_H_

#include "SubstitutionModel.h"

class AbstractSubstitutionModel : public SubstitutionModel
{
	protected:
		const Alphabet * alphabet;
		unsigned int _size;
		Matrix _generator;
		Matrix _leftEigenVectors, _rightEigenVectors;
		Vector _eigenValues;
		Vector _freq;
		mutable ParameterList _parameters;

	public:
		AbstractSubstitutionModel(const Alphabet * alpha);
	
		virtual ~AbstractSubstitutionModel() {}
	
	protected:
		virtual void fillMatrices() = 0;
			
	public:
		const Alphabet * getAlphabet() const;
		double getInitValue(int i, int state) const throw (BadIntException);

		Vector freq() const;
		Matrix getGenerator() const;
		Vector eigenValues() const;
		Matrix leftEigenVector() const;
		Matrix rightEigenVectors() const;
		double freq(int i) const;
		double Qij(int i, int j) const;
	
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
		
};


#endif	//_ABSTRACTSUBSTITUTIONMODEL_H_
