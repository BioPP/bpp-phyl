//
// File: SubstitutionModel.h
// Created by:  <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon May 26 14:52:34 2003
//

#ifndef _SUBSTITUTIONMODEL_H_
#define _SUBSTITUTIONMODEL_H_

#include <cstdlib>
#include <map>
#include <string>

//From Seqlib:
#include <Seq/Alphabet.h>

// From utils:
#include <Utils/Exceptions.h>

// From NumCalc:
#include <NumCalc/Parameter.h>
#include <NumCalc/ParameterList.h>
#include <NumCalc/Parametrizable.h>
#include <NumCalc/VectorTools.h>

// From the MTL:
#include <mtl/matrix.h>

using namespace mtl;

typedef matrix<double>::type Matrix;
typedef Vdouble Vector;

class SubstitutionModel;

/******************************************************************************
 *                        SubstitutionModel exceptions:                       *
 ******************************************************************************/

class SubstitutionModelException : public Exception {

	protected:
		const SubstitutionModel * sm;
			
	public:
		// Class constructor
		SubstitutionModelException(const char *   text, const SubstitutionModel * sm = NULL);
		SubstitutionModelException(const string & text, const SubstitutionModel * sm = NULL);
	
		// Class destructor
		~SubstitutionModelException() throw ();
	public:
		virtual const SubstitutionModel * getSubstitutionModel() const;
};

/*
 * This class is a base class for Markov models.
 */

class SubstitutionModel: public Parametrizable {
	
	public:
		//Destructor:
		virtual ~SubstitutionModel() {};

	public:
		virtual double freq     (int i)                  const = 0;
		virtual double Qij      (int i, int j)           const = 0;
		virtual double Pij_t    (int i, int j, double t) const = 0;
		virtual double dPij_dt  (int i, int j, double t) const = 0;
		virtual double d2Pij_dt2(int i, int j, double t) const = 0;
	
		virtual Vector freq() const = 0;
	
		virtual Matrix getGenerator() const = 0;

		virtual Matrix getPij_t(double t) const = 0;
	
		virtual Matrix getdPij_dt(double t) const = 0;

		virtual Matrix getd2Pij_dt2(double t) const = 0;

		virtual Vector eigenValues() const = 0;

		virtual Matrix leftEigenVector() const = 0;

		virtual Matrix rightEigenVectors() const = 0;

		virtual const Alphabet * getAlphabet() const = 0;
		/* This method is used to initialize likelihoods in reccurence.
		 * Traditionaly, it sends 1 if i = state, 0 else, where
		 * i is one of the possible states if the alphabet allowed in the model
		 * and state is the observed state in the sequence.
		 */
		virtual double getInitValue(int i, int state) const throw (BadIntException) = 0;

		/* Get the name of the model:
		 */
		virtual string getName() const = 0;
	
};

#endif	//_SUBSTITUTIONMODEL_H_
