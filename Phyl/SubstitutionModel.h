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
#include <Seq/SequenceContainer.h>

// From utils:
#include <Utils/Exceptions.h>

// From NumCalc:
#include <NumCalc/Parameter.h>
#include <NumCalc/ParameterList.h>
#include <NumCalc/Parametrizable.h>
#include <NumCalc/VectorTools.h>
#include <NumCalc/Matrix.h>

// From the MTL:
//#include <mtl/matrix.h>
//using namespace mtl;
//typedef matrix<double>::type Matrix;

typedef RowMatrix<double> Mat;
typedef Vdouble Vec;

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

/******************************************************************************
 *                        SubstitutionModel class:                            *
 ******************************************************************************/

/**
 * @brief Interface for all substitution models.
 * 
 * A substitution model is based on a Markov generator \f$Q\f$, the size of
 * which depends on the alphabet used (4 for nucleotides, 20 for proteins, etc.).
 * Each SubstitutionModel will hence include a pointer toward an alphabet, and
 * provides a method to know which alphabet is used (getAlphabet() method).
 *
 * What we want from a substitution model is to compute the probabilities that
 * one particular state \f$i\f$ mutate into state \f$j\f$ after a time \f$t\f$
 * (\f$P_{i,j}(t)\f$).
 * Typically, this is computed using the formula
 * \f[
 * P(t) = e^{t \times Q},
 * \f]
 * where \f$P(t)\f$ is the matrix with all probabilities \f$P_{i,j}(t)\f$.
 * For some models, such \f$P_{i,j}(t)\f$ may be computed analytically.
 * For more complexe models, we need to use a eigen-decomposition of \f$Q\f$:
 * \f[ Q = U \times D \times U^{-1}, \f]
 * where \f$D\f$ is a diagonal matrix.
 * Hence
 * \f[ P(t) = e^{t \times Q} = U \times e^{D} \times U^{-1}, \f]
 * where \f$e^{D}\f$ is a diagonal matrix with all terms equal to exp the terms
 * in \f$D\f$.
 * \f$U\f$ is the matrix of vertical right eigen vectors, and \f$U^{-1}\f$ is the matrix
 * of vertical left eigen vectors.
 * The values on the diagonal of \f$D\f$ are the eigen values of \f$Q\f$.
 * All \f$Q,U,U^{-1}\f$ and \f$D\f$ (its diagonal) may be retrieved from the
 * class.
 *
 * Moreover, for reversible models, we can write:
 * \f[ Q = S \times \pi, \f]
 * where \f$S\f$ is a symetric matrix and \f$\pi\f$ the diagonal matrix with
 * all equilibrium frequencies.
 * The frequences may be retrieved as a vector by the getFrequencies() method
 * or individually by the freq() method.
 * The \f$S\f$ matrix may be obtained by the getExchangeabilityMatrix().
 *
 * Moreover, the equilibrium frequencies may also be retrieved, and first and
 * second order derivatives of \f$P(t)\f$ according to \f$t\f$.
 * These methods may be useful for optimization processes.
 */

class SubstitutionModel: public virtual Parametrizable {
	
	public:
		//Destructor:
		virtual ~SubstitutionModel() {};

	public:
		
		/**
		 * @brief Get the name of the model.
		 *
		 * @return The name of this model.
		 */
		virtual string getName() const = 0;
	
		virtual double freq(int i) const = 0;

		virtual double Qij(int i, int j) const = 0;

		virtual double Pij_t(int i, int j, double t) const = 0;

		virtual double dPij_dt(int i, int j, double t) const = 0;
		
		virtual double d2Pij_dt2(int i, int j, double t) const = 0;
	
		virtual Vec getFrequencies() const = 0;
		
		virtual Mat getExchangeabilityMatrix() const = 0;

		virtual Mat getGenerator() const = 0;

		virtual Mat getPij_t(double t) const = 0;
	
		virtual Mat getdPij_dt(double t) const = 0;

		virtual Mat getd2Pij_dt2(double t) const = 0;

		virtual Vec eigenValues() const = 0;

		virtual Mat horizontalLeftEigenVectors() const = 0;

		virtual Mat verticalRightEigenVectors() const = 0;

		virtual const Alphabet * getAlphabet() const = 0;

		/**
		 * This method is used to initialize likelihoods in reccurence.
		 * Traditionaly, it sends 1 if i = state, 0 else, where
		 * i is one of the possible states of the alphabet allowed in the model
		 * and state is the observed state in the sequence.
		 *
		 * @param i one of the possible states of the alphabet.
		 * @param state An observed state in the sequence.
		 * @return 1 or 0 dpeendeing if two states are compatible.
		 * @throw BadIntException if states are not allowed in the associated alphabet.
		 */
		virtual double getInitValue(int i, int state) const throw (BadIntException) = 0;

		/**
		 * @brief Set equilibrium frequencies equal to the frequencies estimated
		 * from the data.
		 */
		virtual void setFreqFromData(const SequenceContainer & data) = 0;

};

#endif	//_SUBSTITUTIONMODEL_H_
