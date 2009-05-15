//
// File: SubstitutionModel.h
// Created by: Julien Dutheil
// Created on: Mon May 26 14:52:34 2003
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
#include <NumCalc/ParameterAliasable.h>
#include <NumCalc/VectorTools.h>
#include <NumCalc/Matrix.h>

namespace bpp
{

class SubstitutionModel;

/**
 * @brief Exception that may be thrown by susbstitution models.
 *
 * @see SubstitutionModel
 */
class SubstitutionModelException:
  public Exception
{
	protected:
		const SubstitutionModel * sm;
			
	public:
		SubstitutionModelException(const char *   text, const SubstitutionModel * sm = NULL);
		SubstitutionModelException(const string & text, const SubstitutionModel * sm = NULL);
	
		~SubstitutionModelException() throw ();
		
	public:
		/**
		 * @brief Get the model that throw the exception.
		 *
		 * @return The model that throw the exception.
		 */
		virtual const SubstitutionModel * getSubstitutionModel() const;
};

/**
 * @brief Interface for all substitution models.
 * 
 * A substitution model is based on a Markov generator \f$Q\f$, the size of
 * which depends on the alphabet used (4 for nucleotides, 20 for proteins, etc.).
 * Each SubstitutionModel object hence includes a pointer toward an alphabet,
 * and provides a method to retrieve the alphabet used (getAlphabet() method).
 *
 * What we want from a substitution model is to compute the probabilities of state
 * j at time t geven state j at time 0 (\f$P_{i,j}(t)\f$).
 * Typically, this is computed using the formula
 * \f[
 * P(t) = e^{t \times Q},
 * \f]
 * where \f$P(t)\f$ is the matrix with all probabilities \f$P_{i,j}(t)\f$.
 * For some models, the \f$P_{i,j}(t)\f$'s can be computed analytically.
 * For more complexe models, we need to use a eigen-decomposition of \f$Q\f$:
 * \f[ Q = U^{-1} . D . U, \f]
 * where \f$D = diag(\lambda_i)\f$ is a diagonal matrix.
 * Hence
 * \f[
 * P(t) = e^{t \times Q} = U^{-1} . e^{D \times t} . U,
 * \f]
 * where \f$e^{D \times t} = diag\left(e^{\lambda_i \times t}\right)\f$ is a
 * diagonal matrix with all terms equal to exp the terms in \f$D\f$.
 * \f$U\f$ is the matrix of left eigen vectors (by row), and \f$U^{-1}\f$ is the matrix
 * of right eigen vectors (by column).
 * The values in \f$D\f$ are the eigen values of \f$Q\f$.
 * All \f$Q,U,U^{-1}\f$ and \f$D\f$ (its diagonal) may be retrieved from the
 * class (getEigenValues(), getRowRightEigenVectors() and getColumnLeftEigenVectors()
 * functions).
 *
 * First and second order derivatives of \f$P(t)\f$ with respect to \f$t\f$
 * can also be retrieved.
 * These methods may be useful for optimization processes.
 * Derivatives may be computed analytically, or using the general formulas:
 * \f[
 * \frac{\partial P(t)}{\partial t} = 
 * U^{-1} . diag\left(\lambda_i \times e^{\lambda_i \times t}\right) . U
 * \f]
 * and
 * \f[
 * \frac{\partial^2 P(t)}{\partial t^2} = 
 * U^{-1} . diag\left(\lambda_i^2 \times e^{\lambda_i \times t}\right) . U
 * \f]
 * 
 */

class SubstitutionModel:
  public virtual ParameterAliasable
{
	public:
		SubstitutionModel() {}
		virtual ~SubstitutionModel() {};

#ifndef NO_VIRTUAL_COV
    SubstitutionModel * clone() const = 0;
#endif

	public:
		
		/**
		 * @brief Get the name of the model.
		 *
		 * @return The name of this model.
		 */
		virtual string getName() const = 0;
	
		/**
		 * @return Equilibrium frequency associated to character i.
		 * @see getFrequencies()
		 */
		virtual double freq(int i) const = 0;

		/**
		 * @return The rate of change from state i to state j.
		 */
		virtual double Qij(int i, int j) const = 0;

		/**
		 * @return The probability of change from state i to state j during time t.
		 * @see getPij_t()
		 */
		virtual double Pij_t(int i, int j, double t) const = 0;

		/**
		 * @return The first order derivative of the probability of change from state
		 * i to state j with respect to time t, at time t.
		 * @see getdPij_dt()
		 */
		virtual double dPij_dt(int i, int j, double t) const = 0;
		
		/**
		 * @return The second order derivative of the probability of change from state
		 * i to state j with respect to time t, at time t.
		 * @see getd2Pij_dt2()
		 */
		virtual double d2Pij_dt2(int i, int j, double t) const = 0;
	
		/**
		 * @return A vector of all equilibrium frequencies.
		 * @see freq()
		 */
		virtual const Vdouble & getFrequencies() const = 0;
		
		/**
		 * @return The Markov generator matrix, i.e. all rates of changes from state i
		 * to state j. Usually, the generator is normalized so that
		 * (i) \f$ \forall j; \sum_i Q_{i,j} = 0 \f$, meaning that $\f$ \forall j; Q_{j,j} = -\sum_{i \neq j}Q_{i,j}\f$,
		 * and (ii) \f$ \sum_i Q_{i,i} \times \pi_i = -1\f$.
		 * This means that the mean rate of replacement at equilibrium is 1 and that time \f$t\f$ are measured
		 * in units of expected number of changes per site.
		 * 
		 * See Kosiol and Goldman (2005), Molecular Biology And Evolution 22(2) 193-9.
		 * @see Qij()
		 */ 
		virtual const Matrix<double> & getGenerator() const = 0;

		/**
		 * @return All probabilities of change from state i to state j during time t.
		 * @see Pij_t()
		 */
		virtual const Matrix<double> & getPij_t(double t) const = 0;
	
		/**
		 * @return Get all first order derivatives of the probability of change from state
		 * i to state j with respect to time t, at time t.
		 * @see dPij_dt()
		 */
		virtual const Matrix<double> & getdPij_dt(double t) const = 0;

		/**
		 * @return All second order derivatives of the probability of change from state
		 * i to state j with respect to time t, at time t.
		 * @see d2Pij_dt2()
		 */
		virtual const Matrix<double> & getd2Pij_dt2(double t) const = 0;

		/**
		 * @return A vector with all eigen values of the generator of this model;
		 */
		virtual const Vdouble & getEigenValues() const = 0;

		/**
		 * @return A matrix with left eigen vectors.
		 * Each row in the matrix stands for an eigen vector.
		 */
		virtual const Matrix<double> & getRowLeftEigenVectors() const = 0;

		/**
		 * @return A matrix with right eigen vectors.
		 * Each column in the matrix stands for an eigen vector.
		 */
		virtual const Matrix<double> & getColumnRightEigenVectors() const = 0;

		/**
		 * @return Get the alphabet associated to this model.
		 */
		virtual const Alphabet * getAlphabet() const = 0;

    /**
     * @brief Get the number of states.
     *
     * For most models, this equals the size of the alphabet.
     * 
     * @return The number of different states in the model.
     */
    virtual unsigned int getNumberOfStates() const = 0;

		/**
		 * This method is used to initialize likelihoods in reccursions.
		 * It typically sends 1 if i = state, 0 otherwise, where
		 * i is one of the possible states of the alphabet allowed in the model
		 * and state is the observed state in the considered sequence/site.
		 *
		 * @param i one of the possible states of the alphabet.
		 * @param state An observed state in the sequence/site.
		 * @return 1 or 0 depending if the two states are compatible.
		 * @throw BadIntException if states are not allowed in the associated alphabet.
		 */
		virtual double getInitValue(int i, int state) const throw (BadIntException) = 0;

		/**
		 * @brief Set equilibrium frequencies equal to the frequencies estimated
		 * from the data.
		 *
		 * @param data The sequences to use.
     * @param pseudoCount @f$\psi@f$ A quantity to add to adjust the observed values in order to prevent issues due to missing states on small data set.
     * The corrected frequencies shall be computed as
     * @f[
     * \pi_i = \frac{f_i+\psi}{4\psi + \sum_j f_j}
     * @f]
		 */
		virtual void setFreqFromData(const SequenceContainer & data, unsigned int pseudoCount = 0) = 0;

    /**
     * @brief Get the state in the alphabet corresponding to a given state in the model.
     *
     * In most cases, this method will return i.
     * @param i The state.
     * @return The corresponding state in the alphabet.
     * @see MarkovModulatedMarkovModels
     */
    virtual int getState(int i) const = 0;
		
};






/**
 * @brief Interface for reversible substitution models.
 * 
 * For reversible models,
 * \f[ Q = S . \pi, \f]
 * where \f$S\f$ is a symetric matrix called the exchangeability matrix, and
 * \f$\Pi\f$ the diagonal matrix with all equilibrium frequencies.
 * The frequencies may be retrieved as a vector by the getFrequencies() method
 * or individually by the freq() method.
 * The \f$S\f$ matrix may be obtained by the getExchangeabilityMatrix().
 */
class ReversibleSubstitutionModel:
  public virtual SubstitutionModel
{
  public:
		ReversibleSubstitutionModel() {}
		virtual ~ReversibleSubstitutionModel() {};

#ifndef NO_VIRTUAL_COV
    ReversibleSubstitutionModel * clone() const = 0;
#endif

	public:
	
    /**
		 * @return The matrix of exchangeability terms.
		 * It is recommended that exchangeability matrix be normalized so that the normalized 
		 * generator be obtained directly by the dot product \f$S . \pi\f$.
		 */
		virtual const Matrix<double> & getExchangeabilityMatrix() const = 0;

		/**
		 * @return The exchangeability between state i and state j.
     *
     * By definition Sij(i,j) = Sij(j,i).
		 */
    virtual double Sij(int i, int j) const = 0;

};

} //end of namespace bpp.

#endif	//_SUBSTITUTIONMODEL_H_

