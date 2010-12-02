//
// File: SubstitutionCount.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 11:08 2006
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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

#ifndef _SUBSTITUTIONCOUNT_H_
#define _SUBSTITUTIONCOUNT_H_

#include "../Model/SubstitutionModel.h"

#include <Bpp/Numeric/Matrix/Matrix.h>

//From the STL:
#include <vector>

namespace bpp
{

/**
 * @brief The SubstitutionsCount interface.
 *
 * Provide a method to compute the @f$n_{x,y}(t)@f$ function,
 * namely the number of substitutions on a branch of length @f$t@f$, with initial state @f$x@f$ and final state @f$y@f$.
 * 
 * The new implementation offers to perform several counts simultaneously, distinguishing between different types of substitutions.
 * The types are registered as unsigned int, so that all counts can be retrieved in one go as a vector, the type serving as an indice.
 * The type indice for a given substitution can be obtained using the getSubstitutionType method.
 *
 * @author Julien Dutheil
 *
 * See:
 * Dutheil J, Pupko T, Jean-Marie A, Galtier N.
 * A model-based approach for detecting coevolving positions in a molecule.
 * Mol Biol Evol. 2005 Sep;22(9):1919-28. Epub 2005 Jun 8.
 */
class SubstitutionCount
{
	public:
		SubstitutionCount() {}
		virtual ~SubstitutionCount() {}
	
	public:
		/**
		 * @brief Get the number of susbstitutions on a branch, given the initial and final states, and the branch length.
		 *
		 * @param initialState The initial state.
		 * @param finalState   The final state.
		 * @param length       The length of the branch.
     * @param type         The type of susbstitution to count.
		 * @return The number of substitutions on a branch of specified length and
		 * according to initial and final states.
		 */
		virtual double getNumberOfSubstitutions(unsigned int initialState, unsigned int finalState, double length, unsigned int type) const = 0;
		
		/**
		 * @brief Get the numbers of susbstitutions on a branch, for each initial and final states, and given the branch length.
		 *
		 * @param length       The length of the branch.
     * @param type         The type of susbstitution to count.
		 * @return A matrix with all numbers of substitutions for each initial and final states.
		 */
    virtual Matrix<double>* getAllNumbersOfSubstitutions(double length, unsigned int type) const = 0;

    /**
		 * @brief Get the numbers of susbstitutions on a branch for all types, for an initial and final states, given the branch length.
		 *
		 * @param initialState The initial state.
		 * @param finalState   The final state.
		 * @param length       The length of the branch.
		 * @return A matrix with all numbers of substitutions for each initial and final states.
		 */
    virtual std::vector<double> getNumberOfSubstitutionsForEachType(unsigned int initialState, unsigned int finalState, double length) const = 0;

    /**
     * @brief Get the substitution type according to initial and final states indices.
     *
		 * @param initialState The initial state.
		 * @param finalState   The final state.
     */
    virtual unsigned int getSubstitutionType(unsigned int initialState, unsigned int finalState) const = 0;

    /**
     * @return The number of substitution types supported by this class.
     */
    virtual unsigned int getNumberOfSubstitutionTypes() const = 0;

    /**
     * @brief Set the substitution model associated with this count, if relevent.
     *
     * @param model The substitution model to use with this count.
     */
    virtual void setSubstitutionModel(const SubstitutionModel* model) = 0;
};

} //end of namespace bpp.

#endif //_SUBSTITUTIONCOUNT_H_

