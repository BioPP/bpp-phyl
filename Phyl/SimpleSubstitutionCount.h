//
// File: SimpleSubstitutionCount.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 11:08 2006
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004, 2005, 2006)

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

#ifndef _SIMPLESUBSTITUTIONCOUNT_H_
#define _SIMPLESUBSTITUTIONCOUNT_H_

#include "SubstitutionCount.h"

// From NumCalc:
#include <NumCalc/Matrix.h>

namespace bpp
{

/**
 * @brief Naive substitution count.
 *
 * This subsitution count is defined as follow:
 * - 0 if @f$i = j@f$,
 * - 1 if @f$i \neq j @f$.
 *
 * Reference (for instance):
 * Tufféry P, Darlu P.
 * Exploring a phylogenetic approach for the detection of correlated substitutions in proteins.
 * Mol Biol Evol. 2000 Nov;17(11):1753-9
 */
class SimpleSubstitutionCount:
  public SubstitutionCount
{
  protected:
    const Alphabet * _alphabet;
    
	public:
		SimpleSubstitutionCount(const Alphabet * alphabet): _alphabet(alphabet) {}				
		virtual ~SimpleSubstitutionCount() {}
			
	public:
		double getNumberOfSubstitutions(int initialState, int finalState, double length) const
    {
			return initialState == finalState ? 0. : 1.;
		}

    virtual Matrix<double> * getAllNumbersOfSubstitutions(double length) const
    { 
      unsigned int n = _alphabet->getSize();
      RowMatrix<double> * mat = new RowMatrix<double>(n, n);
      for(unsigned int i = 0; i < n; i++) {
        (*mat)(i,i) = 0.;
        for(unsigned int j = 0; j < i; j++) {
          (*mat)(i,j) = (*mat)(j,i) = 1.;
        }
      }
      return mat;
    }

    void setSubstitutionModel(const SubstitutionModel* model) {}

};

} //end of namespace bpp.

#endif // _SIMPLESUBSTITUTIONCOUNT_H_

