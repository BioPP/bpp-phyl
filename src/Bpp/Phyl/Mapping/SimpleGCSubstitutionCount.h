//
// File: SimpleGCSubstitutionCount.h
// Created by: Julien Dutheil
// Created on: Thu Dec 02 16:23 2010
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

#ifndef _SIMPLEGCSUBSTITUTIONCOUNT_H_
#define _SIMPLEGCSUBSTITUTIONCOUNT_H_

#include "SubstitutionCount.h"

#include <Bpp/Numeric/Matrix/Matrix.h>

namespace bpp
{

/**
 * @brief Simple GC substitution count, ignoring multiple substitutions.
 *
 * This substitution count is defined as follow:
 * - 1 if i is A or T and j is G or C and type is 0,
 * - 1 if i is G or C and j is A or T and type is 1,
 * - 0 in all other cases.
 *
 * @see SimpleSubstitutionCount
 *
 * @author Julien Dutheil
 */
class SimpleGCSubstitutionCount:
  public virtual SubstitutionCount
{
  private:
    const NucleicAlphabet* alphabet_;
    
	public:
		SimpleGCSubstitutionCount(const NucleicAlphabet* alphabet) : alphabet_(alphabet) {}				
		
    SimpleGCSubstitutionCount(const SimpleGCSubstitutionCount& ssc) : alphabet_(ssc.alphabet_) {}				
    
    SimpleGCSubstitutionCount& operator=(const SimpleGCSubstitutionCount& ssc)
    {
      alphabet_ = ssc.alphabet_;
      return *this;
    }				
		
    virtual ~SimpleGCSubstitutionCount() {}
			
	public:
		double getNumberOfSubstitutions(unsigned int initialState, unsigned int finalState, double length, unsigned int type = 0) const
    {
			return initialState == finalState ? 0. : 1.;
		}

    Matrix<double>* getAllNumbersOfSubstitutions(double length, unsigned int type = 0) const;
    
    std::vector<double> getNumberOfSubstitutionsForEachType(unsigned int initialState, unsigned int finalState, double length) const
    {
      std::vector<double> v(1);
      v[0] = (initialState == finalState ? 0. : 1.);
      return v;
    }
    
    unsigned int getSubstitutionType(unsigned int initialState, unsigned int finalState) const throw (Exception) {
      if (initialState == finalState)
        throw Exception("SimpleGCSubstitutionCount::getSubstitutionType. Not a substitution!");
      if ((initialState == 1 || initialState == 2) && (finalState == 0 || finalState == 3))
        return 0;
      if ((initialState == 0 || initialState == 3) && (finalState == 1 || finalState == 2))
        return 1;
      throw Exception("SimpleGCSubstitutionCount::getSubstitutionType. Not a supported substitution! (either G<->C or A<->T, which are ignored.");
    }
    unsigned int getNumberOfSubstitutionTypes() const { return 2; }

    void setSubstitutionModel(const SubstitutionModel* model) {}

};

} //end of namespace bpp.

#endif // _SIMPLEGCSUBSTITUTIONCOUNT_H_

