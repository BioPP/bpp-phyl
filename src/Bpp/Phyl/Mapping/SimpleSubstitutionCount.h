//
// File: SimpleSubstitutionCount.h
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

#ifndef _SIMPLESUBSTITUTIONCOUNT_H_
#define _SIMPLESUBSTITUTIONCOUNT_H_

#include "SubstitutionCount.h"

#include <Bpp/Numeric/Matrix/Matrix.h>

namespace bpp
{

/**
 * @brief Naive substitution count.
 *
 * This substitution count is defined as follow:
 * - 0 if @f$i = j@f$,
 * - 1 if @f$i \neq j @f$.
 *
 * Reference (for instance):
 * Tufféry P, Darlu P.
 * Exploring a phylogenetic approach for the detection of correlated substitutions in proteins.
 * Mol Biol Evol. 2000 Nov;17(11):1753-9
 *
 * @author Julien Dutheil
 */
class SimpleSubstitutionCount:
  public AbstractSubstitutionCount
{
	public:
		SimpleSubstitutionCount(SubstitutionRegister* reg) :
      AbstractSubstitutionCount(reg) {}				
		
    SimpleSubstitutionCount(const SimpleSubstitutionCount& ssc) :
      AbstractSubstitutionCount(ssc) {}				
    
    SimpleSubstitutionCount& operator=(const SimpleSubstitutionCount& ssc)
    {
      AbstractSubstitutionCount::operator=(ssc);
      return *this;
    }				
		
    virtual ~SimpleSubstitutionCount() {}
			
	public:
		double getNumberOfSubstitutions(unsigned int initialState, unsigned int finalState, double length, unsigned int type = 1) const
    {
      return (register_->getType(initialState, finalState) == type ? 1. : 0.);
		}

    Matrix<double>* getAllNumbersOfSubstitutions(double length, unsigned int type = 1) const;
    
    std::vector<double> getNumberOfSubstitutionsForEachType(unsigned int initialState, unsigned int finalState, double length) const
    {
      std::vector<double> v(getNumberOfSubstitutionTypes());
      for (unsigned int t = 1; t <= getNumberOfSubstitutionTypes(); ++t) {
        v[t - 1] = (register_->getType(initialState, finalState) == t ? 1. : 0.);
      }
      return v;
    }
    
    void setSubstitutionModel(const SubstitutionModel* model) {}

};

/**
 * @brief Labelling substitution count.
 *
 * This substitution count return a distinct number for each possible mutation.
 * - 0 if @f$i = j@f$,
 * - @f$a(i,j)@f$ if @f$i \neq j @f$, where 'a' is an index giving a unique value for each combination of i and j.
 */
class LabelSubstitutionCount:
  public AbstractSubstitutionCount
{
  private:
    LinearMatrix<double> label_;
    
	public:
		LabelSubstitutionCount(const Alphabet* alphabet);

    virtual ~LabelSubstitutionCount() {}
			
	public:
		double getNumberOfSubstitutions(unsigned int initialState, unsigned int finalState, double length, unsigned int type = 1) const
    {
			return label_(initialState, finalState);
		}

    Matrix<double>* getAllNumbersOfSubstitutions(double length, unsigned int type = 1) const
    {
      return dynamic_cast<Matrix<double>*>(label_.clone());
    }
    
    std::vector<double> getNumberOfSubstitutionsForEachType(unsigned int initialState, unsigned int finalState, double length) const
    {
      std::vector<double> v(1);
      v[0] = label_(initialState, finalState);
      return v;
    }

    void setSubstitutionModel(const SubstitutionModel* model) {}

};

} //end of namespace bpp.

#endif // _SIMPLESUBSTITUTIONCOUNT_H_

