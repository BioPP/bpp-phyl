//
// File: AnalyticalSubstitutionCount.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 11:21 2006
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

#ifndef _ANALYTICALSUBSTITUTIONCOUNT_H_
#define _ANALYTICALSUBSTITUTIONCOUNT_H_

#include "SubstitutionCount.h"
#include "../Model/SubstitutionModel.h"

namespace bpp
{

/**
 * @brief Analytical estimate of the substitution count.
 *
 * This method uses Laplace transforms, as described in 
 * Dutheil J, Pupko T, Jean-Marie A, Galtier N.
 * A model-based approach for detecting coevolving positions in a molecule.
 * Mol Biol Evol. 2005 Sep;22(9):1919-28.
 *
 * @author Julien Dutheil
 */
class AnalyticalSubstitutionCount:
  public AbstractSubstitutionCount
{
	private:
		const SubstitutionModel* model_;
		int cutOff_;
		mutable double currentLength_;
		mutable RowMatrix<double> m_;
	
	public:
		AnalyticalSubstitutionCount(const SubstitutionModel* model, int cutOff) :
      AbstractSubstitutionCount(new TotalSubstitutionRegister(model->getAlphabet())),
      model_        (model),
      cutOff_       (cutOff),
      currentLength_(-1),
      m_            (model->getNumberOfStates(), model->getNumberOfStates())
    {}
	
    AnalyticalSubstitutionCount(const AnalyticalSubstitutionCount& asc) :
      AbstractSubstitutionCount(asc),
      model_        (asc.model_),
      cutOff_       (asc.cutOff_),
      currentLength_(asc.currentLength_),
      m_            (asc.m_)
    {}
				
	  AnalyticalSubstitutionCount& operator=(const AnalyticalSubstitutionCount& asc)
    {
      AbstractSubstitutionCount::operator=(asc);
      model_         = asc.model_;
      cutOff_        = asc.cutOff_;
      currentLength_ = asc.currentLength_;
      m_             = asc.m_;
      return *this;
    }
				
		virtual ~AnalyticalSubstitutionCount() {}
			
	public:
		double getNumberOfSubstitutions(unsigned int initialState, unsigned int finalState, double length, unsigned int type = 1) const;
    Matrix<double>* getAllNumbersOfSubstitutions(double length, unsigned int type = 1) const;
    std::vector<double> getNumberOfSubstitutionsForEachType(unsigned int initialState, unsigned int finalState, double length) const
    {
      std::vector<double> v(0);
      v[0] = getNumberOfSubstitutions(initialState, finalState, length, 0);
      return v;
    }

    void setSubstitutionModel(const SubstitutionModel* model);

  protected:
    void computeCounts(double length) const;
};

} //end of namespace bpp.

#endif //_ANALYTICALSUBSTITUTIONCOUNT_H_

