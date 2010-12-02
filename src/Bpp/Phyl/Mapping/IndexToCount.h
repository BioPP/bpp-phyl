//
// File: IndexToCount.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 15:12 2006
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

#ifndef _INDEXTOCOUNT_H_
#define _INDEXTOCOUNT_H_

#include "SubstitutionCount.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/StateProperties/AlphabetIndex2.h>

namespace bpp
{

/**
 * @brief Naive substitution count weighted by amino-acids properties.
 *
 * This class uses a AlphabetIndex2 object to weight substitutions.
 */
class IndexToCount:
  public SubstitutionCount
{
	private:
		const AlphabetIndex2<double> *dist_;
		bool ownDist_;
	
	public:
		IndexToCount(const AlphabetIndex2<double>* ai2, bool ownDistance) :
			dist_(ai2), ownDist_(ownDistance)
    {}

    IndexToCount(const IndexToCount& index) :
      dist_(index.dist_), ownDist_(index.ownDist_)
    {
      if (ownDist_) dist_ = dynamic_cast<AlphabetIndex2<double>*>(index.dist_->clone());
    }

    IndexToCount& operator=(const IndexToCount& index)
    {
      ownDist_ = index.ownDist_;
      if (ownDist_) dist_ = dynamic_cast<AlphabetIndex2<double>*>(index.dist_->clone());
      else dist_ = index.dist_;
      return *this;
    }
		
		virtual ~IndexToCount()
    {
			if (ownDist_) delete dist_;
		}
			
	public:

		double getNumberOfSubstitutions(unsigned int initialState, unsigned int finalState, double length, unsigned int type = 0) const
    {
			return dist_->getIndex(initialState, finalState);
		}

		Matrix<double>* getAllNumbersOfSubstitutions(double length, unsigned int type = 0) const
    {
      return dist_->getIndexMatrix();
    }

    std::vector<double> getNumberOfSubstitutionsForEachType(unsigned int initialState, unsigned int finalState, double length) const
    {
      std::vector<double> v(0);
      v[0] = getNumberOfSubstitutions(initialState, finalState, length, 0);
      return v;
    }
    
    unsigned int getSubstitutionType(unsigned int initialState, unsigned int finalState) const { return 0; }
    unsigned int getNumberOfSubstitutionTypes() const { return 1; }

    void setSubstitutionModel(const SubstitutionModel* model) {}

	public:
		const AlphabetIndex2<double>* getAlphabetIndex2() const { return dist_; }
};

} //end of namespace bpp.

#endif //_INDEXTOCOUNT_H_

