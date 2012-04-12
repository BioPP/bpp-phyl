//
// File: WeightedSubstitutionCount.h
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

#ifndef _WEIGHTEDSUBSTITUTIONCOUNT_H_
#define _WEIGHTEDSUBSTITUTIONCOUNT_H_

#include "SubstitutionCount.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/StateProperties/AlphabetIndex2.h>

namespace bpp
{

/**
 * @brief Weight substitution counts according to state properties.
 *
 * This class uses a ComprehensiveSubstitutionRegister to map all types of substitutions,
 * and multiplies each counts by a weight specified via an AlphabetIndex2 object.
 * Weighted counts are then summed according the specified substitution register.
 */
class WeightedSubstitutionCount:
  public AbstractSubstitutionCount
{
	private:
    SubstitutionCount* subCount_;
		const AlphabetIndex2<double>* dist_;
		bool ownDist_;
	
	public:
		WeightedSubstitutionCount(SubstitutionCount* subCount, SubstitutionRegister* reg, const AlphabetIndex2<double>* ai2, bool ownDistance) :
      AbstractSubstitutionCount(reg), 
      subCount_(subCount),
      dist_(ai2),
      ownDist_(ownDistance)
    {
      if (!subCount)
        throw NullPointerException("WeightedSubstitutionCount, constructor: input pointer should not be null.");
      subCount_->setSubstitutionRegister(new ComprehensiveSubstitutionRegister(subCount->getAlphabet(), false));
    }

    WeightedSubstitutionCount(const WeightedSubstitutionCount& index) :
      AbstractSubstitutionCount(index),
      subCount_(index.subCount_),
      dist_(index.dist_),
      ownDist_(index.ownDist_)
    {
      subCount_ = dynamic_cast<SubstitutionCount*>(index.subCount_->clone());
      if (ownDist_)
        dist_ = dynamic_cast<AlphabetIndex2<double>*>(index.dist_->clone());
    }

    WeightedSubstitutionCount& operator=(const WeightedSubstitutionCount& index)
    {
      AbstractSubstitutionCount::operator=(index);

      subCount_ = dynamic_cast<SubstitutionCount*>(index.subCount_->clone());
      
      ownDist_ = index.ownDist_;
      if (ownDist_) dist_ = dynamic_cast<AlphabetIndex2<double>*>(index.dist_->clone());
      else dist_ = index.dist_;
      
      return *this;
    }
		
		virtual ~WeightedSubstitutionCount()
    {
      delete subCount_;
			if (ownDist_)
        delete dist_;
		}
		
    WeightedSubstitutionCount* clone() const { return new WeightedSubstitutionCount(*this); }

	public:
		double getNumberOfSubstitutions(unsigned int initialState, unsigned int finalState, double length, unsigned int type = 1) const;

		Matrix<double>* getAllNumbersOfSubstitutions(double length, unsigned int type = 1) const;

    std::vector<double> getNumberOfSubstitutionsForEachType(unsigned int initialState, unsigned int finalState, double length) const;
    
    void setSubstitutionModel(const SubstitutionModel* model) { subCount_->setSubstitutionModel(model); }

		const SubstitutionRegister* getSubstitutionRegister() const { return subCount_->getSubstitutionRegister(); }
		
  public:
    const AlphabetIndex2<double>* getAlphabetIndex2() const { return dist_; }

  protected:
    void substitutionRegisterHasChanged() {}
};

} //end of namespace bpp.

#endif //_INDEXTOCOUNT_H_

