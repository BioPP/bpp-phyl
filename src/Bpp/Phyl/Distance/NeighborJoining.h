//
// File: NeighborJoining.h
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Thu jun 23 10:39 2005
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _NEIGHBORJOINING_H_
#define _NEIGHBORJOINING_H_

#include "AbstractAgglomerativeDistanceMethod.h"

namespace bpp
{

/**
 * @brief The neighbor joining distance method.
 *
 * Reference:
 * N Saitou and M Nei (1987), _Molecular Biology and Evolution_ 4(4) 406-25.
 */ 
class NeighborJoining :
  public AbstractAgglomerativeDistanceMethod
{
	protected:
    std::vector<double> sumDist_;
    bool positiveLengths_;
		
	public:
		NeighborJoining() :
      AbstractAgglomerativeDistanceMethod(),
      sumDist_(),
      positiveLengths_(false)
    {}

		NeighborJoining(const DistanceMatrix& matrix, bool rooted=false, bool positiveLengths=false) throw (Exception) :
      AbstractAgglomerativeDistanceMethod(matrix),
      sumDist_(),
      positiveLengths_(positiveLengths) 
		{
			sumDist_.resize(matrix.size());
			computeTree(rooted);
		}
   
		virtual ~NeighborJoining() {}

    NeighborJoining* clone() const { return new NeighborJoining(*this); }

	public:
    std::string getName() const { return "NJ"; }

		virtual void setDistanceMatrix(const DistanceMatrix& matrix)
		{ 
			AbstractAgglomerativeDistanceMethod::setDistanceMatrix(matrix);
			sumDist_.resize(matrix.size());
		}

    virtual void outputPositiveLengths(bool yn) { positiveLengths_ = yn; }
	
	protected:
		std::vector<unsigned int> getBestPair() throw (Exception);
		std::vector<double> computeBranchLengthsForPair(const std::vector<unsigned int>& pair);
		double computeDistancesFromPair(const std::vector<unsigned int>& pair, const std::vector<double>& branchLengths, unsigned int pos);
		void finalStep(int idRoot);	

};

} //end of namespace bpp.

#endif //_NEIGHBORJOINING_H_

