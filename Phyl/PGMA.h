//
// File: PGMA.h
// Created by: Julien Dutheil
// Created on: Mon jul 11 11:41 2005
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

#ifndef _PGMA_H_
#define _PGMA_H_

#include "AbstractAgglomerativeDistanceMethod.h"
#include "Tree.h"
#include "TreeTemplate.h"

namespace bpp
{

/**
 * @brief Inner data structure for WPGMA and UPGMA distance methods.
 */
struct PGMAInfos
{
	unsigned int numberOfLeaves;
	double time;
};

/**
 * @brief Compute WPGMA and UPGMA trees from a distance matrix.
 * 
 * WPGMA = Weighted pair group method using arithmetic averaging,
 * is equivalent to the average linkegage hierarchical clustering method.
 * The distance between two taxa is the average distance between all individuals in each taxa.
 * The unweighted version (named UPGMA), uses a weighted average, with the number of individuals in a group as a weight.
 */
class PGMA:
  public AbstractAgglomerativeDistanceMethod
{
	protected:
		bool _weighted;
		
	public:
		PGMA(bool weighted = true): _weighted(weighted) {}

		PGMA(const DistanceMatrix & matrix, bool weighted = true, bool verbose = true) throw (Exception):
      AbstractAgglomerativeDistanceMethod(matrix, verbose), _weighted(weighted)
		{
			computeTree(true);
		}
		virtual ~PGMA() {}

	public:
		void setDistanceMatrix(const DistanceMatrix & matrix)
		{ 
			AbstractAgglomerativeDistanceMethod::setDistanceMatrix(matrix);
		}

		TreeTemplate<Node> * getTree() const;
		
		void setWeighted(bool weighted) { _weighted = weighted; }
		bool isWeighted() const { return _weighted; }
	
	protected:
		vector<unsigned int> getBestPair() throw (Exception);
		vector<double> computeBranchLengthsForPair(const vector<unsigned int> & pair);
		double computeDistancesFromPair(const vector<unsigned int> & pair, const vector<double> & branchLengths, unsigned int pos);
		void finalStep(int idRoot);	
		virtual Node * getLeafNode(int id, const string & name);
		virtual Node * getParentNode(int id, Node * son1, Node * son2);

};

} //end of namespace bpp.

#endif // _PGMA_H_

