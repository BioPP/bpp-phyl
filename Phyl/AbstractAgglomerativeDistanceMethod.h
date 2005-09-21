//
// File: AbstractAgglomerativeDistanceMethod.h
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Wed jun 22 10:00 2005
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

#include "AgglomerativeDistanceMethod.h"
#include "DistanceMatrix.h"
#include "Node.h"
#include "TreeTemplate.h"

// From the STL:
#include <map>
using namespace std;

class AbstractAgglomerativeDistanceMethod : public virtual AgglomerativeDistanceMethod
{
	protected:
		DistanceMatrix _matrix;
		Tree * _tree;
		/**
		 * @brief key: indice of the node in the distance matrix.
		 *        value: number of species represented by this node.
		 */
	
		map<unsigned int, Node *> _currentNodes;
	
	public:
		AbstractAgglomerativeDistanceMethod(): _matrix(0), _tree(NULL) {}
		AbstractAgglomerativeDistanceMethod(const DistanceMatrix & matrix): _matrix(matrix) {}
		virtual ~AbstractAgglomerativeDistanceMethod();

	public:
		void setDistanceMatrix(const DistanceMatrix & matrix);

#if defined(VIRTUAL_COV)
		TreeTemplate<Node> * 
#else
		Tree *
#endif
		getTree() const;
		
		void computeTree(bool rooted);
	

	protected:
		virtual vector<unsigned int> getBestPair() = 0;
		virtual vector<double> computeBranchLengthsForPair(const vector<unsigned int> & pair) = 0;
		virtual double computeDistancesFromPair(const vector<unsigned int> & pair, const vector<double> & branchLengths, unsigned int pos) = 0;
		virtual void finalStep(int idRoot) = 0;
		virtual Node * getLeafNode(int id, const string & name);
		virtual Node * getParentNode(int id, Node * son1, Node * son2);
		
};

