//
// File: AbstractTreeParsimonyData.h
// Created by: Julien Dutheil
// Created on: Tue Jan 09 17:15 2007
// From file AbstractTreeParsimonyScore.h
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

#ifndef _ABSTRACTTREEPARSIMONYDATA_H_
#define _ABSTRACTTREEPARSIMONYDATA_H_

#include "TreeParsimonyData.h"

namespace bpp
{

/**
 * @brief Partial implementation of the TreeParsimonyData interface.
 *
 * This data structure provides a simple compression, by performing and storing computations
 * only one time per identical sites.
 *
 * The compression is achieved by the TreeParsimonyScore object.
 * The correspondance between sites in the dataset and the arrays in the structures is given
 * by the _rootPatternLinks array: the array indice for site @f$i@f$ if given by:
 * @code
 * _rootPatternLinks[i]
 * @endcode
 *
 * Finally, the _rootWeights array gives for each array position, the number of sites with this
 * pattern.
 * The global parsimony score is then given by the sum of all scores for each array position,
 * weighted by the corresponding number of sites.
 */
class AbstractTreeParsimonyData:
  public TreeParsimonyData
{
	protected:
		vector<unsigned int> _rootPatternLinks;
		vector<unsigned int> _rootWeights;
		TreeTemplate<Node> * _tree;
    
  public:
    AbstractTreeParsimonyData(): _rootPatternLinks(), _rootWeights(), _tree(NULL) {}
    virtual ~AbstractTreeParsimonyData() {}
    
	public:
		unsigned int getRootArrayPosition(const unsigned int site) const
		{
			return _rootPatternLinks[site];
		}
		unsigned int getWeight(unsigned int pos) const
		{ 
			return _rootWeights[pos];
		}
		const TreeTemplate<Node> * getTree() const { return _tree; }  
		TreeTemplate<Node> * getTree() { return _tree; }
};

} //end of namespace bpp.

#endif //_ABSTRACTTREEPARSIMONYDATA_H_

