//
// File: AbstractTreeParsimonyScore.h
// Created by: Julien Dutheil
// Created on: Thu Jul 28 17:25 2005
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

#ifndef _ABSTRACTTREEPARSIMONYSCORE_H_
#define _ABSTRACTTREEPARSIMONYSCORE_H_

#include "TreeParsimonyScore.h"
#include "Node.h"

// From SeqLib:
#include <Seq/SiteContainer.h>

/**
 * @brief TreeParsimonyScore partial data structure.
 *
 * Stores inner computation for a given node.
 *
 * @see TreeParsimonyData
 */
class TreeParsimonyNodeData
{
	public:
		TreeParsimonyNodeData() {}
		~TreeParsimonyNodeData() {}

	public:
		/**
		 * @brief Get the node associated to this data structure.
		 *
		 * @return The node associated to this structure.
		 */
		virtual const Node * getNode() const = 0;
};

/**
 * @brief TreeParsimonyScore data structure.
 *
 * Stores all the inner computations:
 * - subtree scores and ancestral states for each node,
 * - correspondance between sites in the dataset and array indices.
 *
 * @see TreeParsimonyNodeData
 */
class TreeParsimonyData
{
	public:
		virtual const TreeTemplate<Node> * getTree() const = 0;  
		virtual TreeTemplate<Node> * getTree() = 0;
		virtual unsigned int getArrayPosition(const Node * parent, const Node * son, unsigned int currentPosition) const = 0;
		virtual unsigned int getRootArrayPosition(const unsigned int site) const = 0;
		virtual TreeParsimonyNodeData & getNodeData(const Node * node) = 0;
		virtual const TreeParsimonyNodeData & getNodeData(const Node * node) const = 0;
};

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
class AbstractTreeParsimonyData : public TreeParsimonyData
{
	protected:
		vector<unsigned int> _rootPatternLinks;
		vector<unsigned int> _rootWeights;
		TreeTemplate<Node> * _tree;

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

/**
 * @brief Low-level implementation of interface TreeParsimonyScore.
 */
class AbstractTreeParsimonyScore :
	public virtual TreeParsimonyScore
{
	protected:
		TreeTemplate<Node> * _tree;
		const SiteContainer * _data;
		const Alphabet * _alphabet;
		unsigned int _nbStates;
		
	public:
		AbstractTreeParsimonyScore(
			TreeTemplate<Node> & tree,
			const SiteContainer & data,
			bool verbose)
			throw (Exception);

		virtual ~AbstractTreeParsimonyScore() {}

	public:
		virtual Tree * getTree() { return _tree; }
		virtual const Tree * getTree() const { return _tree; }
		virtual vector<unsigned int> getScoreForEachSite() const;

};

#endif //_ABSTRACTTREEPARSIMONYSCORE_H_

