//
// File: DRTreeParsimonyScore.h
// Created by: Julien Dutheil
// Created on: Thu Jul 28 18:31 2005
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

#ifndef _DRTREEPARSIMONYSCORE_H_
#define _DRTREEPARSIMONYSCORE_H_

#include "AbstractTreeParsimonyScore.h"

// From the STL:
#include <bitset>
using namespace std;

typedef bitset<20> Bitset;

/**
 * @brief Parsimony data structure for a node.
 * 
 * This class is for use with the DRTreeParsimonyData class.
 * 
 * We store for each neighbor node
 * - a vector of bitsets,
 * - a vector of score for the corresponding subtree.
 *
 * @see DRTreeParsimonyData
 */
class DRTreeParsimonyNodeData :
	public TreeParsimonyNodeData
{
	protected:
		mutable map<const Node *, vector<Bitset> > _nodeBitsets;
		mutable map<const Node *, vector<unsigned int> > _nodeScores;
		const Node * _node;

	public:
		const Node * getNode() const { return _node; }
		void setNode(const Node * node) { _node = node; }

		vector<Bitset> & getBitsetsArrayForNeighbor(const Node * neighbor)
		{
			return _nodeBitsets[neighbor];
		}
		const vector<Bitset> & getBitsetsArrayForNeighbor(const Node * neighbor) const
		{
			return _nodeBitsets[neighbor];
		}
		vector<unsigned int> & getScoresArrayForNeighbor(const Node * neighbor)
		{
			return _nodeScores[neighbor];
		}
		const vector<unsigned int> & getScoresArrayForNeighbor(const Node * neighbor) const
		{
			return _nodeScores[neighbor];
		}
};

/**
 * @brief Parsimony data structure for a leaf.
 * 
 * This class is for use with the DRTreeParsimonyData class.
 * 
 * We store the vector of bitsets associated to this leaf.
 *
 * @see DRTreeParsimonyData
 */
class DRTreeParsimonyLeafData :
	public TreeParsimonyNodeData
{
	protected:
		mutable vector<Bitset> _leafBitsets;
		const Node * _leaf;

	public:
		const Node * getNode() const { return _leaf; }
		void setNode(const Node * node) { _leaf = node; }

		vector<Bitset> & getBitsetsArray()
		{
			return _leafBitsets;
		}
		const vector<Bitset> & getBitsetsArray() const
		{
			return _leafBitsets;
		}
};

/**
 * @brief Parsimony data structure for double-recursive (DR) algorithm.
 *
 * States are coded using bitsets for faster computing (@see AbstractTreeParsimonyData).
 * For each inner node in the tree, we store a DRTreeParsimonyNodeData object in _nodeData.
 * For each leaf node in the tree, we store a DRTreeParsimonyLeafData object in _leafData.
 * 
 * The dataset is first compressed, removing all identical sites.
 * The resulting dataset is stored in _shrunkData.
 * The corresponding positions are stored in _rootPatternLinks, inherited from AbstractTreeParsimonyData.
 */
class DRTreeParsimonyData :
	public virtual AbstractTreeParsimonyData
{
	protected:
		TreeTemplate<Node> * _tree;
		mutable map<const Node *, DRTreeParsimonyNodeData> _nodeData;
		mutable map<const Node *, DRTreeParsimonyLeafData> _leafData;
		mutable vector<Bitset> _rootBitsets;
		mutable vector<unsigned int> _rootScores;
		SiteContainer * _shrunkData;
		unsigned int _nbSites; 
		unsigned int _nbStates;
		unsigned int _nbDistinctSites; 

	public:
		DRTreeParsimonyData(TreeTemplate<Node> & tree) : _tree(& tree) {}

	public:
		DRTreeParsimonyNodeData & getNodeData(const Node * node)
		{ 
			return _nodeData[node];
		}
		const DRTreeParsimonyNodeData & getNodeData(const Node * node) const
		{ 
			return _nodeData[node];
		}
		
		DRTreeParsimonyLeafData & getLeafData(const Node * node)
		{ 
			return _leafData[node];
		}
		const DRTreeParsimonyLeafData & getLeafData(const Node * node) const
		{ 
			return _leafData[node];
		}
		
		vector<Bitset> & getBitsetsArray(const Node * node, const Node * neighbor)
		{ 
			return _nodeData[node].getBitsetsArrayForNeighbor(neighbor);
		}
		const vector<Bitset> & getBitsetsArray(const Node * node, const Node * neighbor) const
		{ 
			return _nodeData[node].getBitsetsArrayForNeighbor(neighbor);
		}
		
		vector<unsigned int> & getScoresArray(const Node * node, const Node * neighbor)
		{ 
			return _nodeData[node].getScoresArrayForNeighbor(neighbor);
		}
		const vector<unsigned int> & getScoresArray(const Node * node, const Node * neighbor) const 
		{
			return _nodeData[node].getScoresArrayForNeighbor(neighbor);
		}

		unsigned int getArrayPosition(const Node* parent, const Node* son, unsigned int currentPosition) const
		{ 
			return currentPosition;
		}

		const TreeTemplate<Node> * getTree() const { return _tree; } 

		unsigned int getNumberOfDistinctSites() const { return _nbDistinctSites; }
		unsigned int getNumberOfSites() const { return _nbSites; }
		unsigned int getNumberOfStatees() const { return _nbStates; }
		
		void init(const SiteContainer & sites) throw (Exception);
		void init(const Node * node, const SiteContainer & sites) throw (Exception);
};

/**
 * @brief Double recursive implementation of interface TreeParsimonyScore.
 *
 * Uses a DRTreeParsimonyData object for data storage.
 */
class DRTreeParsimonyScore :
	public virtual AbstractTreeParsimonyScore
{
	protected:
		DRTreeParsimonyData *_parsimonyData;
		vector<unsigned int> _rootScores;
		vector<Bitset>       _rootBitsets;
		unsigned int         _nbDistinctSites;
			
	public:
		DRTreeParsimonyScore(
			TreeTemplate<Node> & tree,
			const SiteContainer & data,
			bool verbose = true)
			throw (Exception);
				
		virtual ~DRTreeParsimonyScore();

	public:

	protected:
		/**
		 * @brief Compute all scores.
		 *
		 * Call the computeScoresPreorder and computeScoresPostorder methods, and then initialize _rootBitsets and _rootScores.
		 */
		virtual void computeScores();
		/**
		 * @brief Compute scores (preorder algorithm).
		 */
		virtual void computeScoresPreorder(const Node *);
		/**
		 * @brief Compute scores (postorder algorithm).
		 */
		virtual void computeScoresPostorder(const Node *);

	public:

		unsigned int getScore() const;
		unsigned int getScoreForSite(unsigned int site) const;
		
		/**
		 * @brief Compute bitsets and scores for each site for a node, in postorder.
		 *
		 * @param pData    The node data to use.
		 * @param rBitsets The bitset array where to store the resulting bitsets.
		 * @param rScores  The score array where to write the resulting scores.
		 */
		static void computeScoresPostorderForNode(const DRTreeParsimonyNodeData & pData, vector<Bitset> & rBitsets, vector<unsigned int> & rScores);
		
		/**
		 * @brief Compute bitsets and scores for each site for a node, in preorder.
		 *
		 * @param pData    The node data to use.
		 * @param source   The node where we are coming from.
		 * @param rBitsets The bitset array where to store the resulting bitsets.
		 * @param rScores  The score array where to write the resulting scores.
		 */
		static void computeScoresPreorderForNode(const DRTreeParsimonyNodeData & pData, const Node * source, vector<Bitset> & rBitsets, vector<unsigned int> & rScores);

		/**
		 * @brief Compute bitsets and scores for each site for a node, in all directions.
		 *
		 * @param pData    The node data to use.
		 * @param rBitsets The bitset array where to store the resulting bitsets.
		 * @param rScores  The score array where to write the resulting scores.
		 */
		static void computeScoresForNode(const DRTreeParsimonyNodeData & pData, vector<Bitset> & rBitsets, vector<unsigned int> & rScores); 

};

#endif //_ABSTRACTTREEPARSIMONYSCORE_H_

