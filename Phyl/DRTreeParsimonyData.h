//
// File: DRTreeParsimonyData.h
// Created by: Julien Dutheil
// Created on: Tue Jan 09 17:15 2007
// From file DRTreeParsimonyScore.h
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

#ifndef _DRTREEPARSIMONYDATA_H_
#define _DRTREEPARSIMONYDATA_H_

#include "AbstractTreeParsimonyData.h"

// From SeqLib
#include <Seq/SiteContainer.h>

// From the STL:
#include <bitset>

using namespace std;

namespace bpp
{

typedef bitset<20> Bitset;

/**
 * @brief Parsimony data structure for a node.
 * 
 * This class is for use with the DRTreeParsimonyData class.
 * 
 * Store for each neighbor node
 * - a vector of bitsets,
 * - a vector of score for the corresponding subtree.
 *
 * @see DRTreeParsimonyData
 */
class DRTreeParsimonyNodeData :
	public TreeParsimonyNodeData
{
	protected:
		mutable map<int, vector<Bitset> > _nodeBitsets;
		mutable map<int, vector<unsigned int> > _nodeScores;
		const Node * _node;

  public:
#ifndef NO_VIRTUAL_COV
    DRTreeParsimonyNodeData*
#else
    Clonable*
#endif
    clone() const { return new DRTreeParsimonyNodeData(*this); }

  public:
		const Node * getNode() const { return _node; }
		
    void setNode(const Node & node) { _node = &node; }

		vector<Bitset> & getBitsetsArrayForNeighbor(int neighborId)
		{
			return _nodeBitsets[neighborId];
		}
		const vector<Bitset> & getBitsetsArrayForNeighbor(int neighborId) const
		{
			return _nodeBitsets[neighborId];
		}
		vector<unsigned int> & getScoresArrayForNeighbor(int neighborId)
		{
			return _nodeScores[neighborId];
		}
		const vector<unsigned int> & getScoresArrayForNeighbor(int neighborId) const
		{
			return _nodeScores[neighborId];
		}

		bool isNeighbor(int neighborId) const
		{
			return _nodeBitsets.find(neighborId) != _nodeBitsets.end();
		}

		void eraseNeighborArrays()
		{
			_nodeBitsets.erase(_nodeBitsets.begin(), _nodeBitsets.end());
			_nodeScores.erase(_nodeScores.begin(), _nodeScores.end());
		}
};

/**
 * @brief Parsimony data structure for a leaf.
 * 
 * This class is for use with the DRTreeParsimonyData class.
 * 
 * Store the vector of bitsets associated to a leaf.
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
#ifndef NO_VIRTUAL_COV
    DRTreeParsimonyLeafData*
#else
    Clonable*
#endif
    clone() const { return new DRTreeParsimonyLeafData(*this); }

	public:
		const Node * getNode() const { return _leaf; }
		void setNode(const Node & node) { _leaf = &node; }

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
	public AbstractTreeParsimonyData
{
	protected:
		mutable map<int, DRTreeParsimonyNodeData> _nodeData;
		mutable map<int, DRTreeParsimonyLeafData> _leafData;
		mutable vector<Bitset> _rootBitsets;
		mutable vector<unsigned int> _rootScores;
		SiteContainer * _shrunkData;
		unsigned int _nbSites; 
		unsigned int _nbStates;
		unsigned int _nbDistinctSites; 

	public:
		DRTreeParsimonyData(TreeTemplate<Node> & tree):
      _nodeData(), _leafData(), _rootBitsets(), _rootScores(), _shrunkData(NULL),
      _nbSites(0), _nbStates(0), _nbDistinctSites(0)
    {
      _tree = &tree;
    }

    DRTreeParsimonyData(const DRTreeParsimonyData & data);
    
    DRTreeParsimonyData & operator=(const DRTreeParsimonyData & data);
		
    virtual ~DRTreeParsimonyData() { delete _shrunkData; }

#ifndef NO_VIRTUAL_COV
    DRTreeParsimonyData*
#else
    Clonable*
#endif
    clone() const { return new DRTreeParsimonyData(*this); }

	public:
    void setTree(TreeTemplate<Node> & tree)
    { 
      _tree = &tree;
      for(map<int, DRTreeParsimonyNodeData>::iterator it = _nodeData.begin(); it != _nodeData.end(); it++)
      {
        int id = it->second.getNode()->getId();
        it->second.setNode(*_tree->getNode(id));
      }
      for(map<int, DRTreeParsimonyLeafData>::iterator it = _leafData.begin(); it != _leafData.end(); it++)
      {
        int id = it->second.getNode()->getId();
        it->second.setNode(*_tree->getNode(id));
      }
    }

		DRTreeParsimonyNodeData & getNodeData(int nodeId)
		{ 
			return _nodeData[nodeId];
		}
		const DRTreeParsimonyNodeData & getNodeData(int nodeId) const
		{ 
			return _nodeData[nodeId];
		}
		
		DRTreeParsimonyLeafData & getLeafData(int nodeId)
		{ 
			return _leafData[nodeId];
		}
		const DRTreeParsimonyLeafData & getLeafData(int nodeId) const
		{ 
			return _leafData[nodeId];
		}
		
		vector<Bitset> & getBitsetsArray(int nodeId, int neighborId)
		{ 
			return _nodeData[nodeId].getBitsetsArrayForNeighbor(neighborId);
		}
		const vector<Bitset> & getBitsetsArray(int nodeId, int neighborId) const
		{ 
			return _nodeData[nodeId].getBitsetsArrayForNeighbor(neighborId);
		}
		
		vector<unsigned int> & getScoresArray(int nodeId, int neighborId)
		{ 
			return _nodeData[nodeId].getScoresArrayForNeighbor(neighborId);
		}
		const vector<unsigned int> & getScoresArray(int nodeId, int neighborId) const 
		{
			return _nodeData[nodeId].getScoresArrayForNeighbor(neighborId);
		}

		unsigned int getArrayPosition(int parentId, int sonId, unsigned int currentPosition) const
		{ 
			return currentPosition;
		}

    vector<Bitset> & getRootBitsets() { return _rootBitsets; }
    const vector<Bitset> & getRootBitsets() const { return _rootBitsets; }
    const Bitset & getRootBitset(unsigned int i) const { return _rootBitsets[i]; }

		vector<unsigned int> & getRootScores() { return _rootScores; }
		const vector<unsigned int> & getRootScores() const { return _rootScores; }
		unsigned int getRootScore(unsigned int i) const { return _rootScores[i]; }

		unsigned int getNumberOfDistinctSites() const { return _nbDistinctSites; }
		unsigned int getNumberOfSites() const { return _nbSites; }
		unsigned int getNumberOfStates() const { return _nbStates; }
		
		void init(const SiteContainer & sites) throw (Exception);
		void reInit() throw (Exception);

	protected:
		void init(const Node * node, const SiteContainer & sites) throw (Exception);
		void reInit(const Node * node) throw (Exception);
		
};

} //end of namespace bpp.

#endif //_DRTREEPARSIMONYDATA_H_

