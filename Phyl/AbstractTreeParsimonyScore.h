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

class TreeParsimonyNodeData
{
	public:
		TreeParsimonyNodeData() {}
		~TreeParsimonyNodeData() {}

	public:
		virtual const Node * getNode() const = 0;
};

class TreeParsimonyData
{
	public:
		virtual unsigned int getArrayPosition(const Node * parent, const Node * son, unsigned int currentPosition) const = 0;
		virtual unsigned int getRootArrayPosition(const unsigned int site) const = 0;
		virtual TreeParsimonyNodeData & getNodeData(const Node * node) = 0;
		virtual const TreeParsimonyNodeData & getNodeData(const Node * node) const = 0;
};

class AbstractTreeParsimonyData
{
	protected:
		vector<unsigned int> _rootPatternLinks;
		vector<unsigned int> _rootWeights;

	public:
		unsigned int getRootArrayPosition(const unsigned int site) const { return _rootPatternLinks[site]; }
		unsigned int getWeight(unsigned int pos) const { return _rootWeights[pos]; }
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


