//
// File: TreeExceptions.cpp
// Created by: Julien Dutheil
// Created on: Mon Nov  3 17:04:46 2003
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

Julien.Dutheil@univ-montp2.fr

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

#include "TreeExceptions.h"
#include "Node.h"
#include "Tree.h"

#include <Utils/TextTools.h>

using namespace bpp;

/******************************************************************************/

NodeException::NodeException(const string & text, const Node * node) :
	Exception("NodeException: " + text + (node != NULL ? "(id:" + TextTools::toString(node->getId()) + ")" : "")),
	_node(node), _nodeId(node->getId()) {};
NodeException::NodeException(const string & text, int nodeId) :
	Exception("NodeException: " + text + "(id:" + TextTools::toString(nodeId) + ")"),
	_node(NULL), _nodeId(nodeId) {};
		
const Node * NodeException::getNode() const { return _node; }
int NodeException::getNodeId() const { return _nodeId; }

/******************************************************************************/

NodeNotFoundException::NodeNotFoundException(const string & text, const string & id) :
	Exception("NodeNotFoundException: " + text + "(" + id + ")"), _id(id) {}

NodeNotFoundException::NodeNotFoundException(const string & text, int id) :
	Exception("NodeNotFoundException: " + text + "(" + TextTools::toString(id) + ")")
{
  _id = TextTools::toString(id);
}
			
/******************************************************************************/

TreeException::TreeException(const string & text, const Tree * tree) :
			Exception("TreeException: " + text + (tree != NULL ? "(" + tree -> getName() + ")" : "")),
			_tree(tree) {};

const Tree * TreeException::getTree() const { return _tree; }

/******************************************************************************/

UnrootedTreeException::UnrootedTreeException(const string & text, const Tree * tree) :
			TreeException("UnrootedTreeException: " + text, tree) {}

/******************************************************************************/

