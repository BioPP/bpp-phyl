//
// File: TreeIterator.h
// Created by: Keren Halabi
// Created on: Thu Jul 5 14:03:18 2018
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

#ifndef _TREEITERATORS_H
#define _TREEITERATORS_H

#include "TreeTemplate.h"
#include "Node.h"

// From the STL:
#include <string>
#include <vector>
#include <map>

/**
 * @brief The phylogenetic tree iterator class.
 *
 * This class is part of the object implementation of phylogenetic trees.
 *
 * The class offers iterators for traversing the nodes in a tree in 3 possible orders:
 * Post-Order
 * Pre-Order
 * In-Order
 *
 * For more information on using trees in BIo++,
 * @see Node
 * @see NodeTemplate
 * @see TreeTools
 */

using namespace std;
namespace bpp
{

    class TreeIterator                               // abstract class from which each iterator type inherits
    {
        protected:
        const TreeTemplate<Node>& tree_;             // The tree to iterate over its nodes. Can't be const because user should be allowed to edit the nodes of the tree during the traversal.
        const Node* currNode_;                       // A pointer to the current node visited by the iterator. Can't be const because user should be allowed to edit the nodes of the tree during the traversal.
        map<int, bool> nodeToVisited_;               // a map that matches to each node id of booleans that for each node id states if it is visited or not
        map<int, bool> nodeToSonVisited_;            // a map that matches to each node id if a son of his was visited or not
        map<int, size_t> nodeToLastVisitedSonIndex_; // a map that matches to each node id the index of its last visited son visited son

        public:
        /* constructors and destructors */
        explicit TreeIterator(const TreeTemplate<Node>& tree):
            tree_(tree),
            currNode_(tree.getRootNode()), // The pointer to the initial node is initialized as 0 since it's actual assignment depends on which traversal is chosen
            nodeToVisited_(),
            nodeToSonVisited_(),
            nodeToLastVisitedSonIndex_()
            { init(); }

        explicit TreeIterator(const TreeIterator& tree_iterator):
        tree_(tree_iterator.tree_),
        currNode_(tree_.getRootNode()),
        nodeToVisited_(),
        nodeToSonVisited_(),
        nodeToLastVisitedSonIndex_()
        {}

        TreeIterator& operator=(const TreeIterator& tree_iterator);

        virtual ~TreeIterator() {}        // must be virtual to assume that upon deletion, the destructor of any inheriting class is called as well (see https://www.geeksforgeeks.org/virtual-destructor/)

        void init();
                                         // function to initialize nodes properties
        /* iterating functions */
        const Node* begin();
        virtual const Node* next() = 0;                     // Set as virtual because the function should be implemented separately in each iterator type
        TreeIterator& operator++();
        const Node* end(){ return NULL; }
    };


    class PostOrderTreeIterator: public TreeIterator
    {
        public:

        /* constructors and destructors */
        explicit PostOrderTreeIterator(const TreeTemplate<Node>& tree):
        TreeIterator(tree) {
            currNode_ = tree_.getNodes()[0]; // Get the leftmost leaf of the tree
        }

        ~PostOrderTreeIterator() {};           // Inherited from TreeIterator

        const Node* getLeftMostPredecessor(const Node* startNode);
        const Node* next();
    };


    class PreOrderTreeIterator: public TreeIterator
    {
        public:

        /* constructors and destructors */
        explicit PreOrderTreeIterator(const TreeTemplate<Node>& tree):
        TreeIterator(tree) {
            currNode_ = tree_.getRootNode();
        }

        ~PreOrderTreeIterator() {};           // Inherited from TreeIterator

        const Node* next();
    };


    class InOrderTreeIterator: public TreeIterator
    {
        public:

        /* constrcutors and destrcutors */
        explicit InOrderTreeIterator(const TreeTemplate<Node>& tree):
        TreeIterator(tree) {
            currNode_ = tree_.getNodes()[0];  // Get the leftmost leaf of the tree
        }

        ~InOrderTreeIterator() {};           // Inherited from TreeIterator

        const Node* doStep(const Node* node);
        const Node* next();
    };


} // end of namespace bpp.

#endif // _TREEITERATORS_H
