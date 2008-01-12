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

#ifndef _ABSTRACTAGGLOMERATIVEDISTANCEMETHOD_H_
#define _ABSTRACTAGGLOMERATIVEDISTANCEMETHOD_H_

#include "AgglomerativeDistanceMethod.h"
#include "Node.h"
#include "TreeTemplate.h"

// From the STL:
#include <map>

using namespace std;

namespace bpp
{

/**
 * @brief Partial implementation of the AgglomerativeDistanceMethod interface.
 *
 * This class provides a DistanceMatrix object for computations, and a map
 * with pivot indices and a pointer toward the corresponding subtree.
 *
 * Several methods, commons to several algorithm are provided.
 */
class AbstractAgglomerativeDistanceMethod:
  public virtual AgglomerativeDistanceMethod
{
	protected:
		DistanceMatrix _matrix;
		Tree * _tree;

 		map<unsigned int, Node *> _currentNodes;
    bool _verbose;
	
	public:
		AbstractAgglomerativeDistanceMethod(): _matrix(0), _tree(NULL), _currentNodes(), _verbose(true) {}
		AbstractAgglomerativeDistanceMethod(bool verbose): _matrix(0), _tree(NULL), _currentNodes(), _verbose(verbose) {}
		AbstractAgglomerativeDistanceMethod(const DistanceMatrix & matrix, bool verbose = true): _matrix(0), _tree(NULL), _currentNodes(), _verbose(verbose)
    {
      setDistanceMatrix(matrix);
    }
		virtual ~AbstractAgglomerativeDistanceMethod()
    {
      delete _tree;
    }
    AbstractAgglomerativeDistanceMethod(const AbstractAgglomerativeDistanceMethod & a): _matrix(a._matrix), _tree(NULL), _currentNodes()
    {
      // Hard copy of inner tree:
      if(a._tree != NULL) _tree = new TreeTemplate<Node>(* a._tree);
    }
    AbstractAgglomerativeDistanceMethod & operator=(const AbstractAgglomerativeDistanceMethod & a)
    {
      _matrix = a._matrix;
      // Hard copy of inner tree:
      if(a._tree != NULL) _tree = new TreeTemplate<Node>(* a._tree);
      else _tree = NULL;
      return *this;
    }

	public:
		virtual void setDistanceMatrix(const DistanceMatrix & matrix);

    /**
     * @brief Get the computed tree, if there is one.
     *
     * @return A copy of the computed tree if there is one, NULL otherwise.
     */
    virtual
#if defined(NO_VIRTUAL_COV)
		Tree *
#else
		TreeTemplate<Node> * 
#endif
		getTree() const
    {
    	//Node * root = TreeTools::cloneSubtree<Node>(* dynamic_cast<TreeTemplate<Node> *>(_tree) -> getRootNode());
	    //return new TreeTemplate<Node>(* root);
      return _tree == NULL ? NULL : new TreeTemplate<Node>(*_tree);
    }
		
    /**
     * @brief Compute the tree corresponding to the distance matrix.
     *
     * This method implements the following algorithm:
     * 1) Build all leaf nodes (getLeafNode method)
     * 2) Get the best pair to agglomerate (getBestPair method)
     * 3) Compute the branch lengths for this pair (computeBranchLengthsForPair method)
     * 4) Build the parent node of the pair (getParentNode method)
     * 5) For each remaining node, update distances from the pair (computeDistancesFromPair method)
     * 6) Return to step 2 while there are more than 3 remaining nodes.
     * 7) Perform the final step, and send a rooted or unrooted tree.
     * 
     * @param rooted Tell if the final tree must be rooted or not.
     */
		virtual void computeTree(bool rooted) throw (Exception);

    void setVerbose(bool yn) { _verbose = yn; }
    bool isVerbose() const { return _verbose; }

	protected:
    /**
     * @name Specific methods.
     *
     * @{
     */

    /**
     * @brief Get the best pair of nodes to agglomerate.
     *
     * Define the criterion to chose the next pair of nodes to agglomerate.
     * This criterion uses the _matrix distance matrix.
     *
     * @return A size 2 vector with the indices of the nodes.
     * @throw Exception If an error occured.
     */
		virtual vector<unsigned int> getBestPair() throw (Exception) = 0;
		
    /**
     * @brief Compute the branch lengths for two nodes to agglomerate.
     *
     * @code
     * +---l1-----N1
     * |
     * +---l2-----N2
     * @endcode
     * This method compute l1 and l2 given N1 and N2.
     *
     * @param pair The indices of the nodes to be agglomerated.
     * @return A size 2 vector with branch lengths.
     */
    virtual vector<double> computeBranchLengthsForPair(const vector<unsigned int> & pair) = 0;

    /**
     * @brief Actualizes the distance matrix according to a given pair and the corresponding branch lengths.
     *
     * @param pair The indices of the nodes to be agglomerated.
     * @param branchLengths The corresponding branch lengths.
     * @param pos The index of the node whose distance ust be updated.
     * @return The distance between the 'pos' node and the agglomerated pair.
     */
		virtual double computeDistancesFromPair(const vector<unsigned int> & pair, const vector<double> & branchLengths, unsigned int pos) = 0;
		
    /**
     * @brief Method called when there ar eonly three remaining node to agglomerate, and creates the root node of the tree.
     *
     * @param idRoot The id of the root node.
     */
    virtual void finalStep(int idRoot) = 0;

    /**
     * @brief Get a leaf node.
     *
     * Create a new node with the given id and name.
     *
     * @param id The id of the node.
     * @param name The name of the node.
     * @return A pointer toward a new node object.
     */
		virtual Node * getLeafNode(int id, const string & name);

    /**
     * @brief Get an inner node.
     *
     * Create a new node with the given id, and set its sons.
     *
     * @param id The id of the node.
     * @param son1 The first son of the node.
     * @param son2 The second son of the node.
     * @return A pointer toward a new node object.
     */
		virtual Node * getParentNode(int id, Node * son1, Node * son2);
    /** @} */
		
};

} //end of namespace bpp.

#endif //_ABSTRACTAGGLOMERATIVEDISTANCEMETHOD_H_

