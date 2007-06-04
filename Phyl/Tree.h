//
// File: Tree.h
// Created by: Julien Dutheil
// Created on: Thu Mar 13 12:03:18 2003
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

#ifndef _TREE_H_
#define _TREE_H_

#include "TreeExceptions.h"

// From Utils:
#include <Utils/Clonable.h>

// From the STL:
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

using namespace std;

/**
 * @brief Interface for phylogenetic tree objects.
 * 
 */
class Tree:
  public Clonable
{

	public: // Constructors and destructor:
		
		Tree() {}
		virtual ~Tree() {}

    Tree * clone() const = 0;

	public:
		
		/**
		 * @brief Tree name.
		 *
		 * @{
		 */
		virtual string getName() const = 0;
		
		virtual void setName(const string & name) = 0;
		/** @} */
		
		virtual unsigned int getNumberOfLeaves() const = 0;
		
		virtual unsigned int getNumberOfNodes() const = 0;

		virtual vector<double> getBranchLengths() const = 0;

		virtual vector<string> getLeavesNames() const = 0;

		/**
		 * @name Retrieving ids.
		 *
		 * @{
		 */
		virtual int getRootId() const = 0;

		virtual int getLeafId(const string & name) const throw (NodeNotFoundException)= 0;
	
		virtual vector<int> getLeavesId() const = 0;

		virtual vector<int> getNodesId() const = 0;

		virtual vector<int> getSonsId(int parentId) const throw (NodeNotFoundException) = 0;

		virtual int getFatherId(int parentId) const throw (NodeNotFoundException) = 0;

		virtual bool hasFather(int nodeId) const throw (NodeNotFoundException) = 0;
		/** @} */

		/**
		 * @name Dealing with node names.
		 *
		 * @{
		 */
		virtual string getNodeName(int nodeId) const throw (NodeNotFoundException) = 0;
		
		virtual void setNodeName(int nodeId, const string & name) throw (NodeNotFoundException) = 0;
		
		virtual void deleteNodeName(int nodeId) throw (NodeNotFoundException) = 0;
		
		virtual bool hasNodeName(int nodeId) const throw (NodeNotFoundException) = 0;
		/** @} */
		
		/**
		 * @name Several tests.
		 *
		 * @{
		 */
    virtual bool hasNode(int nodeId) const = 0;

		virtual bool isLeaf(int nodeId) const throw (NodeNotFoundException) = 0;

		virtual bool isRoot(int nodeId) const throw (NodeNotFoundException) = 0;
		/** @} */

    /**
     * @name Acting on topology.
     *
     * @{
     */

    /**
     * @brief Swap two son nodes.
     *
     * @param tree The tree.
     * @param nodeId The node.
     * @param i1 First son node index.
     * @param i2 Second son node index.
     * @throw NodeNotFoundException If the node is not found.
     * @throw IndexOutOfBoundsException If one node index is not valid, or if the node
     */
    void swapNodes(const Tree & tree, int nodeId, unsigned int i1=0, unsigned int i2=1) throw (NodeNotFoundException,IndexOutOfBoundsException);
  
    /** @} */

		/**
		 * @name Dealing with branch lengths.
		 *
		 * @{
		 */
		virtual double getDistanceToFather(int nodeId) const throw (NodeNotFoundException) = 0;
		
		virtual void setDistanceToFather(int nodeId, double length) throw (NodeNotFoundException) = 0;
		
		virtual void deleteDistanceToFather(int nodeId) throw (NodeNotFoundException) = 0;
		
		virtual bool hasDistanceToFather(int nodeId) const throw (NodeNotFoundException) = 0;
		/** @} */

		/**
		 * @name Node properties.
		 *
		 * @{
		 */
		virtual bool hasNodeProperty(int nodeId, const string & name) const throw (NodeNotFoundException) = 0;
		
		virtual void setNodeProperty(int nodeId, const string & name, const Clonable & property) throw (NodeNotFoundException) = 0;
				
		virtual Clonable * getNodeProperty(int nodeId, const string & name) throw (NodeNotFoundException) = 0;
				
		virtual const Clonable * getNodeProperty(int nodeId, const string & name) const throw (NodeNotFoundException) = 0;
				
		virtual Clonable * removeNodeProperty(int nodeId, const string & name) throw (NodeNotFoundException) = 0;

    virtual vector<string> getNodePropertyNames(int nodeId) const throw (NodeNotFoundException) = 0;
		/** @} */
		
		/**
		 * @name Branch properties.
		 *
		 * @{
		 */
		virtual bool hasBranchProperty(int nodeId, const string & name) const throw (NodeNotFoundException) = 0;
		
		virtual void setBranchProperty(int nodeId, const string & name, const Clonable & property) throw (NodeNotFoundException) = 0;
				
		virtual Clonable * getBranchProperty(int nodeId, const string & name) throw (NodeNotFoundException) = 0;
				
		virtual const Clonable * getBranchProperty(int nodeId, const string & name) const throw (NodeNotFoundException) = 0;
				
		virtual Clonable * removeBranchProperty(int nodeId, const string & name) throw (NodeNotFoundException) = 0;

    virtual vector<string> getBranchPropertyNames(int nodeId) const throw (NodeNotFoundException) = 0;
		/** @} */

		/**
		 * @brief Change the root node.
		 *
		 * Works on unrooted tree.
		 * If the tree is rooted, the method unroots it first.
		 *
		 * @param nodeId The id of the node that will be the new root.
		 */
		virtual void rootAt(int nodeId) throw (NodeNotFoundException) = 0;

		/**
		 * @brief Root a tree by specifying an outgroup.
		 *
		 * If the tree is rooted, unroot it first, change the root node and then
     * reroot the tree using the previous root id.
		 * If the tree is unrooted, change the root node and then create a new root node.
		 *
		 * @param nodeId The id of the node that will be the new root.
		 */
		virtual void newOutGroup(int nodeId) throw (NodeNotFoundException) = 0;
		
		/**
		 * @brief Tell if the tree is rooted.
		 * 
		 * @return True if the tree is rooted.
		 */
		virtual bool isRooted() const = 0;
		
		/**
		 * @brief Unroot a rooted tree.
		 *
		 * @return True if the tree has been unrooted.
		 * @throw UnrootedTreeException If the tree is already rooted.
		 */
		virtual bool unroot() throw (UnrootedTreeException) = 0;

		/**
		 * @brief Number nodes.
		 */
		virtual void resetNodesId() = 0;
		
		// Works on (multi)furcations:
		
		/**
		 * @brief Tell if the tree is multifurcating.
		 * 
		 * @return True if the tree is multifurcating.
		 */
		virtual bool isMultifurcating() const = 0;
		
		/**
		 * @brief Get all the branch lengths of a tree.
		 *
		 * @return A vector with all branch lengths.
		 * @throw NodeException If a branch length is lacking.
		 */
		virtual vector<double> getBranchLengths() throw (NodeException) = 0;

		/**
		 * @brief Get the total length (sum of all branch lengths) of a tree.
		 *
		 * @return The total length of the subtree.
		 * @throw NodeException If a branch length is lacking.
		 */
		virtual double getTotalLength() throw (NodeException) = 0;

		/**
		 * @brief Set all the branch lengths of a tree.
		 *
		 * @param brLen The branch length to apply.
		 */
		virtual void setBranchLengths(double brLen) = 0;
		
		/**
		 * @brief Give a length to branches that don't have one in a tree.
		 *
		 * @param brLen The branch length to apply.
		 */
		virtual void setVoidBranchLengths(double brLen) = 0;
	
		/**
		 * @brief Scale a given tree.
		 *
		 * Multiply all branch lengths by a given factor.
		 *
		 * @param factor The factor to multiply all branch lengths with.
		 * @throw NodeException If a branch length is lacking.
		 */
		virtual void scaleTree(double factor) throw (NodeException) = 0;

    /**
     * @brief Get an id.
     *
     * @return an unused node id.
     */
    virtual int getNextId() = 0;

};

#endif	//_TREE_H_

