//
// File: Tree.h
// Created by: jdutheil <julien.dutheil@ens-lyon.fr>
// Created on: Thu Mar 13 12:03:18 2003
//

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 


#ifndef _TREE_H_
#define _TREE_H_

#include "TreeExceptions.h"

// From the STL:
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

// From Utils:
#include <Utils/Clonable.h>

using namespace std;

/**
 * @brief <p>The phylogenetic node class.</p>
 * <p>This class is part of the object implementation of phylogenetic trees. Tree are made
 * made of nodes, instances of this class.</p>
 * <p>Since trees are oriented (rooted), each node has one <i>father node</i> and possibly
 * many <i>son nodes</i>. Leaves are nodes without descendant and root is defined has the without
 * father. Inner nodes will generally contain two descendants (the tree is then called
 * <i>bifurcating</i>), but mutlifurcating trees are also allowed with this kind of description.
 * In the rooted case, each inner node also defines a <i>subtree</i>.
 * This allows to work recursively on trees, which is very convenient in most cases.</p>
 * <p>This class is made the more general as possible, while keeping it very simple. It contains:</p>
 * <ul>
 * <li>An identity tag, to identity it in the tree;</li>
 * <li>A name, necessary for leaf nodes, optionnal else;</li>
 * <li>A pointer toward the father node;</li>
 * <li>A vector of pointer toward son nodes;</li>
 * <li>The distance from the father node:</li>
 * <li>A property map, that may contain any information to link to each node, e.g. bootastrap
 * value or GC content.</li>
 * </ul>
 * @see Tree
 */

class Node {
			
	protected:
		
		int                             _id;
		string                        * _name;
		vector<Node *>                  _sons;
		Node                          * _father;
		double                        * _distanceToFather;
		mutable map<string, Clonable *> _properties;

	public: // Class constructor and destructor:
	
		/**	
		 * @brief Build a new void Node object.
		 */
		Node();
			
		/**
		 * @brief Build a new Node with specified id.
		 */
		Node(int id);

		/**
		 * @brief Build a new Node with specified name.
		 */
		Node(const string & name);

		/**
		 * @brief Build a new Node with specified id and name.
		 */
		Node(int id, const string & name);

		/**
		 * @brief Copy constructor.
		 * 
		 * @param node The node to copy.
		 */
		Node(const Node & node);

		/**
		 * @brief Assignation operator.
		 *
		 * @param node the node to copy.
		 * @return A reference toward this node.
		 */
		Node & operator=(const Node & node);

		virtual ~Node();

	public:				
				
		// Identity:
			
		/**
		 * @brief Get the node's id.
		 *
		 * @return The identity tag of this node.
		 */
		virtual int getId() const;
		
		/**
		 * @brief Set this node's id.
		 *
		 * @param id The new identity tag.
		 */
		virtual void setId(int id);

		// Name:
			
		/**
		 * @brief Get the name associated to this node, if there is one, 
		 * otherwise throw a NodeException.
		 * @return The name associated to this node.
		 */
		virtual string getName() const throw (NodeException);
				
		/**
		 * @brief <p>Give a name or update the name associated to the node.</p>
		 * @param name The name to give to the node.
		 */
		virtual void setName(string & name);
		
		/**
		 * @brief <p>Delete the name associated to this node (do nothing if
		 * there is no name).</p>
		 */
		virtual void deleteName();
		
		/**
		 * @brief <p>Tell is this node has a name.</p>
		 * @return True if name != NULL.
		 */
		virtual bool hasName() const;
				
		// Distances:
		
		/**
		 * @brief <p>Get the distance to the father node is there is one,
		 * otherwise throw a NodeException.</p>
		 * @return The distance to the father node.
		 */		 
		virtual double getDistanceToFather() const throw (NodeException);
				
		/**
		 * @brief <p>Set or update the distance toward the father node.</p>
		 * @param distance The new distance to the father node.
		 */
		virtual void setDistanceToFather(double distance);
		
		/**
		 * @brief <p>Delete the distance to the father node.</p>
		 */
		virtual void deleteDistanceToFather();
				
		/**
		 * @brief <p>Tell is this node has a distance to the father.</p>
		 * @return True if distanceToFather != NULL.
		 */
		virtual bool hasDistanceToFather() const;
				
		// Father:
			
		/**
		 * @brief <p>Get the father of this node is there is one.</p>
		 * @return A pointer toward the father node, NULL if there is not.
		 */		 
		virtual const Node * getFather() const;
				
		/**
		 * @brief <p>Get the father of this node is there is one.</p>
		 * @return A pointer toward the father node, NULL if there is not.
		 */		 
		virtual Node * getFather();
				
		/**
		 * @brief <p>Set the father node of this node.</p>
		 * @param node The father node.
		 */
		virtual void setFather(Node & node);
				
		/**
		 * @brief <p>Remove the father of this node.</p>
		 */
		virtual Node * removeFather();
				
		/**
		 * @brief <p>Tell if this node has a father node.</p>
		 */
		virtual bool hasFather() const;

		// Sons:
				 
		unsigned int getNumberOfSons() const;

		virtual const Node * getSon(unsigned int i) const;
				
		virtual Node * getSon(unsigned int i);
				
		virtual void addSon(unsigned int pos, Node & node);

		virtual void addSon(Node & node);
				
		virtual void removeSon(unsigned int i);
				
		virtual void removeSon(Node & node);
		
		virtual void removeSons();
				
		virtual void swap(unsigned int branch1, unsigned int branch2);
		
		/*
		 * Operators:
		 * - a positive value send the corresponding son
		 * - a negative value send the father.
		 */
				 
		Node * operator[](int i);
				
		const Node * operator[](int i) const;
				
		// Properties:
				
		virtual void setProperty(const string & name, Clonable * property);
				
		virtual Clonable * getProperty(const string & name);
				
		virtual const Clonable * getProperty(const string & name) const;
				
		virtual void * removeProperty(const string & name);
				
		virtual bool hasProperty(const string & name) const;
				
		// Equality operator:

		virtual bool operator==(const Node & node) const;
				
		// Tests:

		bool isLeaf() const;
			
	friend class Tree;
};


/**
 * @brief <p>The phylogenetic tree class.</p>
 * <p>This class is part of the object implementation of phylogenetic trees. Tree are made
 * made of nodes, instances of the class Node.</p>
 * <p>Trees are oriented (rooted), i.e. each node has one <i>father node</i> and possibly
 * many <i>son nodes</i>. Leaves are nodes without descendant and root is defined has the without
 * father. Inner nodes will generally contain two descendants (the tree is then called
 * <i>bifurcating</i>), but mutlifurcating trees are also allowed with this kind of description.
 * In the rooted case, each inner node also defines a <i>subtree</i>.
 * This allows to work recursively on trees, which is very convenient in most cases.
 * To deal with non-rooted trees, we place an artificial root at a particular node:
 * hence the root node appears to be trifurcated. This is the way unrooted trees are
 * described in the parenthetic description, the so called Newick format.</p>
 * @see Node
 */
class Tree {

	/**
	 * Fields:
	 */
	protected:
		Node * _root;
		string _name;

	public: // Constructors and destructor:
		
		Tree(); 

		Tree(const Tree & t);

		Tree(Node & root); 

		Tree& operator=(const Tree & t);

		virtual ~Tree();
			
	/**
	 * Methods:
 	 */
	
	public:
		Node * getRootNode();

		const Node * getRootNode() const;

		void setRootNode(Node & root);
	
		string getName() const;
	
		void setName(const string & name);

		unsigned int getNumberOfLeaves() const;

		vector<const Node *> getLeaves() const;

		vector<Node *> getLeaves();

		vector<Node *> getNodes();

		vector<const Node *> getNodes() const;

		vector<double> getBranchLengths() const;

		vector<string> getLeavesNames() const;

		//void setNewOutgroup(Tree::Node & node);

		//The foolowing methods are adapted from the tree class of the SEMPHY library:
		vector<Node *> getPathBetweenAnyTwoNodes(Node & node1, Node & node2) const;

		// Works on root:
		
		/**
		 * @brief <p></p>
		 * @param p_iNewRoot The node to be considered as the new root, i.e.
		 * defining the subtree to be considered as the new outgroup.
		 */
		void rootAt(Node & p_iNewRoot);
		
		/**
		 * @brief <p>Tell if the tree is rooted.</p>
		 * @return True if the tree is rooted.
		 */
		bool isRooted() const;
		
		/**
		 * @brief Unroot a rooted tree.
		 *
		 * @return True if the tree has been unrooted.
		 * @throw UnrootedTreeException If the tree is already rooted.
		 */
		bool unroot() throw (UnrootedTreeException);
		
		/**
		 * @brief Number nodes.
		 */
		void resetNodesId();
		
		// Works on (multi)furcations:
		
		/**
		 * @brief <p>Tell if the tree is multifurcating.</p>
		 * @return True if the tree is multifurcating.
		 */
		bool isMultifurcating() const;
		
		
	protected:
		
		void destroyNode(const Node & node);
	
		Node * cloneSubtree(const Node & node);
		
};

#endif	//_TREE_H_
