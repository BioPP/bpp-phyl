//
// File: Tree.h
// Created by: jdutheil <julien.dutheil@ens-lyon.fr>
// Created on: Thu Mar 13 12:03:18 2003
//

/*
Copyright ou © ou Copr. Julien Dutheil, (16 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant à fournir des classes
pour l'analyse de données phylogénétiques.

Ce logiciel est régi par la licence CeCILL soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accepté les
termes.
*/

/*
Copyright or © or Copr. Julien Dutheil, (November 16, 2004)

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

#ifndef _TREE_H_
#define _TREE_H_

#include "TreeExceptions.h"
#include "TreeTools.h"

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
		Node() : _id(0), _name(NULL), _father(NULL), _distanceToFather(NULL) {}
			
		/**
		 * @brief Build a new Node with specified id.
		 */
		Node(int id) : _id(id), _name(NULL), _father(NULL), _distanceToFather(NULL) {}

		/**
		 * @brief Build a new Node with specified name.
		 */
		Node(const string & name) : _id(0), _name(new string(name)), _father(NULL), _distanceToFather(NULL) {}

		/**
		 * @brief Build a new Node with specified id and name.
		 */
		Node(int id, const string & name) : _id(id), _name(new string(name)), _father(NULL), _distanceToFather(NULL) {}

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

		virtual ~Node() {	delete _name;	delete _distanceToFather;	}

	public:				
				
		/**
		 * @name Identity
		 *
		 * @{
		 */
			
		/**
		 * @brief Get the node's id.
		 *
		 * @return The identity tag of this node.
		 */
		virtual int getId() const { return _id; }
		
		/**
		 * @brief Set this node's id.
		 *
		 * @param id The new identity tag.
		 */
		virtual void setId(int id) { _id = id; }

		/** @} */

		/**
		 * @name Name:
		 *
		 * @{
		 */
			
		/**
		 * @brief Get the name associated to this node, if there is one, 
		 * otherwise throw a NodeException.
		 * 
		 * @return The name associated to this node.
		 */
		virtual string getName() const throw (NodeException)
		{
			if(!hasName()) throw NodeException("Node::getName: no name associated to this node.", this);
			return * _name;
		}
				
		/**
		 * @brief Give a name or update the name associated to the node.
		 * 
		 * @param name The name to give to the node.
		 */
		virtual void setName(string & name) { delete _name; _name = new string(name); }
		
		/**
		 * @brief Delete the name associated to this node (do nothing if there is no name).
		 */
		virtual void deleteName() { delete _name; _name = NULL; }
		
		/**
		 * @brief Tell is this node has a name.
		 * 
		 * @return True if name != NULL.
		 */
		virtual bool hasName() const { return _name != NULL; }
		
		/** @} */
		
		/**
		 * @name Distances:
		 *
		 * @{
		 */
		
		/**
		 * @brief Get the distance to the father node is there is one,
		 * otherwise throw a NodeException.
		 * 
		 * @return The distance to the father node.
		 */		 
		virtual double getDistanceToFather() const throw (NodeException) 
		{
			if(!hasDistanceToFather()) throw NodeException("Node::getDistanceToFather: Node has no distance." , this);
			return * _distanceToFather;
		}
				
		/**
		 * @brief Set or update the distance toward the father node.
		 * 
		 * @param distance The new distance to the father node.
		 */
		virtual void setDistanceToFather(double distance) { delete _distanceToFather; _distanceToFather = new double(distance); }
		
		/**
		 * @brief Delete the distance to the father node.
		 */
		virtual void deleteDistanceToFather() { delete _distanceToFather; _distanceToFather = NULL; }
				
		/**
		 * @brief Tell is this node has a distance to the father.
		 * 
		 * @return True if distanceToFather != NULL.
		 */
		virtual bool hasDistanceToFather() const { return _distanceToFather != NULL; }

		/** @} */
				
		/**
		 * @name Father:
		 *
		 * @{
		 */
			
		/**
		 * @brief Get the father of this node is there is one.
		 * 
		 * @return A pointer toward the father node, NULL if there is not.
		 */		 
		virtual const Node * getFather() const { return _father; }
				
		/**
		 * @brief Get the father of this node is there is one.
		 * 
		 * @return A pointer toward the father node, NULL if there is not.
		 */		 
		virtual Node * getFather() { return _father; }
				
		/**
		 * @brief Set the father node of this node.
		 * 
		 * @param node The father node.
		 */
		virtual void setFather(Node & node) { _father = & node; }
				
		/**
		 * @brief Remove the father of this node.
		 */
		virtual Node * removeFather() { Node * f = _father; _father = NULL; return f; }
				
		/**
		 * @brief Tell if this node has a father node.
		 */
		virtual bool hasFather() const { return _father != NULL; }

		/** @} */

		/**
		 * @name Sons:
		 *
		 * @{
		 */
				 
		virtual unsigned int getNumberOfSons() const { return _sons.size(); }

		virtual const Node * getSon(unsigned int i) const { return _sons[i]; }
				
		virtual Node * getSon(unsigned int i) { return _sons[i]; }
				
		virtual void addSon(unsigned int pos, Node & node)
		{
			_sons.insert(_sons.begin() + pos, & node);
			node._father = this;
		}

		virtual void addSon(Node & node)
		{
			_sons.push_back(& node);
			node._father = this;
		}

		virtual void setSon(unsigned int pos, Node & node) throw (NodeException)
		{
			if(pos >= _sons.size()) throw NodeException("Node::setSon(). Invalid node position.", & node);
			_sons[pos] = & node;
			node._father = this;
		}
				
		virtual void removeSon(unsigned int i)
	  {
			Node * node = _sons[i];
			_sons.erase(_sons.begin() + i);
			node -> removeFather();
		}
		
		virtual void removeSon(Node & node)
		{
			for(unsigned int i = 0; i < _sons.size(); i++) if(* _sons[i] == node) _sons.erase(_sons.begin() + i);
			node.removeFather();
		}
		
		virtual void removeSons() {	while(_sons.size() != 0) removeSon(0); }
				
		virtual void swap(unsigned int branch1, unsigned int branch2);

		/** @} */

		virtual vector<const Node *> getNeighbors() const;
		
		virtual vector<Node *> getNeighbors();
		
		/**
		 * @name Operators:
		 * - a positive value send the corresponding son
		 * - a negative value send the father.
		 *
		 * @{
		 */
				 
		Node * operator[](int i) { return (i < 0) ? _father : _sons[i]; }
				
		const Node * operator[](int i) const { return (i < 0) ? _father : _sons[i]; }
		
		/** @} */
		
		/**
		 * @name Properties:
		 *
		 * @{
		 */
				
		virtual void setProperty(const string & name, Clonable * property) { _properties[name] = property; }
				
		virtual Clonable * getProperty(const string & name) { return _properties[name]; }
				
		virtual const Clonable * getProperty(const string & name) const { return const_cast<const Clonable *>(_properties[name]); }
				
		virtual void * removeProperty(const string & name)
		{
			void * removed = _properties[name];
			_properties.erase(name);
			return removed;
		}	
				
		virtual bool hasProperty(const string & name) const { return _properties.find(name) != _properties.end(); }
		
		/** @} */
		
		// Equality operator:

		virtual bool operator==(const Node & node) const { return _id == node._id; }	
				
		// Tests:

		bool isLeaf() const { return _sons.size() == 0; }
			
	friend class Tree<Node>;
};

template<class NodeInfos>
class NodeTemplate : public Node {
	
	protected:

		NodeInfos * _infos;

	public:
		/**	
		 * @brief Build a new void NodeTemplate object.
		 */
		NodeTemplate() : Node() {}
			
		/**
		 * @brief Build a new NodeTemplate with specified id.
		 */
		NodeTemplate(int id) : Node(id) {}

		/**
		 * @brief Build a new NodeTemplate with specified name.
		 */
		NodeTemplate(const string & name) : Node(name) {}

		/**
		 * @brief Build a new NodeTemplate with specified id and name.
		 */
		NodeTemplate(int id, const string & name) : Node(id, name) {}

		/**
		 * @brief Copy constructor.
		 * 
		 * @param node The node to copy.
		 */
		NodeTemplate(const NodeTemplate<NodeInfos> & node)
		{
			_id               = node._id;
			_name             = node.hasName() ? new string(* node._name) : NULL;
			_father           = node._father;
			_distanceToFather = node.hasDistanceToFather() ? new double(* node._distanceToFather) : NULL;
			_sons             = node._sons;
			for(map<string, Clonable *>::iterator i = node._properties.begin(); i != node._properties.end(); i++)
				_properties[i -> first] = i -> second -> clone();
			_infos            = node._infos;
		}

		/**
		 * @brief Copy constructor with a different template.
		 * 
		 * @param node The node to copy.
		 */
		template<class AnotherNodeInfos>
		NodeTemplate(const NodeTemplate<AnotherNodeInfos> & node)
		{
			_id               = node._id;
			_name             = node.hasName() ? new string(* node._name) : NULL;
			_father           = node._father;
			_distanceToFather = node.hasDistanceToFather() ? new double(* node._distanceToFather) : NULL;
			_sons             = node._sons;
			for(map<string, Clonable *>::iterator i = node._properties.begin(); i != node._properties.end(); i++)
				_properties[i -> first] = i -> second -> clone();
			_infos            = node._infos; // This operator must be defined!
		}

		/**
		 * @brief Assignation operator.
		 *
		 * @param node the node to copy.
		 * @return A reference toward this node.
		 */
		NodeTemplate<NodeInfos> & operator=(const NodeTemplate<NodeInfos> & node)
		{
			_id               = node._id;
			_name             = node.hasName() ? new string(* node._name) : NULL;
			_father           = node._father;
			_distanceToFather = node.hasDistanceToFather() ? new double(* node._distanceToFather) : NULL;
			_sons             = node._sons;
			_infos            = node._infos;
			for(map<string, Clonable *>::iterator i = node._properties.begin(); i != node._properties.end(); i++)
				_properties[i -> first] = i -> second -> clone();
			return * this;
		}

		virtual ~NodeTemplate() {}

};


/**
 * @brief The phylogenetic tree class.
 * 
 * This class is part of the object implementation of phylogenetic trees. Tree are made
 * made of nodes, instances of the class Node.
 * 
 * Trees are oriented (rooted), i.e. each node has one <i>father node</i> and possibly
 * many <i>son nodes</i>. Leaves are nodes without descendant and root is defined has the without
 * father. Inner nodes will generally contain two descendants (the tree is then called
 * <i>bifurcating</i>), but mutlifurcating trees are also allowed with this kind of description.
 * In the rooted case, each inner node also defines a <i>subtree</i>.
 * This allows to work recursively on trees, which is very convenient in most cases.
 * To deal with non-rooted trees, we place an artificial root at a particular node:
 * hence the root node appears to be trifurcated. This is the way unrooted trees are
 * described in the parenthetic description, the so called Newick format.
 * 
 * @see Node
 */
template<class N=Node>
class Tree {

	/**
	 * Fields:
	 */
	protected:
		N * _root;
		string _name;

	public: // Constructors and destructor:
		
		Tree() { _root = NULL; }

		Tree(const Tree<N> & t)
		{
			//Perform a hard copy of the nodes:
			_root = TreeTools::cloneSubtree<N>(* t.getRootNode());
		}

		template<class AnotherNodeType>
		Tree(const Tree<AnotherNodeType> & t)
		{
			//Perform a hard copy of the nodes:
			_root = TreeTools::cloneSubtree<AnotherNodeType>(* t.getRootNode());
		}

		Tree(N & root) { _root = &root; }

		Tree<N> & operator=(const Tree<N> & t)
		{
			//Perform a hard copy of the nodes:
			_root = TreeTools::cloneSubtree<N>(* t.getRootNode());
    	return *this;
		}

		virtual ~Tree() { destroyNode(* _root); delete _root; }

			
	/**
	 * Methods:
 	 */
	
	public:
		
		N * getRootNode() { return _root; }

		const N * getRootNode() const { return _root; }

		void setRootNode(N & root) { _root = & root; }
	
		string getName() const { return _name; }
	
		void setName(const string & name) { _name = name; }

		unsigned int getNumberOfLeaves() const { return TreeTools::getNumberOfLeaves(* const_cast<const N *>( _root)); }

		vector<const N *> getLeaves() const { return TreeTools::getLeaves(* const_cast<const N *>(_root)); }

		vector<      N *> getLeaves()       { return TreeTools::getLeaves(* _root); }


		vector<const N *> getNodes() const { return TreeTools::getNodes (* const_cast<const N *>(_root)); }

		vector<      N *> getNodes()       { return TreeTools::getNodes (* _root); }

		vector<double> getBranchLengths() const { return TreeTools::getBranchLengths(* _root); }

		vector<string> getLeavesNames() const { return TreeTools::getLeavesNames(* const_cast<const N *>( _root)); }

		//void setNewOutgroup(Tree::N & node);

		// Works on root:
		
		/**
		 * @brief <p></p>
		 * @param p_iNewRoot The node to be considered as the new root, i.e.
		 * defining the subtree to be considered as the new outgroup.
		 */
		void rootAt(N & p_iNewRoot)
	  {
			if (* _root == p_iNewRoot) return;
			vector<Node *> pathMatrix = TreeTools::getPathBetweenAnyTwoNodes(* _root, p_iNewRoot);
			//pathMatrix size is always bigger than 2.

			for (unsigned int i = 0; i < pathMatrix.size() - 1 ; i++) {
				pathMatrix[i] -> _father = pathMatrix[i + 1];
				pathMatrix[i] -> setDistanceToFather(pathMatrix[i + 1] -> getDistanceToFather());
				vector<Node *>::iterator vec_iter;
				vec_iter = remove(pathMatrix[i] -> _sons.begin(), pathMatrix[i] -> _sons.end(), pathMatrix[i + 1]);
				pathMatrix[i] -> _sons.erase(vec_iter, pathMatrix[i] -> _sons.end()); // pg 1170, primer.
	
				pathMatrix[i+1] -> _sons.push_back(pathMatrix[i + 1] -> getFather());
				pathMatrix[i+1] -> _father = NULL;
			}
			_root = & p_iNewRoot;
		}

		
		/**
		 * @brief Tell if the tree is rooted.
		 * 
		 * @return True if the tree is rooted.
		 */
		bool isRooted() const { return _root -> getNumberOfSons() == 2; }
		
		/**
		 * @brief Unroot a rooted tree.
		 *
		 * @return True if the tree has been unrooted.
		 * @throw UnrootedTreeException If the tree is already rooted.
		 */
		bool unroot() throw (UnrootedTreeException<N>)
		{
			if(!isRooted()) throw UnrootedTreeException<N>("Tree::unroot", this);

    		if(_root -> getNumberOfSons() == 2) {
				N* son1 = _root -> getSon(0);
				N* son2 = _root -> getSon(1);
				if(son1 -> isLeaf() && son2 -> isLeaf()) return false; // We can't unroot a single branch!
					// We manage to have a subtree in position 0:
					if(son1 -> isLeaf()) {
						_root -> swap(0, 1);
						son1 = _root -> getSon(0);
						son2 = _root -> getSon(1);
					}

					// Take care of branch lengths:
					if(son1 -> hasDistanceToFather()) {
						if(son2 -> hasDistanceToFather()) {
						// Both nodes have lengths, we sum them:
						son2 -> setDistanceToFather(son1 -> getDistanceToFather() + son2 -> getDistanceToFather());
					} else {
						// Only node 1 has length, we set it to node 2:
						son2 -> setDistanceToFather(son1 -> getDistanceToFather());
					}
					son1 -> deleteDistanceToFather();
				} // Else node 2 may or may not have a branch length, we do not care!

				// Remove the root:
				_root -> removeSons();
				son1 -> addSon(*son2);
				setRootNode(*son1);
				return true;
			} else return false; // Tree is already rooted.
		}


		
		/**
		 * @brief Number nodes.
		 */
		void resetNodesId()
		{
			vector<Node *> nodes = getNodes();
			for(unsigned int i = 0; i < nodes.size(); i++) nodes[i] -> setId(i);
		}
		
		// Works on (multi)furcations:
		
		/**
		 * @brief Tell if the tree is multifurcating.
		 * 
		 * @return True if the tree is multifurcating.
		 */
		bool isMultifurcating() const
		{
			bool b = false;
			for(unsigned int i = 0; i < _root -> getNumberOfSons(); i++) {
				b = b || TreeTools::isMultifurcating(* _root -> getSon(i));
			}
			return b;
		}
		
		/**
		 * @brief Get all the branch lengths of a tree.
		 *
		 * @param tree The tree.
		 * @return A vector with all branch lengths.
		 * @throw NodeException If a branch length is lacking.
		 */
		Vdouble getBranchLengths() throw (NodeException)
		{
			Vdouble brLen(1);
			for(unsigned int i = 0; i < _root -> getNumberOfSons(); i++) {
				Vdouble sonBrLen = getBranchLengths(* _root -> getSon(i));
				for(unsigned int j = 0; j < sonBrLen.size(); j++) brLen.push_back(sonBrLen[j]);
			}
			return brLen;
		}

		/**
		 * @brief Get the total length (sum of all branch lengths) of a tree.
		 *
		 * @param tree The tree.
		 * @return The total length of the subtree.
		 * @throw NodeException If a branch length is lacking.
		 */
		double getTotalLength() throw (NodeException)
		{
			return TreeTools::getTotalLength(*_root);
		}

		/**
		 * @brief Set all the branch lengths of a tree.
		 *
		 * @param node  The node.
		 * @param brLen The branch length to apply.
		 */
		void setBranchLengths(double brLen)
		{
			TreeTools::setBranchLengths(*_root, brLen);
		}
		
		/**
		 * @brief Give a length to branches that don't have one in a tree.
		 *
		 * @param node  The node.
		 * @param brLen The branch length to apply.
		 */
		void setVoidBranchLengths(double brLen)
		{
			TreeTools::setVoidBranchLengths(*_root, brLen);
		}
	
		/**
		 * @brief Scale a given tree.
		 *
		 * Multiply all branch lengths by a given factor.
		 *
		 * @param tree   The tree to scale.
		 * @param factor The factor to multiply all branch lengths with.
		 * @throw NodeException If a branch length is lacking.
		 */
		void scaleTree(double factor) throw (NodeException)
		{
			TreeTools::scaleTree(* _root, factor);
		}

	protected:
		
		void destroyNode(const N & node)
		{
			if(!node.isLeaf()) {
				for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
					destroyNode(* node[i]);
					delete node[i];
				}
			}
		}
		
};

#endif	//_TREE_H_
