//
// File: TreeTemplate.h
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

#ifndef _TREETEMPLATE_H_
#define _TREETEMPLATE_H_

#include "TreeExceptions.h"
#include "TreeTools.h"
#include "Tree.h"

// From the STL:
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

using namespace std;

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
 * To clone a tree from from another tree with a different template,
 * consider using the TreeTools::cloneSutree<N>() method:
 * <code>
 * Tree * t = new Tree<Node>(...)
 * NodeTemplate<int> * newRoot = TreeTools::cloneSubtree< NodeTemplate<int> >(* (t -> getRootNode()))
 * Tree< NodeTemplate<int> > * tt = new Tree< NodeTemplate<int> >(* newRoot);
 * </code>
 * 
 * @see Node
 * @see NodeTemplate
 * @see TreeTools
 */
template<class N>
class TreeTemplate: public virtual Tree
{

	/**
	 * Fields:
	 */
	protected:
		N * _root;
		string _name;
    int _nextNodeId;

	public: // Constructors and destructor:
		
		TreeTemplate() { _root = NULL; }

		TreeTemplate(const TreeTemplate<N> & t)
		{
			//Perform a hard copy of the nodes:
			_root = TreeTools::cloneSubtree<N>(* t.getRootNode());
		}

		TreeTemplate(N & root) { _root = &root; }

		TreeTemplate<N> & operator=(const TreeTemplate<N> & t)
		{
			//Perform a hard copy of the nodes:
			_root = TreeTools::cloneSubtree<N>(* t.getRootNode());
    	return *this;
		}

		virtual ~TreeTemplate() { destroyNode(* _root); delete _root; }

			
	/**
	 * Methods:
 	 */
	
	public:
		
		string getName() const { return _name; }
	
		void setName(const string & name) { _name = name; }

		int getRootId() const { return _root -> getId(); }

		unsigned int getNumberOfLeaves() const { return TreeTools::getNumberOfLeaves(* _root); }
		
		unsigned int getNumberOfNodes() const { return TreeTools::getNumberOfNodes(* _root); }
		
		int getLeafId(const string & name) const throw (NodeNotFoundException) { return TreeTools::getLeafId(* _root, name); }
		
		vector<int> getLeavesId() const { return TreeTools::getLeavesId(* _root); }

		vector<int> getNodesId() const { return TreeTools::getNodesId(* _root); }

		vector<double> getBranchLengths() const { return TreeTools::getBranchLengths(* _root); }

		vector<string> getLeavesNames() const { return TreeTools::getLeavesNames(* const_cast<const N *>( _root)); }

		vector<int> getSonsId(int parentId) const throw (NodeNotFoundException)	{ return getNode(parentId) -> getSonsId(); }

		int getFatherId(int parentId) const throw (NodeNotFoundException) { return getNode(parentId) -> getFatherId(); }

		bool hasFather(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId) -> hasFather(); }
	
		string getNodeName(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId) -> getName(); }

		bool hasNodeName(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId) -> hasName(); }
		
		void setNodeName(int nodeId, const string & name) throw (NodeNotFoundException) { getNode(nodeId) -> setName(name); }
		
		void deleteNodeName(int nodeId) throw (NodeNotFoundException) { return getNode(nodeId) -> deleteName(); }

		bool isLeaf(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId) -> isLeaf(); }

		bool isRoot(int nodeId) const throw (NodeNotFoundException) { return TreeTools::isRoot(* getNode(nodeId)); }

		double getDistanceToFather(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId) -> getDistanceToFather(); }
		
		void setDistanceToFather(int nodeId, double length) throw (NodeNotFoundException) { getNode(nodeId) -> setDistanceToFather(length); }
		
		void deleteDistanceToFather(int nodeId) throw (NodeNotFoundException) { getNode(nodeId) -> deleteDistanceToFather(); }
		
		bool hasDistanceToFather(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId) -> hasDistanceToFather(); }

		bool hasProperty(int nodeId, const string & name) const throw (NodeNotFoundException) { return getNode(nodeId) -> hasProperty(name); }
	
		void setProperty(int nodeId, const string & name, Clonable * property) throw (NodeNotFoundException) { getNode(nodeId) -> setProperty(name, property); }
				
		Clonable * getProperty(int nodeId, const string & name) throw (NodeNotFoundException) { return getNode(nodeId) -> getProperty(name); }
				
		const Clonable * getProperty(int nodeId, const string & name) const throw (NodeNotFoundException) { return getNode(nodeId) -> getProperty(name); }
				
		Clonable * removeProperty(int nodeId, const string & name) throw (NodeNotFoundException) { return getNode(nodeId) -> removeProperty(name); }
		
    void rootAt(int nodeId) throw (NodeNotFoundException)	{	rootAt(* getNode(nodeId)); }
		
		void newOutGroup(int nodeId) throw (NodeNotFoundException) {	newOutGroup(* getNode(nodeId)); }

		bool isRooted() const { return _root -> getNumberOfSons() == 2; }
		
		bool unroot() throw (UnrootedTreeException)
		{
			if(!isRooted()) throw UnrootedTreeException("Tree::unroot", this);

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
        delete _root;
				return true;
			} else return false; // Tree is already unrooted.
		}

		void resetNodesId()
		{
			vector<N *> nodes = getNodes();
      _nextNodeId = 0;
			for(unsigned int i = 0; i < nodes.size(); i++) {
        nodes[i] -> setId(i);
        _nextNodeId++;
      }
		}
		
		bool isMultifurcating() const
		{
			bool b = false;
			for(unsigned int i = 0; i < _root -> getNumberOfSons(); i++) {
				b = b || TreeTools::isMultifurcating(* _root -> getSon(i));
			}
			return b;
		}
		
		vector<double> getBranchLengths() throw (NodeException)
		{
			Vdouble brLen(1);
			for(unsigned int i = 0; i < _root -> getNumberOfSons(); i++) {
				Vdouble sonBrLen = TreeTools::getBranchLengths(* _root -> getSon(i));
				for(unsigned int j = 0; j < sonBrLen.size(); j++) brLen.push_back(sonBrLen[j]);
			}
			return brLen;
		}

		double getTotalLength() throw (NodeException)
		{
			double length = 0;
			for(unsigned int i = 0; i < _root->getNumberOfSons(); i++) {
        length += TreeTools::getTotalLength(*_root->getSon(i));
      }
      return length;
		}

		void setBranchLengths(double brLen)
		{
			for(unsigned int i = 0; i < _root->getNumberOfSons(); i++) {
        TreeTools::setBranchLengths(*_root->getSon(i), brLen);
      }
		}
		
		void setVoidBranchLengths(double brLen)
		{
			for(unsigned int i = 0; i < _root->getNumberOfSons(); i++) {
			  TreeTools::setVoidBranchLengths(*_root->getSon(i), brLen);
      }
		}
	
		void scaleTree(double factor) throw (NodeException)
		{
			for(unsigned int i = 0; i < _root->getNumberOfSons(); i++) {
			  TreeTools::scaleTree(*_root->getSon(i), factor);
      }
		}

    int getNextId() { return _nextNodeId++; }

		/**
		 * @name Specific methods
		 *
		 * @{
		 */
		virtual void setRootNode(N & root) { _root = & root; }
	
		virtual N * getRootNode() { return _root; }
		
		virtual const N * getRootNode() const { return _root; }
		
		virtual vector<const N *> getLeaves() const { return TreeTools::getLeaves(* const_cast<const N *>(_root)); }

		virtual vector<N *> getLeaves() { return TreeTools::getLeaves(* _root); }
	
		virtual vector<const N *> getNodes() const { return TreeTools::getNodes(* const_cast<const N *>(_root)); }

		virtual vector<N *> getNodes() { return TreeTools::getNodes(* _root); }
		
		virtual vector<const N *> getInnerNodes() const { return TreeTools::getInnerNodes(* const_cast<const N *>(_root)); }

		virtual vector<N *> getInnerNodes() { return TreeTools::getInnerNodes(* _root); }
	
		virtual N * getNode(int id) throw (NodeNotFoundException, Exception)
		{
			vector<N *> nodes = TreeTools::searchNodeWithId<N>(*_root, id);
			if(nodes.size() > 1) throw Exception("TemplateTree::rootAt(): Non-unique id!");
			if(nodes.size() == 0) throw NodeNotFoundException("TemplateTree::rootAt(): Node with id not found.", "" + id);
			return nodes[0];
		}
		
		virtual const N * getNode(int id) const throw (NodeNotFoundException, Exception)
		{
			vector<const N *> nodes = TreeTools::searchNodeWithId<const N>(*_root, id);
			if(nodes.size() > 1) throw Exception("TemplateTree::rootAt(): Non-unique id!");
			if(nodes.size() == 0) throw NodeNotFoundException("TemplateTree::rootAt(): Node with id not found.", "" + id);
			return nodes[0];
		}

		void rootAt(N & newRoot)
	  {
			if (* _root == newRoot) return;
			if(isRooted()) unroot();
			vector<Node *> path = TreeTools::getPathBetweenAnyTwoNodes(* _root, newRoot);

			for (unsigned int i = 0; i < path.size() - 1 ; i++) {
				//pathMatrix[i] -> _father = pathMatrix[i + 1];
				//pathMatrix[i] -> setDistanceToFather(pathMatrix[i + 1] -> getDistanceToFather());
				//typename vector<Node *>::iterator vec_iter;
				//vec_iter = remove(pathMatrix[i] -> _sons.begin(), pathMatrix[i] -> _sons.end(), pathMatrix[i + 1]);
				//pathMatrix[i] -> _sons.erase(vec_iter, pathMatrix[i] -> _sons.end()); // pg 1170, primer.	
				//pathMatrix[i+1] -> _sons.push_back(pathMatrix[i + 1] -> getFather());
				//pathMatrix[i+1] -> _father = NULL;
				path[i] -> removeSon(* path[i+1]);
        if(path[i+1]->hasDistanceToFather()) path[i]->setDistanceToFather(path[i+1]->getDistanceToFather());
        else path[i]->deleteDistanceToFather();
				path[i+1] -> addSon(* path[i]);
			}
      newRoot.deleteDistanceToFather();
			_root = & newRoot;
		}

		void newOutGroup(N & outGroup)
		{
			if(* _root == outGroup) return;
			int rootId;
			if(isRooted()) {
				rootId = getRootId();
				unroot();
			} else {
				rootId = getNextId();
			}
			rootAt(* outGroup.getFather());
			N* oldRoot = _root;
			oldRoot -> removeSon(outGroup);
			_root = new N();
			_root -> setId(rootId);
			_root -> addSon(*oldRoot);
			_root -> addSon(outGroup);
      // Check lengths:
      if(outGroup.hasDistanceToFather()) {
        double l = outGroup.getDistanceToFather() / 2.;
        outGroup.setDistanceToFather(l);
        oldRoot->setDistanceToFather(l);
      }
		}

		/** @} */
	
		
	protected:
		
		virtual void destroyNode(const N & node)
		{
			for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
				const N * son = node.getSon(i);
				destroyNode(* son);
				delete son;
			}
		}
		
};

#endif	//_TREETEMPLATE_H_

