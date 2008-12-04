//
// File: TreeTemplate.h
// Created by: Julien Dutheil
//             Celine Scornavacca
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
#include "TreeTemplateTools.h"
#include "Tree.h"

// From the STL:
#include <string>
#include <vector>
#include <map>

using namespace std;

namespace bpp
{

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
 * @code
 * Tree * t = new Tree<Node>(...)
 * NodeTemplate<int> * newRoot = TreeTools::cloneSubtree< NodeTemplate<int> >(* (t -> getRootNode()))
 * Tree< NodeTemplate<int> > * tt = new Tree< NodeTemplate<int> >(* newRoot);
 * @endcode
 *
 * The getNextId() method sends a id value which is not used in the tree.
 * In the current implementation, it uses the TreeTools::getMPNUId() method.
 * This avoids to use duplicated ids, but is time consuming.
 * In most cases, it is of better efficiency if the user deal with the ids himself, by using the Node::setId() method.
 * The TreeTools::getMaxId() method may also prove useful in this respect.
 * The resetNodesId() method can also be used to re-initialize all ids.
 * 
 * @see Node
 * @see NodeTemplate
 * @see TreeTools
 */
template<class N>
class TreeTemplate:
  public Tree
{
	/**
	 * Fields:
	 */
	protected:
		N * _root;
		string _name;

	public: // Constructors and destructor:
		
		TreeTemplate(): _root(NULL), _name() {}

		TreeTemplate(const TreeTemplate<N> & t):
      _root(NULL), _name(t._name)
		{
			//Perform a hard copy of the nodes:
			_root = TreeTemplateTools::cloneSubtree<N>(* t.getRootNode());
		}

		TreeTemplate(const Tree & t):
      _root(NULL), _name(t.getName())
		{
			//Create new nodes from an existing tree:
			_root = TreeTemplateTools::cloneSubtree<N>(t, t.getRootId());
		}

		TreeTemplate(N & root): _root(&root), _name() {}

		TreeTemplate<N> & operator=(const TreeTemplate<N> & t)
		{
			//Perform a hard copy of the nodes:
			if(_root) { destroySubtree(*_root); delete _root; }
      _root = TreeTemplateTools::cloneSubtree<N>(* t.getRootNode());
      _name = t._name;
    	return *this;
		}

		virtual ~TreeTemplate()
    {
      destroySubtree(* _root);
      delete _root;
    }

    TreeTemplate<N> * clone() const { return new TreeTemplate<N>(*this); }
			
	/**
	 * Methods:
 	 */
	
	public:
		
		string getName() const { return _name; }
	
		void setName(const string & name) { _name = name; }

		int getRootId() const { return _root->getId(); }

		unsigned int getNumberOfLeaves() const { return TreeTemplateTools::getNumberOfLeaves(* _root); }
		
		unsigned int getNumberOfNodes() const { return TreeTemplateTools::getNumberOfNodes(* _root); }
		
		int getLeafId(const string & name) const throw (NodeNotFoundException) { return TreeTemplateTools::getLeafId(* _root, name); }
		
		vector<int> getLeavesId() const { return TreeTemplateTools::getLeavesId(* _root); }

		vector<int> getNodesId() const { return TreeTemplateTools::getNodesId(* _root); }

		vector<double> getBranchLengths() const { return TreeTemplateTools::getBranchLengths(* _root); }

		vector<string> getLeavesNames() const { return TreeTemplateTools::getLeavesNames(* const_cast<const N *>( _root)); }

		vector<int> getSonsId(int parentId) const throw (NodeNotFoundException)	{ return getNode(parentId)->getSonsId(); }

		int getFatherId(int parentId) const throw (NodeNotFoundException) { return getNode(parentId)->getFatherId(); }

		bool hasFather(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId)->hasFather(); }
	
		string getNodeName(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId)->getName(); }

		bool hasNodeName(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId)->hasName(); }
		
		void setNodeName(int nodeId, const string & name) throw (NodeNotFoundException) { getNode(nodeId)->setName(name); }
		
		void deleteNodeName(int nodeId) throw (NodeNotFoundException) { return getNode(nodeId)->deleteName(); }

    bool hasNode(int nodeId) const { return TreeTemplateTools::hasNodeWithId(*_root, nodeId); }

		bool isLeaf(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId)->isLeaf(); }

		bool isRoot(int nodeId) const throw (NodeNotFoundException) { return TreeTemplateTools::isRoot(* getNode(nodeId)); }

		double getDistanceToFather(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId)->getDistanceToFather(); }
		
		void setDistanceToFather(int nodeId, double length) throw (NodeNotFoundException) { getNode(nodeId)->setDistanceToFather(length); }
		
		void deleteDistanceToFather(int nodeId) throw (NodeNotFoundException) { getNode(nodeId)->deleteDistanceToFather(); }
		
		bool hasDistanceToFather(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId)->hasDistanceToFather(); }

		bool hasNodeProperty(int nodeId, const string & name) const throw (NodeNotFoundException) { return getNode(nodeId)->hasNodeProperty(name); }
	
		void setNodeProperty(int nodeId, const string & name, const Clonable & property) throw (NodeNotFoundException) { getNode(nodeId)->setNodeProperty(name, property); }
				
		Clonable * getNodeProperty(int nodeId, const string & name) throw (NodeNotFoundException) { return getNode(nodeId)->getNodeProperty(name); }
				
		const Clonable * getNodeProperty(int nodeId, const string & name) const throw (NodeNotFoundException) { return getNode(nodeId)->getNodeProperty(name); }
				
		Clonable * removeNodeProperty(int nodeId, const string & name) throw (NodeNotFoundException) { return getNode(nodeId)->removeNodeProperty(name); }
		
    vector<string> getNodePropertyNames(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId)->getNodePropertyNames(); }
		
		bool hasBranchProperty(int nodeId, const string & name) const throw (NodeNotFoundException) { return getNode(nodeId)->hasBranchProperty(name); }
	
		void setBranchProperty(int nodeId, const string & name, const Clonable & property) throw (NodeNotFoundException) { getNode(nodeId)->setBranchProperty(name, property); }
				
		Clonable * getBranchProperty(int nodeId, const string & name) throw (NodeNotFoundException) { return getNode(nodeId)->getBranchProperty(name); }
				
		const Clonable * getBranchProperty(int nodeId, const string & name) const throw (NodeNotFoundException) { return getNode(nodeId)->getBranchProperty(name); }
				
		Clonable * removeBranchProperty(int nodeId, const string & name) throw (NodeNotFoundException) { return getNode(nodeId)->removeBranchProperty(name); }
		
    vector<string> getBranchPropertyNames(int nodeId) const throw (NodeNotFoundException) { return getNode(nodeId)->getBranchPropertyNames(); }
		
    void rootAt(int nodeId) throw (NodeNotFoundException)	{	rootAt(* getNode(nodeId)); }
		
		void newOutGroup(int nodeId) throw (NodeNotFoundException) {	newOutGroup(* getNode(nodeId)); }

		bool isRooted() const { return _root->getNumberOfSons() == 2; }
		
		bool unroot() throw (UnrootedTreeException)
		{
			if(!isRooted()) throw UnrootedTreeException("Tree::unroot", this);

    	if(_root->getNumberOfSons() == 2)
      {
				N* son1 = _root->getSon(0);
				N* son2 = _root->getSon(1);
				if(son1->isLeaf() && son2->isLeaf()) return false; // We can't unroot a single branch!
				
        // We manage to have a subtree in position 0:
				if(son1->isLeaf())
        {
				  _root->swap(0, 1);
				  son1 = _root->getSon(0);
				  son2 = _root->getSon(1);
				}

				// Take care of branch lengths:
				if(son1->hasDistanceToFather())
        {
					if(son2->hasDistanceToFather())
          {
					  // Both nodes have lengths, we sum them:
					  son2->setDistanceToFather(son1->getDistanceToFather() + son2->getDistanceToFather());
				  }
          else
          {
					  // Only node 1 has length, we set it to node 2:
					  son2->setDistanceToFather(son1->getDistanceToFather());
				  }
				  son1->deleteDistanceToFather();
				} // Else node 2 may or may not have a branch length, we do not care!

				// Remove the root:
				_root->removeSons();
				son1->addSon(*son2);
        delete _root;
				setRootNode(*son1);
				return true;
			} else return false; // Tree is already unrooted.
		}

		void resetNodesId()
		{
			vector<N *> nodes = getNodes();
			for(unsigned int i = 0; i < nodes.size(); i++)
      {
        nodes[i]->setId(i);
      }
		}
		
		bool isMultifurcating() const
		{
			bool b = false;
			for(unsigned int i = 0; i < _root -> getNumberOfSons(); i++)
      {
				b = b || TreeTemplateTools::isMultifurcating(* _root->getSon(i));
			}
			return b;
		}
		
		vector<double> getBranchLengths() throw (NodeException)
		{
			Vdouble brLen(1);
			for(unsigned int i = 0; i < _root->getNumberOfSons(); i++)
      {
				Vdouble sonBrLen = TreeTemplateTools::getBranchLengths(* _root->getSon(i));
				for(unsigned int j = 0; j < sonBrLen.size(); j++) brLen.push_back(sonBrLen[j]);
			}
			return brLen;
		}

		double getTotalLength() throw (NodeException)
		{
      return TreeTemplateTools::getTotalLength(*_root, false);
		}

		void setBranchLengths(double brLen)
		{
			for(unsigned int i = 0; i < _root->getNumberOfSons(); i++)
      {
        TreeTemplateTools::setBranchLengths(*_root->getSon(i), brLen);
      }
		}
		
		void setVoidBranchLengths(double brLen)
		{
			for(unsigned int i = 0; i < _root->getNumberOfSons(); i++)
      {
			  TreeTemplateTools::setVoidBranchLengths(*_root->getSon(i), brLen);
      }
		}
	
		void scaleTree(double factor) throw (NodeException)
		{
			for(unsigned int i = 0; i < _root->getNumberOfSons(); i++)
      {
			  TreeTemplateTools::scaleTree(*_root->getSon(i), factor);
      }
		}

    int getNextId()
    { 
      return TreeTools::getMPNUId(*this, _root->getId());
    }

    void swapNodes(int parentId, unsigned int i1, unsigned int i2) throw (NodeNotFoundException,IndexOutOfBoundsException)
    {
			vector<N *> nodes = TreeTemplateTools::searchNodeWithId<N>(*_root, parentId);
			if(nodes.size() == 0) throw NodeNotFoundException("TreeTemplate:swapNodes(): Node with id not found.", "" + parentId);
      for(unsigned int i = 0; i < nodes.size(); i++) nodes[i]->swap(i1, i2); 
    }
 

		/**
		 * @name Specific methods
		 *
		 * @{
		 */
		virtual void setRootNode(N & root) { _root = & root; }
	
		virtual N * getRootNode() { return _root; }
		
		virtual const N * getRootNode() const { return _root; }
		
		virtual vector<const N *> getLeaves() const { return TreeTemplateTools::getLeaves(* const_cast<const N *>(_root)); }

		virtual vector<N *> getLeaves() { return TreeTemplateTools::getLeaves(* _root); }
	
		virtual vector<const N *> getNodes() const { return TreeTemplateTools::getNodes(* const_cast<const N *>(_root)); }

		virtual vector<N *> getNodes() { return TreeTemplateTools::getNodes(* _root); }
		
		virtual vector<const N *> getInnerNodes() const { return TreeTemplateTools::getInnerNodes(* const_cast<const N *>(_root)); }

		virtual vector<N *> getInnerNodes() { return TreeTemplateTools::getInnerNodes(* _root); }
	
		virtual N * getNode(int id) throw (NodeNotFoundException, Exception)
		{
			vector<N *> nodes;
      TreeTemplateTools::searchNodeWithId<N>(*_root, id, nodes);
			if(nodes.size() > 1) throw Exception("TreeTemplate::getNode(): Non-unique id! (" + TextTools::toString(id) + ").");
			if(nodes.size() == 0) throw NodeNotFoundException("TreeTemplate::getNode(): Node with id not found.", "" + id);
			return nodes[0];
		}
		
		virtual const N * getNode(int id) const throw (NodeNotFoundException, Exception)
		{
			vector<const N *> nodes;
      TreeTemplateTools::searchNodeWithId<const N>(*_root, id, nodes);
			if(nodes.size() > 1) throw Exception("TreeTemplate::getNode(): Non-unique id! (" + TextTools::toString(id) + ").");
			if(nodes.size() == 0) throw NodeNotFoundException("TreeTemplate::getNode(): Node with id not found.", "" + id);
			return nodes[0];
		}

	  virtual N * getNode(const string & name) throw (NodeNotFoundException, Exception)
    {	
		  vector<N *> nodes;
      TreeTemplateTools::searchNodeWithName(*_root, name, nodes);
  		if(nodes.size() > 1)  throw NodeNotFoundException("TreeTemplate::getNode(): Non-unique name.", "" + name);
	  	if(nodes.size() == 0) throw NodeNotFoundException("TreeTemplate::getNode(): Node with name not found.", "" + name);
  		return nodes[0];
	  }

    virtual const N * getNode(const string & name) const throw (NodeNotFoundException, Exception)
    {	
		  vector<const N *> nodes;
      TreeTemplateTools::searchNodeWithName<const N>(*_root, name, nodes);
  		if(nodes.size() > 1)  throw NodeNotFoundException("TreeTemplate::getNode(): Non-unique name.", "" + name);
	  	if(nodes.size() == 0) throw NodeNotFoundException("TreeTemplate::getNode(): Node with name not found.", "" + name);
  		return nodes[0];
	  }

		void rootAt(N & newRoot)
	  {
			if (* _root == newRoot) return;
			if(isRooted()) unroot();
			vector<Node *> path = TreeTemplateTools::getPathBetweenAnyTwoNodes(* _root, newRoot);

			for (unsigned int i = 0; i < path.size() - 1 ; i++)
      {
				//pathMatrix[i] -> _father = pathMatrix[i + 1];
				//pathMatrix[i] -> setDistanceToFather(pathMatrix[i + 1] -> getDistanceToFather());
				//typename vector<Node *>::iterator vec_iter;
				//vec_iter = remove(pathMatrix[i] -> _sons.begin(), pathMatrix[i] -> _sons.end(), pathMatrix[i + 1]);
				//pathMatrix[i] -> _sons.erase(vec_iter, pathMatrix[i] -> _sons.end()); // pg 1170, primer.	
				//pathMatrix[i+1] -> _sons.push_back(pathMatrix[i + 1] -> getFather());
				//pathMatrix[i+1] -> _father = NULL;
				path[i]->removeSon(* path[i+1]);
        if(path[i+1]->hasDistanceToFather()) path[i]->setDistanceToFather(path[i+1]->getDistanceToFather());
        else path[i]->deleteDistanceToFather();
				path[i+1]->addSon(* path[i]);

        vector<string> names = path[i+1]->getBranchPropertyNames();
        for(unsigned int j = 0; j < names.size(); j++)
        {
          path[i]->setBranchProperty(names[j], *path[i+1]->getBranchProperty(names[j]));
        }
        path[i+1]->deleteBranchProperties();
      }
      newRoot.deleteDistanceToFather();
      newRoot.deleteBranchProperties();
			_root = & newRoot;
		}

		void newOutGroup(N & outGroup)
		{
			if(* _root == outGroup) return;
			int rootId;
			if(isRooted())
      {
        for(unsigned int i = 0; i < _root->getNumberOfSons(); i++)
        {
          if(*_root->getSon(i) == outGroup) return; //This tree is already rooted appropriately.
        }
				rootId = getRootId();
				unroot();
			}
      else
      {
				rootId = getNextId();
			}
			rootAt(* outGroup.getFather());
			N* oldRoot = _root;
			oldRoot->removeSon(outGroup);
			_root = new N();
			_root->setId(rootId);
			_root->addSon(*oldRoot);
			_root->addSon(outGroup);
      // Check lengths:
      if(outGroup.hasDistanceToFather())
      {
        double l = outGroup.getDistanceToFather() / 2.;
        outGroup.setDistanceToFather(l);
        oldRoot->setDistanceToFather(l);
      }
		}

		/** @} */
		
	protected:
		
		virtual void destroySubtree(const N & node)
		{
			for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
      {
				const N * son = node.getSon(i);
				destroySubtree(* son);
				delete son;
			}
		}
		
};

} //end of namespace bpp.

#endif	//_TREETEMPLATE_H_

