//
// File: Tree.h
// Created by: Julien Dutheil <julien.dutheil@ens-lyon.fr>
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
 

#include "Tree.h"
#include "TreeTools.h"

//From Utils:
#include <Utils/Exceptions.h>

//from the STL:
#include <algorithm>
#include <iostream>
using namespace std;

/******************************************************************************
 *                               The node class                               *
 ******************************************************************************/

/** Default constructor: ******************************************************/

Node::Node() : _id(0), _name(NULL), _father(NULL), _distanceToFather(NULL) {}

Node::Node(int id) : _id(id), _name(NULL), _father(NULL), _distanceToFather(NULL) {}

Node::Node(const string & name) : _id(0), _name(new string(name)), _father(NULL), _distanceToFather(NULL) {}

Node::Node(int id, const string & name) : _id(id), _name(new string(name)), _father(NULL), _distanceToFather(NULL) {}
	
/** Copy constructor: *********************************************************/
	
Node::Node(const Node & node) {
	_id               = node._id;
	_name             = node.hasName() ? new string(* node._name) : NULL;
	_father           = node._father;
	_distanceToFather = node.hasDistanceToFather() ? new double(* node._distanceToFather) : NULL;
	_sons             = node._sons;
	for(map<string, Clonable *>::iterator i = node._properties.begin(); i != node._properties.end(); i++)
		_properties[i -> first] = i -> second -> clone();
}

/** Assignation operator: *****************************************************/

Node & Node::operator=(const Node & node) {
	_id               = node._id;
	_name             = node.hasName() ? new string(* node._name) : NULL;
	_father           = node._father;
	_distanceToFather = node.hasDistanceToFather() ? new double(* node._distanceToFather) : NULL;
	_sons             = node._sons;
	for(map<string, Clonable *>::iterator i = node._properties.begin(); i != node._properties.end(); i++)
		_properties[i -> first] = i -> second -> clone();
	return * this;
}

/** Destructor: ***************************************************************/
	
Node::~Node() {
	delete _name;
	delete _distanceToFather;	
}
	
/** Identity: *****************************************************************/

int Node::getId() const { return _id; }
				
void Node::setId(int id) { _id = id; }
				
/** Name: *********************************************************************/
			
string Node::getName() const throw (NodeException) {
	if(!hasName()) throw NodeException("Node::getName: no name associated to this node.", this);
	return * _name;
}
				
void Node::setName(string & name) { delete _name; _name = new string(name); }

void Node::deleteName() { delete _name; _name = NULL; }

bool Node::hasName() const { return _name != NULL; }

/** Distances: ****************************************************************/
				 
double Node::getDistanceToFather() const throw (NodeException) {
	if(!hasDistanceToFather()) throw NodeException("Node::getDistanceToFather: Node has no distance." , this);
	return * _distanceToFather;
}
				
void Node::setDistanceToFather(double distance) { delete _distanceToFather; _distanceToFather = new double(distance); }

void Node::deleteDistanceToFather() { delete _distanceToFather; _distanceToFather = NULL; }

bool Node::hasDistanceToFather() const { return _distanceToFather != NULL; }
				
/** Father: *******************************************************************/
			
const Node * Node::getFather() const { return _father; }
				
      Node * Node::getFather()       { return _father; }
				
void Node::setFather(Node & node) { _father = & node; }

Node * Node::removeFather() { Node * f = _father; _father = NULL; return f; }
				
bool Node::hasFather() const      { return _father != NULL; }

/** Sons: *********************************************************************/

Node * Node::operator[](int i) { return (i < 0) ? _father : _sons[i]; }
				
const Node * Node::operator[](int i) const { return (i < 0) ? _father : _sons[i]; }
				
unsigned int Node::getNumberOfSons() const { return _sons.size(); }

const Node * Node::getSon(unsigned int i) const { return _sons[i]; }

Node * Node::getSon(unsigned int i) { return _sons[i]; }

void Node::addSon(Node & node)
{
	_sons.push_back(& node);
	node._father = this;
}

void Node::addSon(unsigned int pos, Node & node)
{
	_sons.insert(_sons.begin() + pos, & node);
	node._father = this;
}

void Node::removeSon(unsigned int i) {
	Node * node = _sons[i];
	_sons.erase(_sons.begin() + i);
	node -> removeFather();
}

void Node::removeSon(Node & node) {
	for(unsigned int i = 0; i < _sons.size(); i++) if(* _sons[i] == node) _sons.erase(_sons.begin() + i);
	node.removeFather();
}

void Node::removeSons() {
	while(_sons.size() != 0) removeSon(0);
}

void Node::swap(unsigned int branch1, unsigned int branch2) {
    Node* node1 = getSon(branch1);
    Node* node2 = getSon(branch2);
    removeSon(*node1);
    removeSon(*node2);
    addSon(branch1, *node2);
    addSon(branch2, *node1);
}

/** Properties: ***************************************************************/
				
void Node::setProperty(const string & name, Clonable * property) { _properties[name] = property; }
				
Clonable * Node::getProperty(const string & name) { return _properties[name]; }
				
const Clonable * Node::getProperty(const string & name) const { return const_cast<const Clonable *>(_properties[name]); }
				
void * Node::removeProperty(const string & name)
{
	void * removed = _properties[name];
	_properties.erase(name);
	return removed;
}
				
bool Node::hasProperty(const string & name) const { return _properties.find(name) != _properties.end(); }
				
/** Equality operator: ********************************************************/

bool Node::operator==(const Node & node) const { return _id == node._id; }				
			
/** Tests: ********************************************************************/

bool Node::isLeaf() const { return _sons.size() == 0; }
	
/******************************************************************************/






/*******************************************************************************
 *                                The tree class                               *
 *******************************************************************************/

Tree::Tree() {
	_root = NULL;
}

Tree::Tree(Node & root) {
	_root = &root;
}

Tree::Tree(const Tree & t) {
	//Perform a hard copy of the nodes:
	_root = cloneSubtree(* t.getRootNode());
}

/******************************************************************************/

Tree& Tree::operator=(const Tree & t) {
	//Perform a hard copy of the nodes:
	_root = cloneSubtree(* t.getRootNode());
    return *this;
}

/******************************************************************************/

void Tree::destroyNode(const Node & node) {
	if(!node.isLeaf()) {
		for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
			destroyNode(* node[i]);
			delete node[i];
		}
	}
}

/******************************************************************************/

Tree::~Tree() {
	destroyNode(* _root);
	delete _root;
}

/******************************************************************************/

Node * Tree::cloneSubtree(const Node & node) {
	//First we copy this node using default copy constuctor:
	Node * clone = new Node(node);
	//Remove all sons. This is possible since Tree is a friend of Node:
	//clone -> sons.resize(0);
	//Now we perform a hard copy:
	for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
		//clone -> addSon(* cloneSubtree(* node[i]));
		clone -> _sons[i] = cloneSubtree(* node[i]);
	}
	return clone;
}

/******************************************************************************/

Node * Tree::getRootNode() { return _root; }

const Node * Tree::getRootNode() const { return _root; }

void Tree::setRootNode(Node & root) { _root = & root; }

/******************************************************************************/

string Tree::getName() const { return _name; }
	
void Tree::setName(const string & name) { _name = name; }

/******************************************************************************/

vector<const Node *> Tree::getLeaves() const { return TreeTools::getLeaves(* const_cast<const Node *>(_root)); }

vector<      Node *> Tree::getLeaves()       { return TreeTools::getLeaves(* _root); }

vector<const Node *> Tree::getNodes () const { return TreeTools::getNodes (* const_cast<const Node *>(_root)); }

vector<      Node *> Tree::getNodes ()       { return TreeTools::getNodes (* _root); }

/******************************************************************************/

unsigned int Tree::getNumberOfLeaves() const { return TreeTools::getNumberOfLeaves(* const_cast<const Node *>( _root)); }

vector<string> Tree::getLeavesNames() const { return TreeTools::getLeavesNames(* const_cast<const Node *>( _root)); }

vector<double> Tree::getBranchLengths() const { return TreeTools::getBranchLengths(* _root); }

/******************************************************************************
 *                                 Edition tools                              *
 ******************************************************************************/

vector<Node *> Tree::getPathBetweenAnyTwoNodes(Node & node1, Node & node2) const {

	vector<Node *> path;
	vector<Node *> pathMatrix1;
	vector<Node *> pathMatrix2;

	Node * nodeUp = & node1;
	while(nodeUp != _root)	{
		pathMatrix1.push_back(nodeUp);
		nodeUp = nodeUp -> getFather();
	}
	pathMatrix1.push_back(_root);

	nodeUp = & node2;
	while(nodeUp != _root)	{
		pathMatrix2.push_back(nodeUp);
		nodeUp = nodeUp -> getFather();
	}
	pathMatrix2.push_back(_root);

	int tmp1 = pathMatrix1.size() - 1;
	int tmp2 = pathMatrix2.size() - 1;

	while((tmp1 >= 0) && (tmp2 >= 0)) {
		if (pathMatrix1[tmp1] != pathMatrix2[tmp2]) break;
		tmp1--; tmp2--;
	}

	for (int y = 0; y <= tmp1; ++y) path.push_back(pathMatrix1[y]);
	path.push_back(pathMatrix1[tmp1 + 1]); // pushing once, the Node that was common to both.
	for (int j = tmp2; j >= 0; --j) {
		path.push_back(pathMatrix2[j]);
	}
	return path;
}

/******************************************************************************/

void Tree::rootAt(Node & p_iNewRoot) {
	if (* _root == p_iNewRoot) return;
	vector<Node *> pathMatrix = getPathBetweenAnyTwoNodes(* _root, p_iNewRoot);
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

/******************************************************************************/

bool Tree::isRooted() const { return _root -> getNumberOfSons() == 2; }
		
/******************************************************************************/

bool Tree::unroot() throw (UnrootedTreeException)
{
	if(!isRooted()) throw UnrootedTreeException("Tree::unroot", this);

    if(_root -> getNumberOfSons() == 2) {
		Node* son1 = _root -> getSon(0);
		Node* son2 = _root -> getSon(1);
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

/******************************************************************************/

void Tree::resetNodesId() {
	vector<Node *> nodes = getNodes();
	for(unsigned int i = 0; i < nodes.size(); i++) nodes[i] -> setId(i);
}

/******************************************************************************/

bool Tree::isMultifurcating() const {
	bool b = false;
	for(unsigned int i = 0; i < _root -> getNumberOfSons(); i++) {
		b = b || TreeTools::isMultifurcating(* _root -> getSon(i));
	}
	return b;
}

/******************************************************************************/

/*
void Tree::setNewOutgroup(Tree::Node & node) {
	//cout << "New outgroup " << node.getId() << endl;
	if(node == * _root) return; //can't reroot with root!
	Node * father = node.getFather();
	if(father == _root) return; //nothing to do, already rooted here!
	Node * grandFather = father -> getFather();
	if(grandFather != _root) setNewOutgroup(* father);//must reroot with father first!
	//cout << "dealing with node " << node.getId() << endl;
	
	father -> eraseSon(node);
	_root  -> eraseSon(* father);
	Node * uncle = (* _root)[0];
	_root  -> eraseSon(0);//remove uncle.
	uncle  -> setDistanceToFather(uncle -> getDistanceToFather() + father -> getDistanceToFather());
	father -> setDistanceToFather(node.getDistanceToFather() / 2);//Place new root at mid length between node and father.
	node   .  setDistanceToFather(node.getDistanceToFather() / 2);
	father -> addSon(* uncle);
	_root  -> addSon(node);
	_root  -> addSon(* father);//We add it again ;-)
}*/
