//
// File: TreeTools.h
// Created by:  <@bogdanof>
// Created on: Wed Aug  6 13:45:28 2003
//

#include "TreeTools.h"

// From Utils:
#include <Utils/TextTools.h>
#include <Utils/Number.h>

// From NumCalc:
#include <NumCalc/RandomTools.h>

// From the STL:
#include <iostream>
#include <sstream>

using namespace std;

/******************************************************************************/

TreeTools::~TreeTools() {}

/******************************************************************************/

vector<const Node *> TreeTools::getLeaves(const Node & node)
{
	vector<const Node *> leaves;
	if(node.isLeaf()) {
		leaves.push_back(& node);
	} else {
		for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
			vector<const Node *> sonLeaves = getLeaves(* node[i]);
			for(unsigned int j = 0; j < sonLeaves.size(); j++) {
				leaves.push_back(sonLeaves[j]);
			}
		}
	}
	return leaves;
}

/******************************************************************************/

vector<Node *> TreeTools::getLeaves(Node & node)
{
	vector<Node *> leaves;
	if(node.isLeaf()) {
		leaves.push_back(& node);
	} else {
		for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
			vector<Node *> sonLeaves = getLeaves(* node[i]);
			for(unsigned int j = 0; j < sonLeaves.size(); j++) {
				leaves.push_back(sonLeaves[j]);
			}
		}
	}
	return leaves;
}

/******************************************************************************/

vector<const Node *> TreeTools::getNodes(const Node & node)
{
	vector<const Node *> nodes;
	for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
		vector<const Node *> sonNodes = getNodes(* node.getSon(i));
		for(unsigned int j = 0; j < sonNodes.size(); j++) {
			nodes.push_back(sonNodes[j]);
		}
	}
	nodes.push_back(& node);
	return nodes;
}

/******************************************************************************/

vector<Node *> TreeTools::getNodes(Node & node)
{
	vector<Node *> nodes;
	for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
		vector<Node *> sonNodes = getNodes(* node[i]);
		for(unsigned int j = 0; j < sonNodes.size(); j++) {
			nodes.push_back(sonNodes[j]);
		}
	}
	nodes.push_back(& node);
	return nodes;
}

/******************************************************************************/

bool TreeTools::isRoot(const Node & node) { return !node.hasFather(); }

/******************************************************************************/

unsigned int TreeTools::getNumberOfLeaves(const Node & node)
{
	unsigned int nbLeaves = 0;
	if(node.isLeaf()) {
		nbLeaves++;
	} else {
		for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
			nbLeaves += getNumberOfLeaves(* node[i]);
		}
	}
	return nbLeaves;
}

/******************************************************************************/

vector<string> TreeTools::getLeavesNames(const Node & node)
{
	vector<string> names;
	if(node.isLeaf()) {
    names.push_back(node.getName());
  } else {
		for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
			vector<string> subNames = getLeavesNames(* node.getSon(i));
			for(unsigned int j = 0; j < subNames.size(); j++) names.push_back(subNames[j]);
		}
	}
	return names;	 
}

/******************************************************************************/

unsigned int TreeTools::getDepth(const Node & node)
{
	unsigned int d = 0;
	for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
		unsigned int c = getDepth(* node[i]) + 1;
		if( c > d) d = c;
	}
	return d;
}

/******************************************************************************/

TreeTools::Element TreeTools::getElement(string elt) throw (IOException)
{
	Element element;
	element.length    = NULL; //default
	element.bootstrap = NULL; //default
	
	int colon = elt.rfind(':');
	try {
		string elt2;
		if(colon >= 0 && colon < (int)elt.size()) {
			//this is an element with length:
			elt2 = elt.substr(0, colon);
			istringstream iss(elt.substr(colon + 1));
			double length;
			iss >> length;
			element.length = new double(length);
		} else {
			//this is an element without length;
			elt2 = elt;
		}
	
		unsigned int  lastP = elt2.rfind(')');
		unsigned int firstP = elt2.find('(');
		if(firstP > elt2.size()) {
			//This is a leaf:
			element.content = elt2;
		} else {
			//This is a node:
			if(lastP < firstP) throw IOException("Invalid format: bad closing parenthesis in " + elt2);
			element.content = elt2.substr(firstP + 1, lastP - firstP - 1);
			string bootstrap = elt2.substr(lastP + 1);
			//cout << "ELEMENT: BOOTSTRAP: " << bootstrap << endl;
			if(!TextTools::isEmpty(bootstrap)) {
				//this a node with a bootstrap value:
				istringstream iss(bootstrap);
				double bootstrapValue;
				iss >> bootstrapValue;
				element.bootstrap = new double(bootstrapValue);
			}
		}
	} catch(exception e) {
		throw IOException("Bad tree description: " + elt);
	}
	return element;
}	

/******************************************************************************/

class NodeTokenizer
{
	protected:
		vector<string> tokens;
		mutable unsigned int currentPosition;
	
	public:
		NodeTokenizer(const string & description) throw (IOException)
		{
			//cout << "NODETOENIZER: " << description << endl;
			unsigned int tokCount = 0;
			int parCount = 0;
			unsigned int i;
			for(i = 0; i < description.size(); i++) {
				if(description[i] == '(') parCount++; //Another open parenthesis
				if(description[i] == ')') parCount--; //Another close parenthesis
				if(parCount < 0) throw IOException("Invalid tree description: closing parenthesis with no opening one, in " + description);
				if(description[i] == ',' && parCount == 0) {
					//New token found:
					//cout << "NODETOENIZER: NEWTOKEN " << description.substr(tokCount, i - tokCount - 1) << endl;
					tokens.push_back(description.substr(tokCount, i - tokCount));
					tokCount = i + 1;
				}					
			}
			//Add last token:
			//cout << "NODETOENIZER: NEWTOKEN " << description.substr(tokCount) << endl;
			tokens.push_back(description.substr(tokCount));
			
			currentPosition = 0;
		}
		
	public:
		string next() const { string s = tokens[currentPosition]; currentPosition++; return s; }
		bool hasNext() const { return currentPosition < tokens.size(); }
};

/******************************************************************************/
Node * TreeTools::parenthesisToNode(const string & description)
{
	//cout << "NODE: " << description << endl;
	Element elt = getElement(description);

	//New node:
	Node * node = new Node();
	if(elt.length != NULL) {
		node -> setDistanceToFather(* elt.length);
		//cout << "NODE: LENGTH: " << * elt.length << endl;
	}
	if(elt.bootstrap != NULL) {
		node -> setProperty(BOOTSTRAP, new Number<double>(* elt.bootstrap));
		//cout << "NODE: BOOTSTRAP: " << * elt.bootstrap << endl;
	}
	
	NodeTokenizer nt(elt.content);
	vector<string> elements;
	while(nt.hasNext()) {
		elements.push_back(nt.next());
	}

	if(elements.size() == 1) {
		//This is a leaf:
		//cout << "NODE: LEAF: " << elements[0] << endl;
		string name = TextTools::removeSurroundingWhiteSpaces(elements[0]);
		node -> setName(name);
	} else {
		//This is a node:
		for(unsigned int i = 0; i < elements.size(); i++) {
			//cout << "NODE: SUBNODE: " << i << ", " << elements[i] << endl;
			Node * son = parenthesisToNode(elements[i]);
			node -> addSon(* son);
		}
	}
	//Must delete the element:
	delete elt.length;
	delete elt.bootstrap;
	return node;
}

/******************************************************************************/

Tree * TreeTools::parenthesisToTree(const string & description)
{
	int lastP  = description.rfind(')');
	int firstP = description.find('(');
	string content = description.substr(firstP + 1, lastP - firstP - 1);
	//cout << "TREE: " << content << endl;
	//New root node:
	Node * node = new Node();
	
	NodeTokenizer nt(content);
	vector<string> elements;
	while(nt.hasNext()) {
		elements.push_back(nt.next());
	}

	if(elements.size() == 1) {
		//This is a leaf:
		node -> setName(elements[0]);
	} else {
		//This is a node:
		for(unsigned int i = 0; i < elements.size(); i++) {
			Node * son = parenthesisToNode(elements[i]);
			node -> addSon(* son);
		}
	}
	Tree * tree = new Tree();
	tree -> setRootNode(* node);
	tree -> resetNodesId();
	return tree;
}

/******************************************************************************/

string TreeTools::nodeToParenthesis(const Node & node) {
	ostringstream s;
	if(node.isLeaf()) {
		s << node.getName();
	} else {
		s << "(";
		s << nodeToParenthesis(* node[0]);
		for(unsigned int i = 1; i < node.getNumberOfSons(); i++) {
			s << "," << nodeToParenthesis(* node[i]);
		}
		s << ")";
	}
	if(node.hasProperty(BOOTSTRAP)) s << * (const double *)(node.getProperty(BOOTSTRAP));
	if(node.hasDistanceToFather()) s << ":" << node.getDistanceToFather();
	return s.str();	
}

/******************************************************************************/

string TreeTools::treeToParenthesis(const Tree & tree) {
	ostringstream s;
	s << "(";
	const Node * node = tree.getRootNode();
	if(node -> isLeaf()) {
		s << node -> getName();
	} else {
		s << nodeToParenthesis(* node -> getSon(0));
		for(unsigned int i = 1; i < node -> getNumberOfSons(); i++) {
			s << "," << nodeToParenthesis(* node -> getSon(i));
		}
	}
	s << ");" << endl;
	return s.str();	
}

/******************************************************************************/

bool TreeTools::isMultifurcating(const Node & node) {
	if(node.getNumberOfSons() > 2) return true;
	else {
		bool b = false;
		for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
			b = b || isMultifurcating(* node.getSon(i));
		}
		return b;
	}		
}

/******************************************************************************/

Vdouble TreeTools::getBranchLengths(const Node & node) throw (NodeException)
{
	Vdouble brLen(1);
	if(node.hasDistanceToFather()) brLen[0] = node.getDistanceToFather();
	else throw NodeException("TreeTools::getbranchLengths(). No branch length.", &node);
	for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
		Vdouble sonBrLen = getBranchLengths(* node.getSon(i));
		for(unsigned int j = 0; j < sonBrLen.size(); j++) brLen.push_back(sonBrLen[j]);
	}
	return brLen;
}

/******************************************************************************/

Vdouble TreeTools::getBranchLengths(const Tree & tree) throw (NodeException)
{
	Vdouble brLen(1);
	const Node * root = tree.getRootNode();
	for(unsigned int i = 0; i < root -> getNumberOfSons(); i++) {
		Vdouble sonBrLen = getBranchLengths(* root -> getSon(i));
		for(unsigned int j = 0; j < sonBrLen.size(); j++) brLen.push_back(sonBrLen[j]);
	}
	return brLen;
}

/******************************************************************************/

double TreeTools::getTotalLength(const Node & node) throw (NodeException)
{
	if(!node.hasDistanceToFather()) throw NodeException("TreeTools::getTotalLength(). No branch length.", &node);
	double length = node.getDistanceToFather();
	for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
		length += getTotalLength(* node.getSon(i));
	}
	return length;
}

/******************************************************************************/

double TreeTools::getTotalLength(const Tree & tree) throw (NodeException)
{
	return getTotalLength(* tree.getRootNode());
}

/******************************************************************************/

void TreeTools::setBranchLengths(Node & node, double brLen) {
	node.setDistanceToFather(brLen);
	for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
		setBranchLengths(* node.getSon(i), brLen);
	}
}

/******************************************************************************/

void TreeTools::setBranchLengths(Tree & tree, double brLen) {
	setBranchLengths(* tree.getRootNode(), brLen);
}

/******************************************************************************/

void TreeTools::setVoidBranchLengths(Node & node, double brLen) {
	if(!node.hasDistanceToFather()) node.setDistanceToFather(brLen);
	for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
		setVoidBranchLengths(* node.getSon(i), brLen);
	}
}

/******************************************************************************/

void TreeTools::setVoidBranchLengths(Tree & tree, double brLen) {
	setVoidBranchLengths(* tree.getRootNode(), brLen);
}

/******************************************************************************/

void TreeTools::scaleTree(Node & node, double factor) throw (NodeException)
{
	if(!node.hasDistanceToFather()) throw NodeException("TreeTools::scaleTree(). Branch with no length", &node);
	node.setDistanceToFather(node.getDistanceToFather() * factor);
	for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
		scaleTree(* node.getSon(i), factor);
	}
}
		
/******************************************************************************/

void TreeTools::scaleTree(Tree & tree, double factor) throw (NodeException)
{
	scaleTree(* tree.getRootNode(), factor);
}

/******************************************************************************/

Tree * TreeTools::getRandomTree(vector<string> & leavesNames)
{
  if(leavesNames.size() == 0) return NULL; // No taxa.
  // This vector will contain all nodes.
  // Start with all leaves, and then group nodes randomly 2 by 2.
  // Att the end, contains only the root node of the tree.
	vector<Node *> nodes(leavesNames.size());
	// Create all leaves nodes:
	for(unsigned int i = 0; i < leavesNames.size(); i++) {
		nodes[i] = new Node(leavesNames[i]);
	}
	// Now group all nodes:
	while(nodes.size() > 1) {
		// Select random nodes:
		int pos1 = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(nodes.size());
		Node * node1 = nodes[pos1];
		nodes.erase(nodes.begin() + pos1);
		int pos2 = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(nodes.size());
		Node * node2 = nodes[pos2];
		nodes.erase(nodes.begin() + pos2);
		// Add new node:
		Node * parent = new Node();
		parent -> addSon(* node1);
		parent -> addSon(* node2);
		nodes.push_back(parent);
	}
  // Return tree with last node as root node:
  return new Tree(* nodes[0]);
}

/******************************************************************************/

// @TODO
//void TreeTools::displayTree(const Node * node)
//{
//
//}

/******************************************************************************/

string TreeTools::BOOTSTRAP = "bootstrap";

/******************************************************************************/
