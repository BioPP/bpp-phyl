//
// File: TreeTemplateTools.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct  13 13:00 2006
// From file TreeTools.cpp
// Created on: Wed Aug  6 13:45:28 2003
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

#include "TreeTemplateTools.h"
#include "TreeTemplate.h"

// From Utils:
#include <Utils/Number.h>
#include <Utils/BppString.h>
#include <Utils/StringTokenizer.h>
#include <Utils/TextTools.h>

// From NumCalc:
#include <NumCalc/RandomTools.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <sstream>

using namespace std;

/******************************************************************************/

bool TreeTemplateTools::isMultifurcating(const Node & node)
{
  if(node.getNumberOfSons() > 2) return true;
  else {
    bool b = false;
    for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
    {
      b = b || isMultifurcating(* node.getSon(i));
    }
    return b;
  }    
}

/******************************************************************************/

unsigned int TreeTemplateTools::getNumberOfLeaves(const Node & node)
{
  unsigned int nbLeaves = 0;
  if(node.isLeaf())
  {
    nbLeaves++;
  } 
  for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    nbLeaves += getNumberOfLeaves(* node[i]);
  }
  return nbLeaves;
}

/******************************************************************************/

unsigned int TreeTemplateTools::getNumberOfNodes(const Node & node)
{
  unsigned int nbNodes = 1;
  for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    nbNodes += getNumberOfNodes(* node[i]);
  }
  return nbNodes;
}

/******************************************************************************/

vector<string> TreeTemplateTools::getLeavesNames(const Node & node)
{
  vector<string> names;
  if(node.isLeaf())
  {
    names.push_back(node.getName());
  }
  for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    vector<string> subNames = getLeavesNames(* node.getSon(i));
    for(unsigned int j = 0; j < subNames.size(); j++) names.push_back(subNames[j]);
  }
  return names;   
}

/******************************************************************************/

unsigned int TreeTemplateTools::getDepth(const Node & node)
{
  unsigned int d = 0;
  for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    unsigned int c = getDepth(* node[i]) + 1;
    if( c > d) d = c;
  }
  return d;
}

/******************************************************************************/

double TreeTemplateTools::getHeight(const Node & node) throw (NodeException)
{
  double d = 0;
  for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    const Node * son = node[i];
    double dist = 0;
    if(son->hasDistanceToFather()) dist = son->getDistanceToFather();
    else throw NodeException("Node without branch length.", son);
    double c = getHeight(* son) + dist;
    if(c > d) d = c;
  }
  return d;
}

/******************************************************************************/

double TreeTemplateTools::getHeights(const Node & node, map<const Node *, double> & heights) throw (NodeException)
{
  double d = 0;
  for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    const Node * son = node[i];
    double dist = 0;
    if(son->hasDistanceToFather()) dist = son->getDistanceToFather();
    else throw NodeException("Node without branch length.", son);
    double c = getHeights(* son, heights) + dist;
    if(c > d) d = c;
  }
  heights[&node] = d;
  return d;
}

/******************************************************************************/

Node * TreeTemplateTools::parenthesisToNode(const string & description, bool bootstrap, const string & propertyName)
{
  //cout << "NODE: " << description << endl;
  TreeTools::Element elt = TreeTools::getElement(description);

  //New node:
  Node * node = new Node();
  if(!TextTools::isEmpty(elt.length))
  {
    node->setDistanceToFather(TextTools::toDouble(elt.length));
    //cout << "NODE: LENGTH: " << * elt.length << endl;
  }
  if(!TextTools::isEmpty(elt.bootstrap))
  {
    if(bootstrap)
    {
      node->setBranchProperty(TreeTools::BOOTSTRAP, Number<double>(TextTools::toDouble(elt.bootstrap)));
      //cout << "NODE: BOOTSTRAP: " << * elt.bootstrap << endl;
    }
    else
    {
      node->setBranchProperty(propertyName, String(elt.bootstrap));
    }
  }
  
  NodeTokenizer nt(elt.content);
  vector<string> elements;
  while(nt.hasNext())
  {
    elements.push_back(nt.next());
  }

  if(elements.size() == 1)
  {
    //This is a leaf:
    //cout << "NODE: LEAF: " << elements[0] << endl;
    string name = TextTools::removeSurroundingWhiteSpaces(elements[0]);
    node->setName(name);
  }
  else
  {
    //This is a node:
    for(unsigned int i = 0; i < elements.size(); i++)
    {
      //cout << "NODE: SUBNODE: " << i << ", " << elements[i] << endl;
      Node * son = parenthesisToNode(elements[i], bootstrap, propertyName);
      node->addSon(* son);
    }
  }
  return node;
}

/******************************************************************************/

TreeTemplate<Node> * TreeTemplateTools::parenthesisToTree(const string & description, bool bootstrap, const string & propertyName) throw (Exception)
{
  string::size_type lastP  = description.rfind(')');
  if(lastP == string::npos)
    throw Exception("TreeTemplateTools::parenthesisToTree(). Bad format: no closing parenthesis found.");
  string::size_type firstP = description.find('(');
  if(firstP == string::npos)
    throw Exception("TreeTemplateTools::parenthesisToTree(). Bad format: no opening parenthesis found.");
  string::size_type semi = description.rfind(';');
  if(semi == string::npos)
    throw Exception("TreeTemplateTools::parenthesisToTree(). Bad format: no semi-colon found.");
  if(lastP <= firstP)
    throw Exception("TreeTemplateTools::parenthesisToTree(). Bad format: closing parenthesis before opening parenthesis.");
  string content = description.substr(firstP + 1, lastP - firstP - 1);
  string element = semi == string::npos ? description.substr(lastP + 1) : description.substr(lastP + 1, semi - lastP - 1);
  //cout << "TREE: " << content << endl;
  //New root node:
  Node * node = new Node();
  
  NodeTokenizer nt(content);
  vector<string> elements;
  while(nt.hasNext())
  {
    elements.push_back(nt.next());
  }

  if(elements.size() == 1)
  {
    //This is a leaf:
    node->setName(elements[0]);
  }
  else
  {
    //This is a node:
    for(unsigned int i = 0; i < elements.size(); i++)
    {
      Node * son = parenthesisToNode(elements[i], bootstrap, propertyName);
      node->addSon(* son);
    }
    if(! TextTools::isEmpty(element))
    {
      StringTokenizer st(element, ":");
      string lengthS = "";
      string bootstrapS = "";
      if(st.numberOfRemainingTokens() == 1)
      { 
        bootstrapS=st.nextToken();
      }
      else
      {
        bootstrapS=st.nextToken();
        lengthS=st.nextToken();
      }
      if(!TextTools::isEmpty(lengthS))
      {
        node->setDistanceToFather(TextTools::toDouble(lengthS));
        //cout << "NODE: LENGTH: " << * elt.length << endl;
      }
      if(!TextTools::isEmpty(bootstrapS))
      {
        if(bootstrap)
        {
          node->setBranchProperty(TreeTools::BOOTSTRAP, Number<double>(TextTools::toDouble(bootstrapS)));
          //cout << "NODE: BOOTSTRAP: " << * elt.bootstrap << endl;
        }
        else
        {
          node->setBranchProperty(propertyName, String(bootstrapS));
        }
      }
    }
  }
  TreeTemplate<Node> * tree = new TreeTemplate<Node>();
  tree->setRootNode(* node);
  tree->resetNodesId();
  return tree;
}

/******************************************************************************/

string TreeTemplateTools::nodeToParenthesis(const Node & node, bool writeId)
{
  ostringstream s;
  if(node.isLeaf())
  {
    s << node.getName();
  }
  else
  {
    s << "(";
    s << nodeToParenthesis(* node[0], writeId);
    for(unsigned int i = 1; i < node.getNumberOfSons(); i++)
    {
      s << "," << nodeToParenthesis(* node[i], writeId);
    }
    s << ")";
  }
  if(writeId)
  {
    if(node.isLeaf()) s << "_";
    s << node.getId();
  }
  else
  {
    if(node.hasBranchProperty(TreeTools::BOOTSTRAP))
      s << (dynamic_cast<const Number<double> *>(node.getBranchProperty(TreeTools::BOOTSTRAP))->getValue());
  }
  if(node.hasDistanceToFather()) s << ":" << node.getDistanceToFather();
  return s.str();  
}

/******************************************************************************/

string TreeTemplateTools::nodeToParenthesis(const Node & node, bool bootstrap, const string & propertyName)
{
  ostringstream s;
  if(node.isLeaf())
  {
    s << node.getName();
  }
  else
  {
    s << "(";
    s << nodeToParenthesis(* node[0], bootstrap, propertyName);
    for(unsigned int i = 1; i < node.getNumberOfSons(); i++)
    {
      s << "," << nodeToParenthesis(* node[i], bootstrap, propertyName);
    }
    s << ")";
  
    if(bootstrap)
    {
      if(node.hasBranchProperty(TreeTools::BOOTSTRAP))
        s << (dynamic_cast<const Number<double> *>(node.getBranchProperty(TreeTools::BOOTSTRAP))->getValue());
    }
    else
    {
      if(node.hasBranchProperty(propertyName))
        s << *(dynamic_cast<const String *>(node.getBranchProperty(propertyName)));
    }
  }
  if(node.hasDistanceToFather()) s << ":" << node.getDistanceToFather();
  return s.str();  
}

/******************************************************************************/

string TreeTemplateTools::treeToParenthesis(const TreeTemplate<Node> & tree, bool writeId)
{
  ostringstream s;
  s << "(";
  const Node * node = tree.getRootNode();
  if(node->isLeaf())
  {
    s << node->getName();
    for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
    {
      s << "," << nodeToParenthesis(* node->getSon(i), writeId);
    }
  }
  else
  {
    s << nodeToParenthesis(* node->getSon(0), writeId);
    for(unsigned int i = 1; i < node->getNumberOfSons(); i++)
    {
      s << "," << nodeToParenthesis(* node->getSon(i), writeId);
    }
  }
  s << ");" << endl;
  return s.str();  
}

/******************************************************************************/

string TreeTemplateTools::treeToParenthesis(const TreeTemplate<Node> & tree, bool bootstrap, const string & propertyName)
{
  ostringstream s;
  s << "(";
  const Node * node = tree.getRootNode();
  if(node->isLeaf())
  {
    s << node -> getName();
    for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
    {
      s << "," << nodeToParenthesis(* node->getSon(i), bootstrap, propertyName);
    }
  }
  else
  {
    s << nodeToParenthesis(* node->getSon(0), bootstrap, propertyName);
    for(unsigned int i = 1; i < node->getNumberOfSons(); i++)
    {
      s << "," << nodeToParenthesis(* node->getSon(i), bootstrap, propertyName);
    }
  }
  s << ")";
  if(bootstrap)
  {
    if(node->hasBranchProperty(TreeTools::BOOTSTRAP))
      s << (dynamic_cast<const Number<double> *>(node->getBranchProperty(TreeTools::BOOTSTRAP))->getValue());
  }
  else
  {
    if(node->hasBranchProperty(propertyName))
      s << *(dynamic_cast<const String *>(node->getBranchProperty(propertyName)));
  }
  s << ";" << endl;
  return s.str();  
}

/******************************************************************************/

Vdouble TreeTemplateTools::getBranchLengths(const Node & node) throw (NodeException)
{
  Vdouble brLen(1);
  if(node.hasDistanceToFather()) brLen[0] = node.getDistanceToFather();
  else throw NodeException("TreeTools::getbranchLengths(). No branch length.", &node);
  for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    Vdouble sonBrLen = getBranchLengths(* node.getSon(i));
    for(unsigned int j = 0; j < sonBrLen.size(); j++) brLen.push_back(sonBrLen[j]);
  }
  return brLen;
}

/******************************************************************************/

double TreeTemplateTools::getTotalLength(const Node & node, bool includeAncestor) throw (NodeException)
{
  if(includeAncestor && !node.hasDistanceToFather()) throw NodeException("TreeTools::getTotalLength(). No branch length.", &node);
  double length = includeAncestor ? node.getDistanceToFather() : 0;
  for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    length += getTotalLength(* node.getSon(i), true);
  }
  return length;
}

/******************************************************************************/

void TreeTemplateTools::setBranchLengths(Node & node, double brLen)
{
  node.setDistanceToFather(brLen);
  for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    setBranchLengths(* node.getSon(i), brLen);
  }
}

/******************************************************************************/

void TreeTemplateTools::deleteBranchLengths(Node & node)
{
  node.deleteDistanceToFather();
  for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    deleteBranchLengths(* node.getSon(i));
  }
}

/******************************************************************************/

void TreeTemplateTools::setVoidBranchLengths(Node & node, double brLen)
{
  if(!node.hasDistanceToFather()) node.setDistanceToFather(brLen);
  for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    setVoidBranchLengths(* node.getSon(i), brLen);
  }
}

/******************************************************************************/

void TreeTemplateTools::scaleTree(Node & node, double factor) throw (NodeException)
{
  if(node.hasFather())
  {
    if(!node.hasDistanceToFather()) throw NodeException("TreeTemplateTools::scaleTree(). Branch with no length", &node);
    node.setDistanceToFather(node.getDistanceToFather() * factor);
  }
  for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    scaleTree(* node.getSon(i), factor);
  }
}
    
/******************************************************************************/

TreeTemplate<Node> * TreeTemplateTools::getRandomTree(vector<string> & leavesNames)
{
  if(leavesNames.size() == 0) return NULL; // No taxa.
  // This vector will contain all nodes.
  // Start with all leaves, and then group nodes randomly 2 by 2.
  // Att the end, contains only the root node of the tree.
  vector<Node *> nodes(leavesNames.size());
  // Create all leaves nodes:
  for(unsigned int i = 0; i < leavesNames.size(); i++)
  {
    nodes[i] = new Node(leavesNames[i]);
  }
  // Now group all nodes:
  while(nodes.size() > 1)
  {
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
  TreeTemplate<Node> * tree = new TreeTemplate<Node>(* nodes[0]);
  tree->resetNodesId();
  return tree;
}

/******************************************************************************/

vector<Node *> TreeTemplateTools::getPathBetweenAnyTwoNodes(Node & node1, Node & node2, bool includeAncestor)
{
  vector<Node *> path;
  vector<Node *> pathMatrix1;
  vector<Node *> pathMatrix2;

  Node * nodeUp = & node1;
  while(nodeUp -> hasFather())  { // while(nodeUp != root)
    pathMatrix1.push_back(nodeUp);
    nodeUp = nodeUp -> getFather();
  }
  pathMatrix1.push_back(nodeUp); // The root.

  nodeUp = & node2;
  while(nodeUp -> hasFather())  {
    pathMatrix2.push_back(nodeUp);
    nodeUp = nodeUp -> getFather();
  }
  pathMatrix2.push_back(nodeUp); // The root.
  // Must check that the two nodes have the same root!!!

  int tmp1 = pathMatrix1.size() - 1;
  int tmp2 = pathMatrix2.size() - 1;

  while((tmp1 >= 0) && (tmp2 >= 0)) {
    if (pathMatrix1[tmp1] != pathMatrix2[tmp2]) break;
    tmp1--; tmp2--;
  }

  for (int y = 0; y <= tmp1; ++y) path.push_back(pathMatrix1[y]);
  if(includeAncestor) path.push_back(pathMatrix1[tmp1 + 1]); // pushing once, the Node that was common to both.
  for (int j = tmp2; j >= 0; --j) {
    path.push_back(pathMatrix2[j]);
  }
  return path;
}

/******************************************************************************/

vector<const Node *> TreeTemplateTools::getPathBetweenAnyTwoNodes(const Node & node1, const Node & node2, bool includeAncestor)
{
  vector<const Node *> path;
  vector<const Node *> pathMatrix1;
  vector<const Node *> pathMatrix2;

  const Node * nodeUp = & node1;
  while(nodeUp -> hasFather())  { // while(nodeUp != root)
    pathMatrix1.push_back(nodeUp);
    nodeUp = nodeUp -> getFather();
  }
  pathMatrix1.push_back(nodeUp); // The root.

  nodeUp = & node2;
  while(nodeUp -> hasFather())  {
    pathMatrix2.push_back(nodeUp);
    nodeUp = nodeUp -> getFather();
  }
  pathMatrix2.push_back(nodeUp); // The root.
  // Must check that the two nodes have the same root!!!

  int tmp1 = pathMatrix1.size() - 1;
  int tmp2 = pathMatrix2.size() - 1;

  while((tmp1 >= 0) && (tmp2 >= 0)) {
    if (pathMatrix1[tmp1] != pathMatrix2[tmp2]) break;
    tmp1--; tmp2--;
  }

  for (int y = 0; y <= tmp1; ++y) path.push_back(pathMatrix1[y]);
  if(includeAncestor) path.push_back(pathMatrix1[tmp1 + 1]); // pushing once, the Node that was common to both.
  for (int j = tmp2; j >= 0; --j) {
    path.push_back(pathMatrix2[j]);
  }
  return path;
}

/******************************************************************************/

double TreeTemplateTools::getDistanceBetweenAnyTwoNodes(const Node & node1, const Node & node2)
{
  vector<const Node *> path = getPathBetweenAnyTwoNodes(node1, node2, false);
  double d = 0;
  for(unsigned int i = 0; i < path.size(); i++)
  {
    d += path[i]->getDistanceToFather();
  }
  return d;
}
  
/******************************************************************************/

DistanceMatrix * TreeTemplateTools::getDistanceMatrix(const TreeTemplate<Node> & tree)
{
  vector<const Node *> nodes = tree.getLeaves();
  vector<string> names(nodes.size());
  for(unsigned int i = 0; i < nodes.size(); i++)
    names[i] = nodes[i]->getName();
  DistanceMatrix * mat = new DistanceMatrix(names);
  for(unsigned int i = 0; i < nodes.size(); i++)
  {
    (* mat)(i, i) = 0;
    for(unsigned int j = 0; j < i; j++)
    {
      (* mat)(i, j) = (* mat)(j, i) = getDistanceBetweenAnyTwoNodes(*nodes[i], *nodes[j]);
    }
  }
  return mat;
}
/******************************************************************************/

vector<const Node *> TreeTemplateTools::getRemainingNeighbors(const Node * node1, const Node * node2, const Node * node3)
{
  vector<const Node *> neighbors = node1->getNeighbors();
  vector<const Node *> neighbors2;
  for(unsigned int k = 0; k < neighbors.size(); k++)
  {
    const Node * n = neighbors[k];
    if(n != node2 && n != node3) neighbors2.push_back(n);
  }
  return neighbors2;
}

/******************************************************************************/

void TreeTemplateTools::incrementAllIds(Node * node, int increment)
{
  node->setId(node->getId() + increment);
  for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    incrementAllIds(node->getSon(i), increment);
  }
}

/******************************************************************************/

