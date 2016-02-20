//
// File: TreeTemplateTools.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct  13 13:00 2006
// From file TreeTools.cpp
// Created on: Wed Aug  6 13:45:28 2003
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

#include "TreeTemplateTools.h"
#include "TreeTemplate.h"

#include <Bpp/Numeric/Number.h>
#include <Bpp/BppString.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/NestedStringTokenizer.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <sstream>
#include <limits>

using namespace std;

/******************************************************************************/

bool TreeTemplateTools::isMultifurcating(const Node& node)
{
  if (node.getNumberOfSons() > 2)
    return true;
  else
  {
    bool b = false;
    for (size_t i = 0; i < node.getNumberOfSons(); i++)
    {
      b = b || isMultifurcating(*node.getSon(i));
    }
    return b;
  }
}

/******************************************************************************/

size_t TreeTemplateTools::getNumberOfLeaves(const Node& node)
{
  size_t nbLeaves = 0;
  if (node.isLeaf())
  {
    nbLeaves++;
  }
  for (int i = 0; i < static_cast<int>(node.getNumberOfSons()); i++)
  {
    nbLeaves += getNumberOfLeaves(*node[i]);
  }
  return nbLeaves;
}

/******************************************************************************/

vector<string> TreeTemplateTools::getLeavesNames(const Node& node)
{
  vector<string> names;
  if (node.isLeaf())
  {
    names.push_back(node.getName());
  }
  for (size_t i = 0; i < node.getNumberOfSons(); i++)
  {
    vector<string> subNames = getLeavesNames(*node.getSon(i));
    for (size_t j = 0; j < subNames.size(); j++)
    {
      names.push_back(subNames[j]);
    }
  }
  return names;
}

/******************************************************************************/

unsigned int TreeTemplateTools::getDepth(const Node& node)
{
  unsigned int d = 0;
  for (int i = 0; i < static_cast<int>(node.getNumberOfSons()); i++)
  {
    unsigned int c = getDepth(*node[i]) + 1;
    if (c > d)
      d = c;
  }
  return d;
}

/******************************************************************************/

unsigned int TreeTemplateTools::getDepths(const Node& node, map<const Node*, unsigned int>& depths)
{
  unsigned int d = 0;
  for (int i = 0; i < static_cast<int>(node.getNumberOfSons()); i++)
  {
    unsigned int c = getDepths(*node[i], depths) + 1;
    if (c > d)
      d = c;
  }
  depths[&node] = d;
  return d;
}

/******************************************************************************/

double TreeTemplateTools::getHeight(const Node& node)
{
  double d = 0;
  for (int i = 0; i < static_cast<int>(node.getNumberOfSons()); i++)
  {
    const Node* son = node[i];
    double dist = son->getDistanceToFather();
    double c = getHeight(*son) + dist;
    if (c > d)
      d = c;
  }
  return d;
}

/******************************************************************************/

double TreeTemplateTools::getHeights(const Node& node, map<const Node*, double>& heights)
{
  double d = 0;
  for (int i = 0; i < static_cast<int>(node.getNumberOfSons()); i++)
  {
    const Node* son = node[i];
    double dist = son->getDistanceToFather();
    double c = getHeights(*son, heights) + dist;
    if (c > d)
      d = c;
  }
  heights[&node] = d;
  return d;
}

/******************************************************************************/

TreeTemplateTools::Element TreeTemplateTools::getElement(const string& elt) throw (IOException)
{
  Element element;
  element.length    = ""; // default
  element.bootstrap = ""; // default
  element.isLeaf    = false; // default

  size_t colonIndex;
  bool hasColon = false;
  for (colonIndex = elt.size(); colonIndex > 0 && elt[colonIndex] != ')'; colonIndex--)
  {
    if (elt[colonIndex] == ':')
    {
      hasColon = true;
      break;
    }
  }
  try
  {
    string elt2;
    if (hasColon)
    {
      // this is an element with length:
      elt2 = elt.substr(0, colonIndex);
      element.length = TextTools::removeSurroundingWhiteSpaces(elt.substr(colonIndex + 1));
    }
    else
    {
      // this is an element without length;
      elt2 = elt;
    }

    string::size_type lastP = elt2.rfind(')');
    string::size_type firstP = elt2.find('(');
    if (firstP == string::npos)
    {
      // This is a leaf:
      element.content = elt2;
      element.isLeaf = true;
    }
    else
    {
      // This is a node:
      if (lastP < firstP)
        throw IOException("TreeTemplateTools::getElement(). Invalid format: bad closing parenthesis in " + elt2);
      element.content = TextTools::removeSurroundingWhiteSpaces(elt2.substr(firstP + 1, lastP - firstP - 1));
      string bootstrap = TextTools::removeSurroundingWhiteSpaces(elt2.substr(lastP + 1));
      // cout << "ELEMENT: BOOTSTRAP: " << bootstrap << endl;
      if (!TextTools::isEmpty(bootstrap))
      {
        element.bootstrap = bootstrap;
      }
    }
  }
  catch (exception e)
  {
    throw IOException("Bad tree description: " + elt);
  }
  return element;
}

/******************************************************************************/


Node* TreeTemplateTools::parenthesisToNode(const string& description, unsigned int& nodeCounter, bool bootstrap, const string& propertyName, bool withId, bool verbose)
{
  // cout << "NODE: " << description << endl;
  Element elt = getElement(description);

  // New node:
  Node* node = new Node();
  if (!TextTools::isEmpty(elt.length))
  {
    node->setDistanceToFather(TextTools::toDouble(elt.length));
    // cout << "NODE: LENGTH: " << * elt.length << endl;
  }
  if (!TextTools::isEmpty(elt.bootstrap))
  {
    if (withId)
    {
      node->setId(TextTools::toInt(elt.bootstrap));
    }
    else
    {
      if (bootstrap)
      {
        node->setBranchProperty(TreeTools::BOOTSTRAP, Number<double>(TextTools::toDouble(elt.bootstrap)));
        // cout << "NODE: BOOTSTRAP: " << * elt.bootstrap << endl;
      }
      else
      {
        node->setBranchProperty(propertyName, BppString(elt.bootstrap));
      }
    }
  }

  NestedStringTokenizer nt(elt.content, "(", ")", ",");
  vector<string> elements;
  while (nt.hasMoreToken())
  {
    elements.push_back(nt.nextToken());
  }

  if (elt.isLeaf)
  {
    // This is a leaf:
    string name = TextTools::removeSurroundingWhiteSpaces(elements[0]);
    if (withId)
    {
      StringTokenizer st(name, "_", true, true);
      ostringstream realName;
      for (size_t i = 0; i < st.numberOfRemainingTokens() - 1; ++i)
      {
        if (i != 0)
        {
          realName << "_";
        }
        realName << st.getToken(i);
      }
      node->setName(realName.str());
      node->setId(TextTools::toInt(st.getToken(st.numberOfRemainingTokens() - 1)));
    }
    else
    {
      node->setName(name);
    }
  }
  else
  {
    // This is a node:
    for (size_t i = 0; i < elements.size(); i++)
    {
      // cout << "NODE: SUBNODE: " << i << ", " << elements[i] << endl;
      Node* son = parenthesisToNode(elements[i], nodeCounter, bootstrap, propertyName, withId, verbose);
      node->addSon(son);
    }
  }
  nodeCounter++;
  if (verbose)
    ApplicationTools::displayUnlimitedGauge(nodeCounter);
  return node;
}

/******************************************************************************/

TreeTemplate<Node>* TreeTemplateTools::parenthesisToTree(const string& description, bool bootstrap, const string& propertyName, bool withId, bool verbose) throw (Exception)
{
  string::size_type semi = description.rfind(';');
  if (semi == string::npos)
    throw Exception("TreeTemplateTools::parenthesisToTree(). Bad format: no semi-colon found.");
  string content = description.substr(0, semi);
  unsigned int nodeCounter = 0;
  Node* node = parenthesisToNode(content, nodeCounter, bootstrap, propertyName, withId, verbose);
  TreeTemplate<Node>* tree = new TreeTemplate<Node>();
  tree->setRootNode(node);
  if (!withId)
  {
    tree->resetNodesId();
  }
  if (verbose) {
    (*ApplicationTools::message) << " nodes loaded.";
    ApplicationTools::message->endLine();
  }
  return tree;
}

/******************************************************************************/

string TreeTemplateTools::nodeToParenthesis(const Node& node, bool writeId)
{
  ostringstream s;
  if (node.isLeaf())
  {
    s << node.getName();
  }
  else
  {
    s << "(";
    s << nodeToParenthesis(*node[0], writeId);
    for (int i = 1; i < static_cast<int>(node.getNumberOfSons()); i++)
    {
      s << "," << nodeToParenthesis(*node[i], writeId);
    }
    s << ")";
  }
  if (writeId)
  {
    if (node.isLeaf())
      s << "_";
    s << node.getId();
  }
  else
  {
    if (node.hasBranchProperty(TreeTools::BOOTSTRAP))
      s << (dynamic_cast<const Number<double>*>(node.getBranchProperty(TreeTools::BOOTSTRAP))->getValue());
  }
  if (node.hasDistanceToFather())
    s << ":" << node.getDistanceToFather();
  return s.str();
}

/******************************************************************************/

string TreeTemplateTools::nodeToParenthesis(const Node& node, bool bootstrap, const string& propertyName)
{
  ostringstream s;
  if (node.isLeaf())
  {
    s << node.getName();
  }
  else
  {
    s << "(";
    s << nodeToParenthesis(*node[0], bootstrap, propertyName);
    for (int i = 1; i < static_cast<int>(node.getNumberOfSons()); i++)
    {
      s << "," << nodeToParenthesis(*node[i], bootstrap, propertyName);
    }
    s << ")";

    if (bootstrap)
    {
      if (node.hasBranchProperty(TreeTools::BOOTSTRAP))
        s << (dynamic_cast<const Number<double>*>(node.getBranchProperty(TreeTools::BOOTSTRAP))->getValue());
    }
    else
    {
      if (node.hasBranchProperty(propertyName))
      {
        const BppString* ppt = dynamic_cast<const BppString*>(node.getBranchProperty(propertyName));
        if (ppt)
          s << *ppt;
        else
          throw Exception("TreeTemplateTools::nodeToParenthesis. Property should be a BppString.");
      }
    }
  }
  if (node.hasDistanceToFather())
    s << ":" << node.getDistanceToFather();
  return s.str();
}

/******************************************************************************/

string TreeTemplateTools::treeToParenthesis(const TreeTemplate<Node>& tree, bool writeId)
{
  ostringstream s;
  s << "(";
  const Node* node = tree.getRootNode();
  if (node->isLeaf() && node->hasName()) // In case we have a tree like ((A:1.0)); where the root node is an unamed leaf!
  {
    s << node->getName();
    for (size_t i = 0; i < node->getNumberOfSons(); ++i)
    {
      s << "," << nodeToParenthesis(*node->getSon(i), writeId);
    }
  }
  else
  {
    s << nodeToParenthesis(*node->getSon(0), writeId);
    for (size_t i = 1; i < node->getNumberOfSons(); ++i)
    {
      s << "," << nodeToParenthesis(*node->getSon(i), writeId);
    }
  }
  s << ")";
  if (node->hasDistanceToFather())
    s << ":" << node->getDistanceToFather();
  s << ";" << endl;
  return s.str();
}

/******************************************************************************/

string TreeTemplateTools::treeToParenthesis(const TreeTemplate<Node>& tree, bool bootstrap, const string& propertyName)
{
  ostringstream s;
  s << "(";
  const Node* node = tree.getRootNode();
  if (node->isLeaf())
  {
    s << node->getName();
    for (size_t i = 0; i < node->getNumberOfSons(); i++)
    {
      s << "," << nodeToParenthesis(*node->getSon(i), bootstrap, propertyName);
    }
  }
  else
  {
    s << nodeToParenthesis(*node->getSon(0), bootstrap, propertyName);
    for (size_t i = 1; i < node->getNumberOfSons(); i++)
    {
      s << "," << nodeToParenthesis(*node->getSon(i), bootstrap, propertyName);
    }
  }
  s << ")";
  if (bootstrap)
  {
    if (node->hasBranchProperty(TreeTools::BOOTSTRAP))
      s << (dynamic_cast<const Number<double>*>(node->getBranchProperty(TreeTools::BOOTSTRAP))->getValue());
  }
  else
  {
    if (node->hasBranchProperty(propertyName))
    {
      const BppString* ppt = dynamic_cast<const BppString*>(node->getBranchProperty(propertyName));
      if (ppt)
        s << *ppt;
      else
        throw Exception("TreeTemplateTools::nodeToParenthesis. Property should be a BppString.");
    }
  }
  s << ";" << endl;
  return s.str();
}

/******************************************************************************/

Vdouble TreeTemplateTools::getBranchLengths(const Node& node) throw (NodePException)
{
  Vdouble brLen(1);
  brLen[0] = node.getDistanceToFather();
  for (size_t i = 0; i < node.getNumberOfSons(); i++)
  {
    Vdouble sonBrLen = getBranchLengths(*node.getSon(i));
    for (size_t j = 0; j < sonBrLen.size(); j++)
    {
      brLen.push_back(sonBrLen[j]);
    }
  }
  return brLen;
}

/******************************************************************************/

double TreeTemplateTools::getTotalLength(const Node& node, bool includeAncestor) throw (NodePException)
{
  if (includeAncestor && !node.hasDistanceToFather())
    throw NodePException("TreeTools::getTotalLength(). No branch length.", &node);
  double length = includeAncestor ? node.getDistanceToFather() : 0;
  for (size_t i = 0; i < node.getNumberOfSons(); i++)
  {
    length += getTotalLength(*node.getSon(i), true);
  }
  return length;
}

/******************************************************************************/

void TreeTemplateTools::setBranchLengths(Node& node, double brLen)
{
  node.setDistanceToFather(brLen);
  for (size_t i = 0; i < node.getNumberOfSons(); i++)
  {
    setBranchLengths(*node.getSon(i), brLen);
  }
}

/******************************************************************************/

void TreeTemplateTools::deleteBranchLengths(Node& node)
{
  node.deleteDistanceToFather();
  for (size_t i = 0; i < node.getNumberOfSons(); i++)
  {
    deleteBranchLengths(*node.getSon(i));
  }
}

/******************************************************************************/

void TreeTemplateTools::setVoidBranchLengths(Node& node, double brLen)
{
  if (!node.hasDistanceToFather())
    node.setDistanceToFather(brLen);
  for (size_t i = 0; i < node.getNumberOfSons(); i++)
  {
    setVoidBranchLengths(*node.getSon(i), brLen);
  }
}

/******************************************************************************/

void TreeTemplateTools::scaleTree(Node& node, double factor) throw (NodePException)
{
  if (node.hasFather())
  {
    node.setDistanceToFather(node.getDistanceToFather() * factor);
  }
  for (size_t i = 0; i < node.getNumberOfSons(); i++)
  {
    scaleTree(*node.getSon(i), factor);
  }
}

/******************************************************************************/

TreeTemplate<Node>* TreeTemplateTools::getRandomTree(vector<string>& leavesNames, bool rooted)
{
  if (leavesNames.size() == 0)
    return 0;                                // No taxa.
  // This vector will contain all nodes.
  // Start with all leaves, and then group nodes randomly 2 by 2.
  // Att the end, contains only the root node of the tree.
  vector<Node*> nodes(leavesNames.size());
  // Create all leaves nodes:
  for (size_t i = 0; i < leavesNames.size(); ++i)
  {
    nodes[i] = new Node(leavesNames[i]);
  }
  // Now group all nodes:
  while (nodes.size() > (rooted ? 2 : 3))
  {
    // Select random nodes:
    size_t pos1 = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(nodes.size());
    Node* node1 = nodes[pos1];
    nodes.erase(nodes.begin() + static_cast<ptrdiff_t>(pos1));
    size_t pos2 = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(nodes.size());
    Node* node2 = nodes[pos2];
    nodes.erase(nodes.begin() + static_cast<ptrdiff_t>(pos2));
    // Add new node:
    Node* parent = new Node();
    parent->addSon(node1);
    parent->addSon(node2);
    nodes.push_back(parent);
  }
  // Return tree with last node as root node:
  Node* root = new Node();
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    root->addSon(nodes[i]);
  }
  TreeTemplate<Node>* tree = new TreeTemplate<Node>(root);
  tree->resetNodesId();
  return tree;
}

/******************************************************************************/

vector<Node*> TreeTemplateTools::getPathBetweenAnyTwoNodes(Node& node1, Node& node2, bool includeAncestor, bool includeAncestorAtEndOfPath)
{
  vector<Node*> path;
  vector<Node*> pathMatrix1;
  vector<Node*> pathMatrix2;

  Node* nodeUp = &node1;
  while (nodeUp->hasFather())   // while(nodeUp != root)
  {
    pathMatrix1.push_back(nodeUp);
    nodeUp = nodeUp->getFather();
  }
  pathMatrix1.push_back(nodeUp); // The root.

  nodeUp = &node2;
  while (nodeUp->hasFather())
  {
    pathMatrix2.push_back(nodeUp);
    nodeUp = nodeUp->getFather();
  }
  pathMatrix2.push_back(nodeUp); // The root.

  size_t pos1 = pathMatrix1.size() - 1;
  size_t pos2 = pathMatrix2.size() - 1;
  // Must check that the two nodes have the same root!!!
  if (pathMatrix1[pos1] != pathMatrix2[pos2])
    throw Exception("TreeTemplateTools::getPathBetweenAnyTwoNodes(). The two nodes do not have any ancestor in common / do not belong to the same tree.");

  if (pos1 == 0 && pos2 == 0) {
    //Node 1 and 2 are the root node!
    path.push_back(pathMatrix1[0]);
  } else if (pos1 == 0) {
    //Node 1 is the root node
    //Note: we need to use push_back here as the insert method does not work with reverse iterators.
    for (size_t i = (includeAncestorAtEndOfPath ? pathMatrix2.size(): pathMatrix2.size() - 1); i > 0; --i)
      path.push_back(pathMatrix2[i-1]);
  } else if (pos2 == 0) {
    //Node 2 is the root node
    path.insert(path.end(), pathMatrix1.begin(), (includeAncestorAtEndOfPath ? pathMatrix1.end() : --pathMatrix1.end()));
  } else {
    Node* commonAnc = 0;
    while (pathMatrix1[pos1] == pathMatrix2[pos2] && pos1 > 0 && pos2 > 0)
    {
      commonAnc = pathMatrix1[pos1];
      pos1--; pos2--;
    }

    path.insert(path.end(), pathMatrix1.begin(), pathMatrix1.begin() + static_cast<ptrdiff_t>(pos1 + 1));
    if (includeAncestor && commonAnc)
      path.push_back(commonAnc); // pushing once the Node that was common to both.
    // If node1 or node2 is the common ancestor, then commonAnc is null
    // and was added as node1 or node2, respectively, if includeAncestorAtEndOfPath was set to true.
    //Note: we need to use push_back here as the insert method does not work with reverse iterators.
    for (size_t i = pos2 + 1; i > 0; --i)
      path.push_back(pathMatrix2[i-1]);
  }
  return path;
}

/******************************************************************************/

vector<const Node*> TreeTemplateTools::getPathBetweenAnyTwoNodes(const Node& node1, const Node& node2, bool includeAncestor, bool includeAncestorAtEndOfPath)
{
  vector<const Node*> path;
  vector<const Node*> pathMatrix1;
  vector<const Node*> pathMatrix2;

  const Node* nodeUp = &node1;
  while (nodeUp->hasFather())   // while(nodeUp != root)
  {
    pathMatrix1.push_back(nodeUp);
    nodeUp = nodeUp->getFather();
  }
  pathMatrix1.push_back(nodeUp); // The root.

  nodeUp = &node2;
  while (nodeUp->hasFather())
  {
    pathMatrix2.push_back(nodeUp);
    nodeUp = nodeUp->getFather();
  }
  pathMatrix2.push_back(nodeUp); // The root.

  size_t pos1 = pathMatrix1.size() - 1;
  size_t pos2 = pathMatrix2.size() - 1;
  // Must check that the two nodes have the same root!!!
  if (pathMatrix1[pos1] != pathMatrix2[pos2])
    throw Exception("TreeTemplateTools::getPathBetweenAnyTwoNodes(). The two nodes do not have any ancestor in common / do not belong to the same tree.");

  if (pos1 == 0 && pos2 == 0) {
    //Node 1 and 2 are the root node!
    path.push_back(pathMatrix1[0]);
  } else if (pos1 == 0) {
    //Node 1 is the root node
    //Note: we need to use push_back here as the insert method does not work with reverse iterators.
    for (size_t i = (includeAncestorAtEndOfPath ? pathMatrix2.size(): pathMatrix2.size() - 1); i > 0; --i)
      path.push_back(pathMatrix2[i-1]);
  } else if (pos2 == 0) {
    //Node 2 is the root node
    path.insert(path.end(), pathMatrix1.begin(), (includeAncestorAtEndOfPath ? pathMatrix1.end() : --pathMatrix1.end()));
  } else {
    const Node* commonAnc = 0;
    while (pathMatrix1[pos1] == pathMatrix2[pos2] && pos1 > 0 && pos2 > 0)
    {
      commonAnc = pathMatrix1[pos1];
      pos1--; pos2--;
    }

    path.insert(path.end(), pathMatrix1.begin(), pathMatrix1.begin() + static_cast<ptrdiff_t>(pos1 + 1));
    if (commonAnc &&includeAncestor)
        path.push_back(commonAnc); // pushing once the Node that was common to both.
    // If node1 or node2 is the common ancestor, then commonAnc is null
    // and was added as node1 or node2, respectively, if includeAncestorAtEndOfPath was set to true.
    
    //Note: we need to use push_back here as the insert method does not work with reverse iterators.
    for (size_t i = pos2 + 1; i > 0; --i)
      path.push_back(pathMatrix2[i-1]);
  }
  return path;
}

/******************************************************************************/

double TreeTemplateTools::getDistanceBetweenAnyTwoNodes(const Node& node1, const Node& node2)
{
  vector<const Node*> path = getPathBetweenAnyTwoNodes(node1, node2, false, false);
  double d = 0;
  for (size_t i = 0; i < path.size(); ++i)
  {
    d += path[i]->getDistanceToFather();
  }
  return d;
}

/******************************************************************************/

void TreeTemplateTools::processDistsInSubtree_(const Node* node, DistanceMatrix& matrix, vector< std::pair<string, double> >& distsToNodeFather)
{
  distsToNodeFather.clear();

  // node-is-leaf case
  if (node->getNumberOfSons() == 0)
  {
    distsToNodeFather.push_back(make_pair(node->getName(), node->getDistanceToFather()));
    return;
  }

  // For all leaves in node's subtree, get leaf-to-node distances.
  // Leaves are classified upon node's sons.
  map<const Node*, vector< pair<string, double> > > leavesDists;
  for (size_t i = 0; i < node->getNumberOfSons(); ++i)
  {
    const Node* son = node->getSon(i);
    processDistsInSubtree_(son, matrix, leavesDists[son]); // recursivity
  }
  // Write leaf-leaf distances to the distance matrix.
  // Only pairs in which the two leaves belong to different
  // sons are considered.
  for (size_t son1_loc = 0; son1_loc < node->getNumberOfSons(); ++son1_loc)
  {
    for (size_t son2_loc = 0; son2_loc < son1_loc; ++son2_loc)
    {
      const Node* son1 = node->getSon(son1_loc);
      const Node* son2 = node->getSon(son2_loc);

      for (vector< pair<string, double> >::iterator son1_leaf = leavesDists[son1].begin();
           son1_leaf != leavesDists[son1].end();
           ++son1_leaf)
      {
        for (vector< pair<string, double> >::iterator son2_leaf = leavesDists[son2].begin();
             son2_leaf != leavesDists[son2].end();
             ++son2_leaf)
        {
          matrix(son1_leaf->first, son2_leaf->first) =
            matrix(son2_leaf->first, son1_leaf->first) =
              ( son1_leaf->second + son2_leaf->second );
        }
      }
    }
  }

  // node-is-root case
  if (!node->hasFather())
  {
    // node-is-root-and-leaf case
    if (node->isLeaf() )
    {
      string root_name = node->getName();
      for (vector< pair<string, double> >::iterator other_leaf = leavesDists[node->getSon(0)].begin();
           other_leaf != leavesDists[node->getSon(0)].end();
           ++other_leaf)
      {
        matrix(root_name, other_leaf->first) = matrix( other_leaf->first, root_name) = other_leaf->second;
      }
    }

    return;
  }

  // Get distances from node's father to considered leaves
  distsToNodeFather.clear();
  double nodeToFather = node->getDistanceToFather();
  for (map<const Node*, vector<pair<string, double> > >::iterator son = leavesDists.begin(); son != leavesDists.end(); ++son)
  {
    for (vector< pair<string, double> >::iterator leaf = (son->second).begin(); leaf != (son->second).end(); ++leaf)
    {
      distsToNodeFather.push_back(make_pair(leaf->first, (leaf->second + nodeToFather)));
    }
  }
}

DistanceMatrix* TreeTemplateTools::getDistanceMatrix(const TreeTemplate<Node>& tree)
{
  DistanceMatrix* matrix = new DistanceMatrix(tree.getLeavesNames());
  vector< pair<string, double> > distsToRoot;
  processDistsInSubtree_(tree.getRootNode(), *matrix, distsToRoot);
  return matrix;
}

/******************************************************************************/

std::vector<const Node*> TreeTemplateTools::getRemainingNeighbors(const Node* node1, const Node* node2, const Node* node3)
{
  vector<const Node*> neighbors = node1->getNeighbors();
  vector<const Node*> neighbors2;
  for (size_t k = 0; k < neighbors.size(); k++)
  {
    const Node* n = neighbors[k];
    if (n != node2 && n != node3)
      neighbors2.push_back(n);
  }
  return neighbors2;
}

/******************************************************************************/

void TreeTemplateTools::incrementAllIds(Node* node, int increment)
{
  node->setId(node->getId() + increment);
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    incrementAllIds(node->getSon(i), increment);
  }
}

/******************************************************************************/

void TreeTemplateTools::getNodePropertyNames(const Node& node, vector<string>& propertyNames)
{
  VectorTools::extend(propertyNames, node.getNodePropertyNames());
  for (size_t i = 0; i < node.getNumberOfSons(); i++)
  {
    getNodePropertyNames(*node.getSon(i), propertyNames);
  }
}

void TreeTemplateTools::getNodeProperties(const Node& node, const string& propertyName, map<int, const Clonable*>& properties)
{
  if (node.hasNodeProperty(propertyName))
    properties[node.getId()] = node.getNodeProperty(propertyName);
  for (size_t i = 0; i < node.getNumberOfSons(); i++)
  {
    getNodeProperties(*node.getSon(i), propertyName, properties);
  }
}

void TreeTemplateTools::getNodeProperties(Node& node, const string& propertyName, map<int, Clonable*>& properties)
{
  if (node.hasNodeProperty(propertyName))
    properties[node.getId()] = node.getNodeProperty(propertyName);
  for (size_t i = 0; i < node.getNumberOfSons(); i++)
  {
    getNodeProperties(*node.getSon(i), propertyName, properties);
  }
}

/******************************************************************************/

void TreeTemplateTools::getBranchPropertyNames(const Node& node, vector<string>& propertyNames)
{
  VectorTools::extend(propertyNames, node.getBranchPropertyNames());
  for (size_t i = 0; i < node.getNumberOfSons(); i++)
  {
    getBranchPropertyNames(*node.getSon(i), propertyNames);
  }
}

void TreeTemplateTools::getBranchProperties(const Node& node, const string& propertyName, map<int, const Clonable*>& properties)
{
  if (node.hasBranchProperty(propertyName))
    properties[node.getId()] = node.getBranchProperty(propertyName);
  for (size_t i = 0; i < node.getNumberOfSons(); i++)
  {
    getBranchProperties(*node.getSon(i), propertyName, properties);
  }
}

void TreeTemplateTools::getBranchProperties(Node& node, const string& propertyName, map<int, Clonable*>& properties)
{
  if (node.hasBranchProperty(propertyName))
    properties[node.getId()] = node.getBranchProperty(propertyName);
  for (size_t i = 0; i < node.getNumberOfSons(); i++)
  {
    getBranchProperties(*node.getSon(i), propertyName, properties);
  }
}

/******************************************************************************/

bool TreeTemplateTools::haveSameOrderedTopology(const Node& n1, const Node& n2)
{
  if (n1.isLeaf() && n2.isLeaf() && n1.getName() != n2.getName())
    return false;
  size_t nl1 = n1.getNumberOfSons();
  size_t nl2 = n2.getNumberOfSons();
  if (nl1 != nl2)
    return false;

  bool test = true;
  for (size_t i = 0; test && i < n1.getNumberOfSons(); ++i)
  {
    test &= haveSameOrderedTopology(*n1.getSon(i), *n2.getSon(i));
  }
  return test;
}

/******************************************************************************/

TreeTemplateTools::OrderTreeData_ TreeTemplateTools::orderTree_(Node& node, bool downward, bool orderLeaves)
{
  OrderTreeData_ otd;

  if (node.isLeaf() && node.hasFather())
  {
    otd.size = 1;
    otd.firstLeaf = node.getName();
  }
  else
  {
    vector<size_t> nbSons;
    vector<string> firstLeaves;
    for (size_t i = 0; i < node.getNumberOfSons(); i++)
    {
      OrderTreeData_ otdsub = orderTree_(*node.getSon(i), downward, orderLeaves);
      if (i == 0)
        otd.firstLeaf = otdsub.firstLeaf;
      else if (orderLeaves && otdsub.firstLeaf < otd.firstLeaf)
        otd.firstLeaf = otdsub.firstLeaf;
      nbSons.push_back(otdsub.size);
      firstLeaves.push_back(otdsub.firstLeaf);
    }
    otd.size = VectorTools::sum(nbSons);

    // Now swap nodes:
    if (downward)
    {
      for (size_t i = 0; i < nbSons.size() - 1; ++i)
      {
        size_t pos;
        vector<size_t> index = VectorTools::whichMaxAll(nbSons);
        if (index.size() == 1 || !orderLeaves)
        {
          pos = index[0];
        }
        else
        {
          // There are ties to solve:
          vector<string> v;
          for (size_t j = 0; j < index.size(); ++j)
          {
            v.push_back(firstLeaves[index[j]]);
          }
          size_t mx = VectorTools::whichMax(v);
          pos = index[mx];
        }
        if (pos != i)
        {
          node.swap(i, pos);
          nbSons[pos] = nbSons[i];
        }
        nbSons[i] = 0;
      }
    }
    else
    {
      for (size_t i = 0; i < nbSons.size() - 1; ++i)
      {
        size_t pos;
        vector<size_t> index = VectorTools::whichMinAll(nbSons);
        if (index.size() == 1 || !orderLeaves)
        {
          pos = index[0];
        }
        else
        {
          // There are ties to solve:
          vector<string> v;
          for (size_t j = 0; j < index.size(); ++j)
          {
            v.push_back(firstLeaves[index[j]]);
          }
          size_t mx = VectorTools::whichMin(v);
          pos = index[mx];
        }
        if (pos != i)
        {
          node.swap(i, pos);
          nbSons[pos] = nbSons[i];
        }
        nbSons[i] = otd.size + 1;
      }
    }
  }
  return otd;
}

/******************************************************************************/
const short TreeTemplateTools::MIDROOT_VARIANCE = 0;
const short TreeTemplateTools::MIDROOT_SUM_OF_SQUARES = 1;

void TreeTemplateTools::midRoot(TreeTemplate<Node>& tree, short criterion, bool forceBranchRoot)
{
  if (criterion != MIDROOT_VARIANCE && criterion != MIDROOT_SUM_OF_SQUARES)
    throw Exception("TreeTemplateTools::midRoot(). Illegal criterion value '" + TextTools::toString(criterion) + "'");

  if (tree.isRooted())
    tree.unroot();
  Node* ref_root = tree.getRootNode();
  //
  // The bestRoot object records :
  // -- the current best branch : .first
  // -- the current best value of the criterion : .second["value"]
  // -- the best position of the root on the branch : .second["position"]
  //      0 is toward the original root, 1 is away from it
  //
  pair<Node*, map<string, double> > best_root_branch;
  best_root_branch.first = ref_root; // nota: the root does not correspond to a branch as it has no father
  best_root_branch.second ["position"] = -1;
  best_root_branch.second ["score"] = numeric_limits<double>::max();

  // find the best root
  getBestRootInSubtree_(tree, criterion, ref_root, best_root_branch);
  tree.rootAt(ref_root); // back to the original root

  // reroot
  const double pos = best_root_branch.second["position"];
  if (pos < 1e-6 or pos > 1 - 1e-6)
    // The best root position is on a node (this is often the case with the sum of squares criterion)
    tree.rootAt(pos < 1e-6 ? best_root_branch.first->getFather() : best_root_branch.first);
  else
  // The best root position is somewhere on a branch (a new Node is created)
  {
    Node* new_root = new Node();
    new_root->setId( TreeTools::getMPNUId(tree, tree.getRootId()) );

    double root_branch_length = best_root_branch.first->getDistanceToFather();
    Node* best_root_father = best_root_branch.first->getFather();

    best_root_father->removeSon(best_root_branch.first);
    best_root_father->addSon(new_root);
    new_root->addSon(best_root_branch.first);

    new_root->setDistanceToFather(max(pos * root_branch_length, 1e-6));
    best_root_branch.first->setDistanceToFather(max((1 - pos) * root_branch_length, 1e-6));

    // The two branches leaving the root must have the same branch properties
    const vector<string> branch_properties = best_root_branch.first->getBranchPropertyNames();
    for (vector<string>::const_iterator p = branch_properties.begin(); p != branch_properties.end(); ++p)
    {
      new_root->setBranchProperty(*p, *best_root_branch.first->getBranchProperty(*p));
    }

    tree.rootAt(new_root);
  }

  if (forceBranchRoot) // if we want the root to be on a branch, not on a node
  {
    Node* orig_root = tree.getRootNode();
    vector<Node*> root_sons = orig_root->getSons();
    if (root_sons.size() > 2)
    {
      Node* nearest = root_sons.at(0);
      for (vector<Node*>::iterator n = root_sons.begin(); n !=
           root_sons.end(); ++n)
      {
        if ((**n).getDistanceToFather() < nearest->getDistanceToFather())
          nearest = *n;
      }
      const double d = nearest->getDistanceToFather();
      Node* new_root = new Node();
      new_root->setId( TreeTools::getMPNUId(tree, tree.getRootId()) );
      orig_root->removeSon(nearest);
      orig_root->addSon(new_root);
      new_root->addSon(nearest);
      new_root->setDistanceToFather(d / 2.);
      nearest->setDistanceToFather(d / 2.);
      const vector<string> branch_properties = nearest->getBranchPropertyNames();
      for (vector<string>::const_iterator p = branch_properties.begin(); p != branch_properties.end(); ++p)
      {
        new_root->setBranchProperty(*p, *nearest->getBranchProperty(*p));
      }
      tree.rootAt(new_root);
    }
  }
}

/******************************************************************************/

double TreeTemplateTools::getRadius(TreeTemplate<Node>& tree)
{
  TreeTemplateTools::midRoot(tree, MIDROOT_SUM_OF_SQUARES, false);
  Moments_ moments = getSubtreeMoments_(tree.getRootNode());
  double radius = moments.sum / moments.numberOfLeaves;
  return radius;
}

/******************************************************************************/

void TreeTemplateTools::unresolveUncertainNodes(Node& subtree, double threshold, const std::string& property)
{
  for (size_t i = 0; i < subtree.getNumberOfSons(); ++i)
  {
    Node* son = subtree.getSon(i);
    if (son->getNumberOfSons() > 0)
    {
      // Recursion:
      unresolveUncertainNodes(*son, threshold, property);
      // Deal with this node:
      if (son->hasBranchProperty(property))
      {
        double value = dynamic_cast<Number<double>*>(son->getBranchProperty(property))->getValue();
        if (value < threshold)
        {
          // We remove this branch:
          double brlen = son->getDistanceToFather();
          for (size_t j = 0; j < son->getNumberOfSons(); ++j)
          {
            Node* grandSon = son->getSon(j);
            grandSon->setDistanceToFather(grandSon->getDistanceToFather() + brlen);
            subtree.addSon(i, grandSon);
          }
          subtree.removeSon(son);
          delete son;
        }
      }
    }
  }
}

void TreeTemplateTools::getBestRootInSubtree_(TreeTemplate<Node>& tree, short criterion, Node* node, pair<Node*, map<string, double> >& bestRoot)
{
  const vector<Node*> sons = node->getSons(); // copy
  tree.rootAt(node);

  // Try to place the root on each branch downward node
  for (vector<Node*>::const_iterator son = sons.begin(); son != sons.end(); ++son)
  {
    // Compute the moment of the subtree on son's side
    Moments_ son_moment = getSubtreeMoments_(*son);

    // Compute the moment of the subtree on node's side
    tree.rootAt(*son);
    Moments_ node_moment = getSubtreeMoments_(node);
    tree.rootAt(node);

    /*
     * Get the position of the root on this branch that
     * minimizes the root-to-leaves distances variance.
     *
     * This variance can be written in the form A x^2 + B x + C
     */
    double min_criterion_value;
    double best_position; // 0 is toward the root, 1 is away from it

    const TreeTemplateTools::Moments_& m1 = node_moment;
    const TreeTemplateTools::Moments_& m2 = son_moment;
    const double d = (**son).getDistanceToFather();
    const double n1 = m1.numberOfLeaves;
    const double n2 = m2.numberOfLeaves;

    double A = 0, B = 0, C = 0;
    if (criterion == MIDROOT_SUM_OF_SQUARES)
    {
      A = (n1 + n2) * d * d;
      B = 2 * d * (m1.sum - m2.sum) - 2 * n2 * d * d;
      C = m1.squaresSum + m2.squaresSum
          + 2 * m2.sum * d
          + n2 * d * d;
    }
    else if (criterion == MIDROOT_VARIANCE)
    {
      A = 4 * n1 * n2 * d * d;
      B = 4 * d * ( n2 * m1.sum - n1 * m2.sum - d * n1 * n2);
      C = (n1 + n2) * (m1.squaresSum + m2.squaresSum) + n1 * d * n2 * d
          + 2 * n1 * d * m2.sum - 2 * n2 * d * m1.sum
          - (m1.sum + m2.sum) * (m1.sum + m2.sum);
    }

    if (A < 1e-20)
    {
      min_criterion_value = numeric_limits<double>::max();
      best_position = 0.5;
    }
    else
    {
      min_criterion_value = C - B * B / (4 * A);
      best_position = -B / (2 * A);
      if (best_position < 0)
      {
        best_position = 0;
        min_criterion_value = C;
      }
      else if (best_position > 1)
      {
        best_position = 1;
        min_criterion_value = A + B + C;
      }
    }

    // Is this branch is the best seen, update 'bestRoot'
    if (min_criterion_value < bestRoot.second["score"])
    {
      bestRoot.first = *son;
      bestRoot.second["position"] = best_position;
      bestRoot.second["score"] = min_criterion_value;
    }

    // Recurse
    TreeTemplateTools::getBestRootInSubtree_(tree, criterion, *son, bestRoot);
  }
}

/******************************************************************************/

TreeTemplateTools::Moments_ TreeTemplateTools::getSubtreeMoments_(const Node* node)
{
  TreeTemplateTools::Moments_ moments = {0, 0, 0};

  if (node->isLeaf())
  {
    moments.numberOfLeaves = 1;
  }
  else
  {
    const size_t nsons = node->getNumberOfSons();
    for (size_t i = 0; i < nsons; ++i)
    {
      const Node* son = node->getSon(i);
      const TreeTemplateTools::Moments_ son_moments = TreeTemplateTools::getSubtreeMoments_(son);
      const double d = son->getDistanceToFather();
      moments.numberOfLeaves += son_moments.numberOfLeaves;
      moments.sum += son_moments.sum + d * son_moments.numberOfLeaves;
      moments.squaresSum += son_moments.squaresSum + 2 * d * son_moments.sum + son_moments.numberOfLeaves * d * d;
    }
  }

  return moments;
}

/******************************************************************************/
