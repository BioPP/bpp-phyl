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

using namespace std;

/******************************************************************************/

bool TreeTemplateTools::isMultifurcating(const Node& node)
{
  if (node.getNumberOfSons() > 2)
    return true;
  else
  {
    bool b = false;
    for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
    {
      b = b || isMultifurcating(*node.getSon(i));
    }
    return b;
  }
}

/******************************************************************************/

unsigned int TreeTemplateTools::getNumberOfLeaves(const Node& node)
{
  unsigned int nbLeaves = 0;
  if (node.isLeaf())
  {
    nbLeaves++;
  }
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    nbLeaves += getNumberOfLeaves(*node[i]);
  }
  return nbLeaves;
}

/******************************************************************************/

unsigned int TreeTemplateTools::getNumberOfNodes(const Node& node)
{
  unsigned int nbNodes = 1;
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    nbNodes += getNumberOfNodes(*node[i]);
  }
  return nbNodes;
}

/******************************************************************************/

vector<string> TreeTemplateTools::getLeavesNames(const Node& node)
{
  vector<string> names;
  if (node.isLeaf())
  {
    names.push_back(node.getName());
  }
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    vector<string> subNames = getLeavesNames(*node.getSon(i));
    for (unsigned int j = 0; j < subNames.size(); j++)
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
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    unsigned int c = getDepth(*node[i]) + 1;
    if (c > d)
      d = c;
  }
  return d;
}

/******************************************************************************/

double TreeTemplateTools::getHeight(const Node& node)
{
  double d = 0;
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
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
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
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

  unsigned int colonIndex;
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
      element.length = elt.substr(colonIndex + 1);
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
    }
    else
    {
      // This is a node:
      if (lastP < firstP)
        throw IOException("TreeTemplateTools::getElement(). Invalid format: bad closing parenthesis in " + elt2);
      element.content = elt2.substr(firstP + 1, lastP - firstP - 1);
      string bootstrap = elt2.substr(lastP + 1);
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


Node* TreeTemplateTools::parenthesisToNode(const string& description, bool bootstrap, const string& propertyName, bool withId)
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

  if (elements.size() == 1)
  {
    // This is a leaf:
    // cout << "NODE: LEAF: " << elements[0] << endl;
    string name = TextTools::removeSurroundingWhiteSpaces(elements[0]);
    if (withId)
    {
      StringTokenizer st(name, "_", true, true);
      ostringstream realName;
      for (int i = 0; i < static_cast<int>(st.numberOfRemainingTokens()) - 1; i++)
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
    for (unsigned int i = 0; i < elements.size(); i++)
    {
      // cout << "NODE: SUBNODE: " << i << ", " << elements[i] << endl;
      Node* son = parenthesisToNode(elements[i], bootstrap, propertyName, withId);
      node->addSon(son);
    }
  }
  return node;
}

/******************************************************************************/

TreeTemplate<Node>* TreeTemplateTools::parenthesisToTree(const string& description, bool bootstrap, const string& propertyName, bool withId) throw (Exception)
{
  string::size_type lastP  = description.rfind(')');
  if (lastP == string::npos)
    throw Exception("TreeTemplateTools::parenthesisToTree(). Bad format: no closing parenthesis found.");
  string::size_type firstP = description.find('(');
  if (firstP == string::npos)
    throw Exception("TreeTemplateTools::parenthesisToTree(). Bad format: no opening parenthesis found.");
  string::size_type semi = description.rfind(';');
  if (semi == string::npos)
    throw Exception("TreeTemplateTools::parenthesisToTree(). Bad format: no semi-colon found.");
  if (lastP <= firstP)
    throw Exception("TreeTemplateTools::parenthesisToTree(). Bad format: closing parenthesis before opening parenthesis.");
  string content = description.substr(firstP + 1, lastP - firstP - 1);
  string element = (semi == string::npos) ? description.substr(lastP + 1) : description.substr(lastP + 1, semi - lastP - 1);
  // cout << "TREE: " << content << endl;
  // New root node:
  Node* node = new Node();

  NestedStringTokenizer nt(content, "(", ")", ",");
  vector<string> elements;
  while (nt.hasMoreToken())
  {
    elements.push_back(nt.nextToken());
  }

  if (elements.size() == 1)
  {
    // This is a leaf:
    if (withId)
    {
      StringTokenizer st(elements[0], "_", true, true);
      ostringstream realName;
      for (int i = 0; i < (int)st.numberOfRemainingTokens() - 1; i++)
      {
        if (i != 0)
        {
          realName << "_";
        }
        realName << st.getToken(i);
      }
      node->setName(realName.str());
      node->setName(realName.str());
      node->setId(TextTools::toInt(st.getToken(1)));
    }
    else
    {
      node->setName(elements[0]);
    }
  }
  else
  {
    //cout << element << endl;
    // This is a node:
    for (unsigned int i = 0; i < elements.size(); i++)
    {
      Node* son = parenthesisToNode(elements[i], bootstrap, propertyName, withId);
      node->addSon(son);
    }
    if (!TextTools::isEmpty(element))
    {
      StringTokenizer st(element, ":");
      string lengthS = "";
      string bootstrapS = "";
      if (st.numberOfRemainingTokens() == 1)
      {
        if (element[0] == ':')
          lengthS = st.nextToken();
        else 
          bootstrapS = st.nextToken();
      }
      else
      {
        bootstrapS = st.nextToken();
        lengthS = st.nextToken();
      }
      //cout << "L=" << lengthS << "\tB=" << bootstrapS << endl;
      if (!TextTools::isEmpty(lengthS))
      {
        node->setDistanceToFather(TextTools::toDouble(lengthS));
        //cout << "NODE: LENGTH: " << lengthS << endl;
      }
      if (!TextTools::isEmpty(bootstrapS))
      {
        if (withId)
        {
          node->setId(TextTools::toInt(bootstrapS));
        }
        else
        {
          if (bootstrap)
          {
            node->setBranchProperty(TreeTools::BOOTSTRAP, Number<double>(TextTools::toDouble(bootstrapS)));
            // cout << "NODE: BOOTSTRAP: " << * elt.bootstrap << endl;
          }
          else
          {
            node->setBranchProperty(propertyName, BppString(bootstrapS));
          }
        }
      }
    }
  }
  TreeTemplate<Node>* tree = new TreeTemplate<Node>();
  tree->setRootNode(node);
  if (!withId)
  {
    tree->resetNodesId();
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
    for (unsigned int i = 1; i < node.getNumberOfSons(); i++)
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
    for (unsigned int i = 1; i < node.getNumberOfSons(); i++)
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
  if (node->isLeaf())
  {
    s << node->getName();
    for (unsigned int i = 0; i < node->getNumberOfSons(); ++i)
    {
      s << "," << nodeToParenthesis(*node->getSon(i), writeId);
    }
  }
  else
  {
    s << nodeToParenthesis(*node->getSon(0), writeId);
    for (unsigned int i = 1; i < node->getNumberOfSons(); ++i)
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
    for (unsigned int i = 0; i < node->getNumberOfSons(); i++)
    {
      s << "," << nodeToParenthesis(*node->getSon(i), bootstrap, propertyName);
    }
  }
  else
  {
    s << nodeToParenthesis(*node->getSon(0), bootstrap, propertyName);
    for (unsigned int i = 1; i < node->getNumberOfSons(); i++)
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
    if (node->hasBranchProperty(propertyName)) {
      const BppString* ppt =dynamic_cast<const BppString*>(node->getBranchProperty(propertyName));
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
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    Vdouble sonBrLen = getBranchLengths(*node.getSon(i));
    for (unsigned int j = 0; j < sonBrLen.size(); j++)
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
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    length += getTotalLength(*node.getSon(i), true);
  }
  return length;
}

/******************************************************************************/

void TreeTemplateTools::setBranchLengths(Node& node, double brLen)
{
  node.setDistanceToFather(brLen);
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    setBranchLengths(*node.getSon(i), brLen);
  }
}

/******************************************************************************/

void TreeTemplateTools::deleteBranchLengths(Node& node)
{
  node.deleteDistanceToFather();
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    deleteBranchLengths(*node.getSon(i));
  }
}

/******************************************************************************/

void TreeTemplateTools::setVoidBranchLengths(Node& node, double brLen)
{
  if (!node.hasDistanceToFather())
    node.setDistanceToFather(brLen);
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
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
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    scaleTree(*node.getSon(i), factor);
  }
}

/******************************************************************************/

TreeTemplate<Node>* TreeTemplateTools::getRandomTree(vector<string>& leavesNames, bool rooted)
{
  if (leavesNames.size() == 0)
    return 0;  // No taxa.
  // This vector will contain all nodes.
  // Start with all leaves, and then group nodes randomly 2 by 2.
  // Att the end, contains only the root node of the tree.
  vector<Node*> nodes(leavesNames.size());
  // Create all leaves nodes:
  for (unsigned int i = 0; i < leavesNames.size(); ++i)
  {
    nodes[i] = new Node(leavesNames[i]);
  }
  // Now group all nodes:
  while (nodes.size() > (rooted ? 2 : 3))
  {
    // Select random nodes:
    int pos1 = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(nodes.size());
    Node* node1 = nodes[pos1];
    nodes.erase(nodes.begin() + pos1);
    int pos2 = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(nodes.size());
    Node* node2 = nodes[pos2];
    nodes.erase(nodes.begin() + pos2);
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

vector<Node*> TreeTemplateTools::getPathBetweenAnyTwoNodes(Node& node1, Node& node2, bool includeAncestor)
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
  // Must check that the two nodes have the same root!!!

  int tmp1 = pathMatrix1.size() - 1;
  int tmp2 = pathMatrix2.size() - 1;

  while ((tmp1 >= 0) && (tmp2 >= 0))
  {
    if (pathMatrix1[tmp1] != pathMatrix2[tmp2])
      break;
    tmp1--; tmp2--;
  }

  for (int y = 0; y <= tmp1; ++y)
  {
    path.push_back(pathMatrix1[y]);
  }
  if (includeAncestor)
    path.push_back(pathMatrix1[tmp1 + 1]);  // pushing once, the Node that was common to both.
  for (int j = tmp2; j >= 0; --j)
  {
    path.push_back(pathMatrix2[j]);
  }
  return path;
}

/******************************************************************************/

vector<const Node*> TreeTemplateTools::getPathBetweenAnyTwoNodes(const Node& node1, const Node& node2, bool includeAncestor)
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
  // Must check that the two nodes have the same root!!!

  int tmp1 = pathMatrix1.size() - 1;
  int tmp2 = pathMatrix2.size() - 1;

  while ((tmp1 >= 0) && (tmp2 >= 0))
  {
    if (pathMatrix1[tmp1] != pathMatrix2[tmp2])
      break;
    tmp1--; tmp2--;
  }

  for (int y = 0; y <= tmp1; ++y)
  {
    path.push_back(pathMatrix1[y]);
  }
  if (includeAncestor)
    path.push_back(pathMatrix1[tmp1 + 1]);  // pushing once, the Node that was common to both.
  for (int j = tmp2; j >= 0; --j)
  {
    path.push_back(pathMatrix2[j]);
  }
  return path;
}

/******************************************************************************/

double TreeTemplateTools::getDistanceBetweenAnyTwoNodes(const Node& node1, const Node& node2)
{
  vector<const Node*> path = getPathBetweenAnyTwoNodes(node1, node2, false);
  double d = 0;
  for (unsigned int i = 0; i < path.size(); i++)
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
  for (unsigned int i = 0; i < node->getNumberOfSons(); ++i)
  {
    const Node* son = node->getSon(i);
    processDistsInSubtree_(son, matrix, leavesDists[son]); // recursivity
  }
  // Write leaf-leaf distances to the distance matrix.
  // Only pairs in which the two leaves belong to different
  // sons are considered.
  for (unsigned int son1_loc = 0; son1_loc < node->getNumberOfSons(); ++son1_loc)
  {
    for (unsigned int son2_loc = 0; son2_loc < son1_loc; ++son2_loc)
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
  for (unsigned int k = 0; k < neighbors.size(); k++)
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
  for (unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    incrementAllIds(node->getSon(i), increment);
  }
}

/******************************************************************************/

void TreeTemplateTools::getNodePropertyNames(const Node& node, vector<string>& propertyNames)
{
  VectorTools::extend(propertyNames, node.getNodePropertyNames());
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    getNodePropertyNames(*node.getSon(i), propertyNames);
  }
}

void TreeTemplateTools::getNodeProperties(const Node& node, const string& propertyName, map<int, const Clonable*>& properties)
{
  if (node.hasNodeProperty(propertyName))
    properties[node.getId()] = node.getNodeProperty(propertyName);
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    getNodeProperties(*node.getSon(i), propertyName, properties);
  }
}

void TreeTemplateTools::getNodeProperties(Node& node, const string& propertyName, map<int, Clonable*>& properties)
{
  if (node.hasNodeProperty(propertyName))
    properties[node.getId()] = node.getNodeProperty(propertyName);
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    getNodeProperties(*node.getSon(i), propertyName, properties);
  }
}

/******************************************************************************/

void TreeTemplateTools::getBranchPropertyNames(const Node& node, vector<string>& propertyNames)
{
  VectorTools::extend(propertyNames, node.getBranchPropertyNames());
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    getBranchPropertyNames(*node.getSon(i), propertyNames);
  }
}

void TreeTemplateTools::getBranchProperties(const Node& node, const string& propertyName, map<int, const Clonable*>& properties)
{
  if (node.hasBranchProperty(propertyName))
    properties[node.getId()] = node.getBranchProperty(propertyName);
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    getBranchProperties(*node.getSon(i), propertyName, properties);
  }
}

void TreeTemplateTools::getBranchProperties(Node& node, const string& propertyName, map<int, Clonable*>& properties)
{
  if (node.hasBranchProperty(propertyName))
    properties[node.getId()] = node.getBranchProperty(propertyName);
  for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
  {
    getBranchProperties(*node.getSon(i), propertyName, properties);
  }
}

/******************************************************************************/

bool TreeTemplateTools::haveSameOrderedTopology(const Node& n1, const Node& n2)
{
  if (n1.isLeaf() && n2.isLeaf() && n1.getName() != n2.getName())
    return false;
  unsigned int nl1 = n1.getNumberOfSons();
  unsigned int nl2 = n2.getNumberOfSons();
  if (nl1 != nl2)
    return false;

  bool test = true;
  for (unsigned int i = 0; test && i < n1.getNumberOfSons(); ++i)
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
    vector<unsigned int> nbSons;
    vector<string> firstLeaves;
    for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
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

