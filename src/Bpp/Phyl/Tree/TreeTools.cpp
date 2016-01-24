//
// File: TreeTools.cpp
// Created by: Julien Dutheil
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

#include "TreeTools.h"
#include "Tree.h"
#include "BipartitionTools.h"
#include "../Model/Nucleotide/JCnuc.h"
#include "../Distance/DistanceEstimation.h"
#include "../Distance/BioNJ.h"
#include "../Parsimony/DRTreeParsimonyScore.h"
#include "../OptimizationTools.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/BppString.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <sstream>

using namespace std;

/******************************************************************************/

const string TreeTools::BOOTSTRAP = "bootstrap";

/******************************************************************************/

vector<int> TreeTools::getLeavesId(const Tree& tree, int nodeId) throw (NodeNotFoundException)
{
  vector<int> leaves;
  getLeavesId(tree, nodeId, leaves);
  return leaves;
}

void TreeTools::getLeavesId(const Tree& tree, int nodeId, std::vector<int>& leaves) throw (NodeNotFoundException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::getLeavesId", nodeId);
  if (tree.isLeaf(nodeId))
  {
    leaves.push_back(nodeId);
  }
  vector<int> sons = tree.getSonsId(nodeId);
  for (size_t i = 0; i < sons.size(); i++)
  {
    getLeavesId(tree, sons[i], leaves);
  }
}

size_t TreeTools::getNumberOfLeaves(const Tree& tree, int nodeId) throw (NodeNotFoundException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::getNumberOfLeaves", nodeId);

  size_t n = 0;
  if (tree.isLeaf(nodeId))
  {
    ++n;
  }
  vector<int> sons = tree.getSonsId(nodeId);
  for (size_t i = 0; i < sons.size(); ++i)
  {
    n += getNumberOfLeaves(tree, sons[i]);
  }
  return n;
}

/******************************************************************************/

int TreeTools::getLeafId(const Tree& tree, int nodeId, const std::string& name)
throw (NodeNotFoundException)
{
  int* id = NULL;
  searchLeaf(tree, nodeId, name, id);
  if (id == NULL)
    throw NodeNotFoundException("TreeTools::getLeafId().", name);
  else
  {
    int i = *id;
    delete id;
    return i;
  }
}

void TreeTools::searchLeaf(const Tree& tree, int nodeId, const string& name, int*& id)
throw (NodeNotFoundException)
{
  if (tree.isLeaf(nodeId))
  {
    if (tree.getNodeName(nodeId) == name)
    {
      id = new int(nodeId);
      return;
    }
  }
  vector<int> sons;
  for (size_t i = 0; i < sons.size(); i++)
  {
    searchLeaf(tree, nodeId, name, id);
  }
}

/******************************************************************************/

vector<int> TreeTools::getNodesId(const Tree& tree, int nodeId) throw (NodeNotFoundException)
{
  vector<int> nodes;
  getNodesId(tree, nodeId, nodes);
  return nodes;
}

void TreeTools::getNodesId(const Tree& tree, int nodeId, vector<int>& nodes) throw (NodeNotFoundException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::getNodesId", nodeId);
  vector<int> sons = tree.getSonsId(nodeId);
  for (size_t i = 0; i < sons.size(); i++)
  {
    getNodesId(tree, sons[i], nodes);
  }
  nodes.push_back(nodeId);
}

/******************************************************************************/

size_t TreeTools::getDepth(const Tree& tree, int nodeId) throw (NodeNotFoundException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::getDepth", nodeId);
  size_t d = 0;
  vector<int> sons = tree.getSonsId(nodeId);
  for (size_t i = 0; i < sons.size(); i++)
  {
    size_t c = getDepth(tree, sons[i]) + 1;
    if (c > d)
      d = c;
  }
  return d;
}

/******************************************************************************/

size_t TreeTools::getDepths(const Tree& tree, int nodeId, map<int, size_t>& depths) throw (NodeNotFoundException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::getDepth", nodeId);
  size_t d = 0;
  vector<int> sons = tree.getSonsId(nodeId);
  for (size_t i = 0; i < sons.size(); i++)
  {
    size_t c = getDepths(tree, sons[i], depths) + 1;
    if (c > d)
      d = c;
  }
  depths[nodeId] = d;
  return d;
}

/******************************************************************************/

double TreeTools::getHeight(const Tree& tree, int nodeId) throw (NodeNotFoundException, NodeException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::getHeight", nodeId);
  double d = 0;
  vector<int> sons = tree.getSonsId(nodeId);
  for (size_t i = 0; i < sons.size(); i++)
  {
    double dist = 0;
    if (tree.hasDistanceToFather(sons[i]))
      dist = tree.getDistanceToFather(sons[i]);
    else
      throw NodeException("Node without branch length.", sons[i]);
    double c = getHeight(tree, sons[i]) + dist;
    if (c > d)
      d = c;
  }
  return d;
}

/******************************************************************************/

double TreeTools::getHeights(const Tree& tree, int nodeId, map<int, double>& heights) throw (NodeNotFoundException, NodeException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::getHeight", nodeId);
  double d = 0;
  vector<int> sons = tree.getSonsId(nodeId);
  for (size_t i = 0; i < sons.size(); i++)
  {
    double dist = 0;
    if (tree.hasDistanceToFather(sons[i]))
      dist = tree.getDistanceToFather(sons[i]);
    else
      throw NodeException("Node without branch length.", sons[i]);
    double c = getHeights(tree, sons[i], heights) + dist;
    if (c > d)
      d = c;
  }
  heights[nodeId] = d;
  return d;
}

/******************************************************************************/

string TreeTools::nodeToParenthesis(const Tree& tree, int nodeId, bool writeId) throw (NodeNotFoundException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::nodeToParenthesis", nodeId);
  ostringstream s;
  if (tree.isLeaf(nodeId))
  {
    s << tree.getNodeName(nodeId);
  }
  else
  {
    s << "(";
    vector<int> sonsId = tree.getSonsId(nodeId);
    s << nodeToParenthesis(tree, sonsId[0], writeId);
    for (size_t i = 1; i < sonsId.size(); i++)
    {
      s << "," << nodeToParenthesis(tree, sonsId[i], writeId);
    }
    s << ")";
  }
  if (writeId)
  {
    if (tree.isLeaf(nodeId))
      s << "_";
    s << nodeId;
  }
  else
  {
    if (tree.hasBranchProperty(nodeId, BOOTSTRAP))
      s << (dynamic_cast<const Number<double>*>(tree.getBranchProperty(nodeId, BOOTSTRAP))->getValue());
  }
  if (tree.hasDistanceToFather(nodeId))
    s << ":" << tree.getDistanceToFather(nodeId);
  return s.str();
}

/******************************************************************************/

string TreeTools::nodeToParenthesis(const Tree& tree, int nodeId, bool bootstrap, const string& propertyName) throw (NodeNotFoundException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::nodeToParenthesis", nodeId);
  ostringstream s;

  if (tree.isLeaf(nodeId))
  {
    s << tree.getNodeName(nodeId);
  }
  else
  {
    s << "(";
    vector<int> sonsId = tree.getSonsId(nodeId);
    s << nodeToParenthesis(tree, sonsId[0], bootstrap, propertyName);
    for (size_t i = 1; i < sonsId.size(); i++)
    {
      s << "," << nodeToParenthesis(tree, sonsId[i], bootstrap, propertyName);
    }
    s << ")";

    if (bootstrap)
    {
      if (tree.hasBranchProperty(nodeId, BOOTSTRAP))
        s << (dynamic_cast<const Number<double>*>(tree.getBranchProperty(nodeId, BOOTSTRAP))->getValue());
    }
    else
    {
      if (tree.hasBranchProperty(nodeId, propertyName))
        s << dynamic_cast<const BppString*>(tree.getBranchProperty(nodeId, propertyName))->toSTL();
    }
  }
  if (tree.hasDistanceToFather(nodeId))
    s << ":" << tree.getDistanceToFather(nodeId);
  
  return s.str();
}

/******************************************************************************/

string TreeTools::treeToParenthesis(const Tree& tree, bool writeId)
{
  ostringstream s;
  s << "(";
  int rootId = tree.getRootId();
  vector<int> sonsId = tree.getSonsId(rootId);
  if (tree.isLeaf(rootId))
  {
    s << tree.getNodeName(rootId);
    for (size_t i = 0; i < sonsId.size(); i++)
    {
      s << "," << nodeToParenthesis(tree, sonsId[i], writeId);
    }
  }
  else
  {
    if (sonsId.size() > 0)
    {
      s << nodeToParenthesis(tree, sonsId[0], writeId);
      for (size_t i = 1; i < sonsId.size(); i++)
      {
        s << "," << nodeToParenthesis(tree, sonsId[i], writeId);
      }
    }
    // Otherwise, this is an empty tree!
  }
  s << ");" << endl;
  return s.str();
}

/******************************************************************************/

string TreeTools::treeToParenthesis(const Tree& tree, bool bootstrap, const string& propertyName)
{
  ostringstream s;
  s << "(";
  int rootId = tree.getRootId();
  vector<int> sonsId = tree.getSonsId(rootId);
  
  if (tree.isLeaf(rootId))
  {
    s << tree.getNodeName(rootId);
    for (size_t i = 0; i < sonsId.size(); i++)
    {
      s << "," << nodeToParenthesis(tree, sonsId[i], bootstrap, propertyName);
    }
  }
  else
  {
    s << nodeToParenthesis(tree, sonsId[0], bootstrap, propertyName);
    for (size_t i = 1; i < sonsId.size(); i++)
    {
      s << "," << nodeToParenthesis(tree, sonsId[i], bootstrap, propertyName);
    }
  }
  s << ")";
  if (bootstrap)
  {
    if (tree.hasBranchProperty(rootId, BOOTSTRAP))
      s << (dynamic_cast<const Number<double>*>(tree.getBranchProperty(rootId, BOOTSTRAP))->getValue());
  }
  else
  {
    if (tree.hasBranchProperty(rootId, propertyName))
      s << dynamic_cast<const BppString*>(tree.getBranchProperty(rootId, propertyName))->toSTL();
  }
  s << ";" << endl;
  
  return s.str();
}

/******************************************************************************/

vector<int> TreeTools::getPathBetweenAnyTwoNodes(const Tree& tree, int nodeId1, int nodeId2, bool includeAncestor)
throw (NodeNotFoundException)
{
  if (!tree.hasNode(nodeId1))
    throw NodeNotFoundException("TreeTools::getPathBetweenAnyTwoNodes", nodeId1);
  if (!tree.hasNode(nodeId2))
    throw NodeNotFoundException("TreeTools::getPathBetweenAnyTwoNodes", nodeId2);
  vector<int> path;
  vector<int> pathMatrix1;
  vector<int> pathMatrix2;

  int nodeUp = nodeId1;
  while (tree.hasFather(nodeUp))
  {
    pathMatrix1.push_back(nodeUp);
    nodeUp = tree.getFatherId(nodeUp);
  }
  pathMatrix1.push_back(nodeUp); // The root.

  nodeUp = nodeId2;
  while (tree.hasFather(nodeUp))
  {
    pathMatrix2.push_back(nodeUp);
    nodeUp = tree.getFatherId(nodeUp);
  }
  pathMatrix2.push_back(nodeUp); // The root.
  // Must check that the two nodes have the same root!!!

  size_t tmp1 = pathMatrix1.size();
  size_t tmp2 = pathMatrix2.size();

  while ((tmp1 > 0) && (tmp2 > 0))
  {
    if (pathMatrix1[tmp1 - 1] != pathMatrix2[tmp2 - 1])
      break;
    tmp1--; tmp2--;
  }
  // (tmp1 - 1) and (tmp2 - 1) now point toward the first non-common nodes

  for (size_t y = 0; y < tmp1; ++y)
  {
    path.push_back(pathMatrix1[y]);
  }
  if (includeAncestor)
    path.push_back(pathMatrix1[tmp1]);  // pushing once, the Node that was common to both.
  for (size_t j = tmp2; j > 0; --j)
  {
    path.push_back(pathMatrix2[j - 1]);
  }
  return path;
}

/******************************************************************************/

double TreeTools::getDistanceBetweenAnyTwoNodes(const Tree& tree, int nodeId1, int nodeId2)
{
  if (!tree.hasNode(nodeId1))
    throw NodeNotFoundException("TreeTools::getDistanceBetweenAnyTwoNodes", nodeId1);
  if (!tree.hasNode(nodeId2))
    throw NodeNotFoundException("TreeTools::getDistanceBetweenAnyTwoNodes", nodeId2);
  vector<int> path = getPathBetweenAnyTwoNodes(tree, nodeId1, nodeId2, false);
  double d = 0;
  for (size_t i = 0; i < path.size(); i++)
  {
    d += tree.getDistanceToFather(path[i]);
  }
  return d;
}

/******************************************************************************/

Vdouble TreeTools::getBranchLengths(const Tree& tree, int nodeId) throw (NodeNotFoundException, NodeException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::getBranchLengths", nodeId);
  Vdouble brLen(1);
  if (tree.hasDistanceToFather(nodeId))
    brLen[0] = tree.getDistanceToFather(nodeId);
  else
    throw NodeException("TreeTools::getbranchLengths(). No branch length.", nodeId);
  vector<int> sons = tree.getSonsId(nodeId);
  for (size_t i = 0; i < sons.size(); i++)
  {
    Vdouble sonBrLen = getBranchLengths(tree, sons[i]);
    for (size_t j = 0; j < sonBrLen.size(); j++)
    {
      brLen.push_back(sonBrLen[j]);
    }
  }
  return brLen;
}

/******************************************************************************/

double TreeTools::getTotalLength(const Tree& tree, int nodeId, bool includeAncestor) throw (NodeNotFoundException, NodeException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::getTotalLength", nodeId);
  if (includeAncestor && !tree.hasDistanceToFather(nodeId))
    throw NodeException("TreeTools::getTotalLength(). No branch length.", nodeId);
  double length = includeAncestor ? tree.getDistanceToFather(nodeId) : 0;
  vector<int> sons = tree.getSonsId(nodeId);
  for (size_t i = 0; i < sons.size(); i++)
  {
    length += getTotalLength(tree, sons[i], true);
  }
  return length;
}

/******************************************************************************/

void TreeTools::setBranchLengths(Tree& tree, int nodeId, double brLen) throw (NodeNotFoundException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::setBranchLengths", nodeId);
  vector<int> nodes = getNodesId(tree, nodeId);
  for (size_t i = 0; i < nodes.size(); i++)
  {
    tree.setDistanceToFather(nodes[i], brLen);
  }
}

/******************************************************************************/

void TreeTools::setVoidBranchLengths(Tree& tree, int nodeId, double brLen) throw (NodeNotFoundException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::setVoidBranchLengths", nodeId);
  vector<int> nodes = getNodesId(tree, nodeId);
  for (size_t i = 0; i < nodes.size(); i++)
  {
    if (!tree.hasDistanceToFather(nodes[i]))
      tree.setDistanceToFather(nodes[i], brLen);
  }
}

/******************************************************************************/

void TreeTools::scaleTree(Tree& tree, int nodeId, double factor) throw (NodeNotFoundException, NodeException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::scaleTree", nodeId);
  vector<int> nodes = getNodesId(tree, nodeId);
  for (size_t i = 0; i < nodes.size(); i++)
  {
    if (tree.hasFather(nodes[i]))
    {
      if (!tree.hasDistanceToFather(nodes[i]))
        throw NodeException("TreeTools::scaleTree(). Branch with no length", nodes[i]);
      double brLen = tree.getDistanceToFather(nodes[i]) * factor;
      tree.setDistanceToFather(nodes[i], brLen);
    }
  }
}

/******************************************************************************/

size_t TreeTools::initBranchLengthsGrafen(Tree& tree, int nodeId) throw (NodeNotFoundException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::initBranchLengthsGrafen", nodeId);
  vector<int> sons = tree.getSonsId(nodeId);
  vector<size_t> h(sons.size());
  for (size_t i = 0; i < sons.size(); i++)
  {
    h[i] = initBranchLengthsGrafen(tree, sons[i]);
  }
  size_t thish = sons.size() == 0 ? 0 : VectorTools::sum<size_t>(h) + sons.size() - 1;
  for (size_t i = 0; i < sons.size(); i++)
  {
    tree.setDistanceToFather(sons[i], (double)(thish - h[i]));
  }
  return thish;
}

void TreeTools::initBranchLengthsGrafen(Tree& tree)
{
  initBranchLengthsGrafen(tree, tree.getRootId());
}

/******************************************************************************/

void TreeTools::computeBranchLengthsGrafen(
  Tree& tree,
  int nodeId,
  double power,
  double total,
  double& height,
  double& heightRaised)
throw (NodeNotFoundException, NodeException)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::computeBranchLengthsGrafen", nodeId);
  vector<int> sons = tree.getSonsId(nodeId);
  vector<double> hr(sons.size());
  height = 0;
  for (size_t i = 0; i < sons.size(); i++)
  {
    int son = sons[i];
    if (tree.hasDistanceToFather(son))
    {
      double h;
      computeBranchLengthsGrafen(tree, sons[i], power, total, h, hr[i]);
      double d = h + tree.getDistanceToFather(son);
      if (d > height)
        height = d;
    }
    else
      throw NodeException ("TreeTools::computeBranchLengthsGrafen. Branch length lacking.", son);
  }
  heightRaised = pow(height / total, power) * total;
  for (size_t i = 0; i < sons.size(); i++)
  {
    tree.setDistanceToFather(sons[i], heightRaised - hr[i]);
  }
}

void TreeTools::computeBranchLengthsGrafen(Tree& tree, double power, bool init)
throw (NodeException)
{
  int rootId = tree.getRootId();
  if (init)
  {
    initBranchLengthsGrafen(tree);
  }
  // Scale by total heigth:
  double totalHeight = getHeight(tree, rootId);
  double h, hr;
  computeBranchLengthsGrafen(tree, rootId, power, totalHeight, h, hr);
}

/******************************************************************************/

double TreeTools::convertToClockTree(Tree& tree, int nodeId, bool noneg)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::convertToClockTree", nodeId);
  vector<int> sons = tree.getSonsId(nodeId);
  vector<double> h(sons.size());
  // We compute the mean height:
  double l = 0;
  double maxh = -1.;
  for (size_t i = 0; i < sons.size(); i++)
  {
    int son = sons[i];
    if (tree.hasDistanceToFather(son))
    {
      h[i] = convertToClockTree(tree, son);
      if (h[i] > maxh)
        maxh = h[i];
      l += h[i] + tree.getDistanceToFather(son);
    }
    else
      throw NodeException ("TreeTools::convertToClockTree. Branch length lacking.", son);
  }
  if (sons.size() > 0)
    l /= (double)sons.size();
  if (l < maxh)
    l = maxh;
  for (size_t i = 0; i < sons.size(); i++)
  {
    tree.setDistanceToFather(sons[i], l - h[i]);
  }
  return l;
}

/******************************************************************************/

double TreeTools::convertToClockTree2(Tree& tree, int nodeId)
{
  if (!tree.hasNode(nodeId))
    throw NodeNotFoundException("TreeTools::convertToClockTree2", nodeId);
  vector<int> sons = tree.getSonsId(nodeId);
  vector<double> h(sons.size());
  // We compute the mean height:
  double l = 0;
  double maxh = -1.;
  for (size_t i = 0; i < sons.size(); i++)
  {
    int son = sons[i];
    if (tree.hasDistanceToFather(son))
    {
      h[i] = convertToClockTree2(tree, son);
      if (h[i] > maxh)
        maxh = h[i];
      l += h[i] + tree.getDistanceToFather(son);
    }
    else
      throw NodeException ("TreeTools::convertToClockTree2. Branch length lacking.", son);
  }
  if (sons.size() > 0)
    l /= (double)sons.size();
  for (size_t i = 0; i < sons.size(); i++)
  {
    scaleTree(tree, sons[i], h[i] > 0 ? l / h[i] : 0);
  }
  return l;
}

/******************************************************************************/

DistanceMatrix* TreeTools::getDistanceMatrix(const Tree& tree)
{
  vector<string> names = tree.getLeavesNames();
  DistanceMatrix* mat = new DistanceMatrix(names);
  for (size_t i = 0; i < names.size(); i++)
  {
    (*mat)(i, i) = 0;
    for (size_t j = 0; j < i; j++)
    {
      (*mat)(i, j) = (*mat)(j, i) = getDistanceBetweenAnyTwoNodes(tree, tree.getLeafId(names[i]), tree.getLeafId(names[j]));
    }
  }
  return mat;
}

/******************************************************************************/

void TreeTools::midpointRooting(Tree& tree)
{
  throw Exception("TreeTools::midpointRooting(Tree). This function is deprecated, use TreeTemplateTools::midRoot instead!");
  if (tree.isRooted())
    tree.unroot();
  DistanceMatrix* dist = getDistanceMatrix(tree);
  vector<size_t> pos = MatrixTools::whichMax(dist->asMatrix());
  double dmid = (*dist)(pos[0], pos[1]) / 2;
  int id1 = tree.getLeafId(dist->getName(pos[0]));
  int id2 = tree.getLeafId(dist->getName(pos[1]));
  int rootId = tree.getRootId();
  double d1 = getDistanceBetweenAnyTwoNodes(tree, id1, rootId);
  double d2 = getDistanceBetweenAnyTwoNodes(tree, id2, rootId);
  int current = d2 > d1 ? id2 : id1;
  delete dist;
  double l = tree.getDistanceToFather(current);
  double c = l;
  while (c < dmid)
  {
    current = tree.getFatherId(current);
    l = tree.getDistanceToFather(current);
    c += l;
  }
  tree.newOutGroup(current);
  int brother = tree.getSonsId(tree.getRootId())[1];
  if (brother == current)
    brother = tree.getSonsId(tree.getRootId())[0];
  tree.setDistanceToFather(current, l - (c - dmid));
  tree.setDistanceToFather(brother, c - dmid);
}

/******************************************************************************/

int TreeTools::getMaxId(const Tree& tree, int id)
{
  int maxId = id;
  vector<int> sonsId = tree.getSonsId(id);
  for (size_t i = 0; i < sonsId.size(); i++)
  {
    int subMax = getMaxId(tree, sonsId[i]);
    if (subMax > maxId)
      maxId = subMax;
  }
  return maxId;
}

/******************************************************************************/

int TreeTools::getMPNUId(const Tree& tree, int id)
{
  vector<int> ids = getNodesId(tree, id);
  sort(ids.begin(), ids.end());
  // Check if some id is "missing" in the subtree:
  for (size_t i = 0; i < ids.size(); i++)
  {
    if (ids[i] != (int)i)
      return (int)i;
  }
  // Well, all ids are from 0 to n, so send n+1:
  return (int)ids.size();
}

/******************************************************************************/

bool TreeTools::checkIds(const Tree& tree, bool throwException) throw (Exception)
{
  vector<int> ids = tree.getNodesId();
  sort(ids.begin(), ids.end());
  for (size_t i = 1; i < ids.size(); i++)
  {
    if (ids[i] == ids[i - 1])
    {
      if (throwException)
        throw Exception("TreeTools::checkIds. This id is used at least twice: " + TextTools::toString(ids[i]));
      return false;
    }
  }
  return true;
}

/******************************************************************************/

VectorSiteContainer* TreeTools::MRPEncode(const vector<Tree*>& vecTr)
{
  vector<BipartitionList*> vecBipL;
  for (size_t i = 0; i < vecTr.size(); i++)
  {
    vecBipL.push_back(new BipartitionList(*vecTr[i]));
  }

  VectorSiteContainer* cont = BipartitionTools::MRPEncode(vecBipL);

  for (size_t i = 0; i < vecTr.size(); i++)
  {
    delete vecBipL[i];
  }

  return cont;
}

/******************************************************************************/

VectorSiteContainer* TreeTools::MRPEncodeMultilabel(const vector<Tree*>& vecTr)
{
    vector<BipartitionList*> vecBipL;
    for (size_t i = 0; i < vecTr.size(); i++)
    {
        vecBipL.push_back(new BipartitionList(*vecTr[i]));
    }
    
    VectorSiteContainer* cont = BipartitionTools::MRPEncodeMultilabel(vecBipL);
    
    for (size_t i = 0; i < vecTr.size(); i++)
    {
        delete vecBipL[i];
    }
    
    return cont;
}

/******************************************************************************/

bool TreeTools::haveSameTopology(const Tree& tr1, const Tree& tr2)
{
  size_t jj, nbbip;
  BipartitionList* bipL1, * bipL2;
  vector<size_t> size1, size2;

  /* compare sets of leaves */
  if (!VectorTools::haveSameElements(tr1.getLeavesNames(), tr2.getLeavesNames()))
    return false;

  /* construct bipartitions */
  bipL1 = new BipartitionList(tr1, true);
  bipL1->removeTrivialBipartitions();
  bipL1->removeRedundantBipartitions();
  bipL1->sortByPartitionSize();
  bipL2 = new BipartitionList(tr2, true);
  bipL2->removeTrivialBipartitions();
  bipL2->removeRedundantBipartitions();
  bipL2->sortByPartitionSize();

  /* compare numbers of bipartitions */
  if (bipL1->getNumberOfBipartitions() != bipL2->getNumberOfBipartitions())
    return false;
  nbbip = bipL1->getNumberOfBipartitions();

  /* compare partition sizes */
  for (size_t i = 0; i < nbbip; i++)
  {
    size1.push_back(bipL1->getPartitionSize(i));
    size2.push_back(bipL1->getPartitionSize(i));
    if (size1[i] != size2[i])
      return false;
  }

  /* compare bipartitions */
  for (size_t i = 0; i < nbbip; i++)
  {
    for (jj = 0; jj < nbbip; jj++)
    {
      if (size1[i] == size2[jj] && BipartitionTools::areIdentical(*bipL1, i, *bipL2, jj))
        break;
    }
    if (jj == nbbip)
      return false;
  }

  return true;
}

/******************************************************************************/

int TreeTools::robinsonFouldsDistance(const Tree& tr1, const Tree& tr2, bool checkNames, int* missing_in_tr2, int* missing_in_tr1) throw (Exception)
{
  BipartitionList* bipL1, * bipL2;
  size_t i, j;
  vector<size_t> size1, size2;
  vector<bool> bipOK2;


  if (checkNames && !VectorTools::haveSameElements(tr1.getLeavesNames(), tr2.getLeavesNames()))
    throw Exception("Distinct leaf sets between trees ");

  /* prepare things */
  int missing1 = 0;
  int missing2 = 0;

  bipL1 = new BipartitionList(tr1, true);
  bipL1->removeTrivialBipartitions();
  bipL1->sortByPartitionSize();
  bipL2 = new BipartitionList(tr2, true);
  bipL2->removeTrivialBipartitions();
  bipL2->sortByPartitionSize();


  for (i = 0; i < bipL1->getNumberOfBipartitions(); i++)
  {
    size1.push_back(bipL1->getPartitionSize(i));
  }
  for (i = 0; i < bipL2->getNumberOfBipartitions(); i++)
  {
    size2.push_back(bipL2->getPartitionSize(i));
  }

  for (i = 0; i < bipL2->getNumberOfBipartitions(); i++)
  {
    bipOK2.push_back(false);
  }

  /* main loops */

  for (i = 0; i < bipL1->getNumberOfBipartitions(); i++)
  {
    for (j = 0; j < bipL2->getNumberOfBipartitions(); j++)
    {
      if (bipOK2[j])
        continue;
      if (size1[i] == size2[j] && BipartitionTools::areIdentical(*bipL1, i, *bipL2, j))
      {
        bipOK2[j] = true;
        break;
      }
    }
    if (j == bipL2->getNumberOfBipartitions())
      missing2++;
  }

  missing1 = static_cast<int>(bipL2->getNumberOfBipartitions()) - static_cast<int>(bipL1->getNumberOfBipartitions()) + missing2;

  if (missing_in_tr1)
    *missing_in_tr1 = missing1;
  if (missing_in_tr2)
    *missing_in_tr2 = missing2;
  return missing1 + missing2;
}

/******************************************************************************/

BipartitionList* TreeTools::bipartitionOccurrences(const vector<Tree*>& vecTr, vector<size_t>& bipScore)
{
  vector<BipartitionList*> vecBipL;
  BipartitionList* mergedBipL;
  vector<size_t> bipSize;
  size_t nbBip;

  /*  build and merge bipartitions */
  for (size_t i = 0; i < vecTr.size(); i++)
  {
    vecBipL.push_back(new BipartitionList(*vecTr[i]));
  }
  mergedBipL = BipartitionTools::mergeBipartitionLists(vecBipL);
  for (size_t i = 0; i < vecTr.size(); i++)
  {
    delete vecBipL[i];
  }

  mergedBipL->removeTrivialBipartitions();
  nbBip = mergedBipL->getNumberOfBipartitions();
  bipScore.clear();
  for (size_t i = 0; i < nbBip; i++)
  {
    bipSize.push_back(mergedBipL->getPartitionSize(i));
    bipScore.push_back(1);
  }

  /* compare bipartitions */
  for (size_t i = nbBip; i > 0; i--)
  {
    if (bipScore[i - 1] == 0)
      continue;
    for (size_t j = i - 1; j > 0; j--)
    {
      if (bipScore[j - 1] && bipSize[i - 1] == bipSize[j - 1] && mergedBipL->areIdentical(i - 1, j - 1))
      {
        bipScore[i - 1]++;
        bipScore[j - 1] = 0;
      }
    }
  }

  /* keep only distinct bipartitions */
  for (size_t i = nbBip; i > 0; i--)
  {
    if (bipScore[i - 1] == 0)
    {
      bipScore.erase(bipScore.begin() + static_cast<ptrdiff_t>(i - 1));
      mergedBipL->deleteBipartition(i - 1);
    }
  }

  /* add terminal branches */
  mergedBipL->addTrivialBipartitions(false);
  for (size_t i = 0; i < mergedBipL->getNumberOfElements(); i++)
  {
    bipScore.push_back(vecTr.size());
  }

  return mergedBipL;
}

/******************************************************************************/

TreeTemplate<Node>* TreeTools::thresholdConsensus(const vector<Tree*>& vecTr, double threshold, bool checkNames) throw (Exception)
{
  vector<size_t> bipScore;
  vector<string> tr0leaves;
  BipartitionList* bipL;
  double score;

  if (vecTr.size() == 0)
    throw Exception("TreeTools::thresholdConsensus. Empty vector passed");

  /* check names */
  if (checkNames)
  {
    tr0leaves = vecTr[0]->getLeavesNames();
    for (size_t i = 1; i < vecTr.size(); i++)
    {
      if (!VectorTools::haveSameElements(vecTr[i]->getLeavesNames(), tr0leaves))
        throw Exception("TreeTools::thresholdConsensus. Distinct leaf sets between trees");
    }
  }

  bipL = bipartitionOccurrences(vecTr, bipScore);

  for (size_t i = bipL->getNumberOfBipartitions(); i > 0; i--)
  {
    if (bipL->getPartitionSize(i - 1) == 1)
      continue;
    score = static_cast<int>(bipScore[i - 1]) / static_cast<double>(vecTr.size());
    if (score <= threshold && score != 1.)
    {
      bipL->deleteBipartition(i - 1);
      continue;
    }
    if (score > 0.5)
      continue;
    for (size_t j = bipL->getNumberOfBipartitions(); j > i; j--)
    {
      if (!bipL->areCompatible(i - 1, j - 1))
      {
        bipL->deleteBipartition(i - 1);
        break;
      }
    }
  }

  TreeTemplate<Node>* tr = bipL->toTree();
  delete bipL;
  return tr;
}

/******************************************************************************/

TreeTemplate<Node>* TreeTools::fullyResolvedConsensus(const vector<Tree*>& vecTr, bool checkNames)
{
  return thresholdConsensus(vecTr, 0., checkNames);
}

/******************************************************************************/

TreeTemplate<Node>* TreeTools::majorityConsensus(const vector<Tree*>& vecTr, bool checkNames)
{
  return thresholdConsensus(vecTr, 0.5, checkNames);
}

/******************************************************************************/

TreeTemplate<Node>* TreeTools::strictConsensus(const vector<Tree*>& vecTr, bool checkNames)
{
  return thresholdConsensus(vecTr, 1., checkNames);
}

/******************************************************************************/

Tree* TreeTools::MRP(const vector<Tree*>& vecTr)
{
  // matrix representation
  VectorSiteContainer* sites = TreeTools::MRPEncode(vecTr);

  // starting bioNJ tree
  const DNA* alphabet = dynamic_cast<const DNA*>(sites->getAlphabet());
  JCnuc* jc = new JCnuc(alphabet);
  ConstantDistribution* constRate = new ConstantDistribution(1.);
  DistanceEstimation distFunc(jc, constRate, sites, 0, true);
  BioNJ bionjTreeBuilder(false, false);
  bionjTreeBuilder.setDistanceMatrix(*(distFunc.getMatrix()));
  bionjTreeBuilder.computeTree();
  if (ApplicationTools::message)
    ApplicationTools::message->endLine();
  TreeTemplate<Node>* startTree = new TreeTemplate<Node>(*bionjTreeBuilder.getTree());

  // MP optimization
  DRTreeParsimonyScore* MPScore = new DRTreeParsimonyScore(*startTree, *sites, false);
  MPScore = OptimizationTools::optimizeTreeNNI(MPScore, 0);
  delete startTree;
  Tree* retTree = new TreeTemplate<Node>(MPScore->getTree());
  delete MPScore;

  return retTree;
}

/******************************************************************************/

void TreeTools::computeBootstrapValues(Tree& tree, const vector<Tree*>& vecTr, bool verbose, int format)
{
  vector<int> index;
  BipartitionList bpTree(tree, true, &index);
  vector<size_t> occurences;
  BipartitionList* bpList = bipartitionOccurrences(vecTr, occurences);

  vector< Number<double> > bootstrapValues(bpTree.getNumberOfBipartitions());

  for (size_t i = 0; i < bpTree.getNumberOfBipartitions(); i++)
  {
    if (verbose)
      ApplicationTools::displayGauge(i, bpTree.getNumberOfBipartitions() - 1, '=');
    for (size_t j = 0; j < bpList->getNumberOfBipartitions(); j++)
    {
      if (BipartitionTools::areIdentical(bpTree, i, *bpList, j))
      {
        bootstrapValues[i] = format >= 0 ? round(static_cast<double>(occurences[j]) * pow(10., 2 + format) / static_cast<double>(vecTr.size())) / pow(10., format) : static_cast<double>(occurences[j]);
        break;
      }
    }
  }

  for (size_t i = 0; i < index.size(); i++)
  {
    if (!tree.isLeaf(index[i]))
      tree.setBranchProperty(index[i], BOOTSTRAP, bootstrapValues[i]);
  }

  delete bpList;
}

/******************************************************************************/

vector<int> TreeTools::getAncestors(const Tree& tree, int nodeId) throw (NodeNotFoundException)
{
  vector<int> ids;
  int currentId = nodeId;
  while (tree.hasFather(currentId))
  {
    currentId = tree.getFatherId(currentId);
    ids.push_back(currentId);
  }
  return ids;
}

/******************************************************************************/

int TreeTools::getLastCommonAncestor(const Tree& tree, const vector<int>& nodeIds) throw (NodeNotFoundException, Exception)
{
  if (nodeIds.size() == 0)
    throw Exception("TreeTools::getLastCommonAncestor(). You must provide at least one node id.");
  vector< vector<int> > ancestors(nodeIds.size());
  for (size_t i = 0; i < nodeIds.size(); i++)
  {
    ancestors[i] = getAncestors(tree, nodeIds[i]);
    ancestors[i].insert(ancestors[i].begin(), nodeIds[i]);
  }
  int lca = tree.getRootId();
  size_t count = 1;
  for ( ; ; )
  {
    if (ancestors[0].size() <= count)
      return lca;
    int current = ancestors[0][ancestors[0].size() - count - 1];
    for (size_t i = 1; i < nodeIds.size(); i++)
    {
      if (ancestors[i].size() <= count)
        return lca;
      if (ancestors[i][ancestors[i].size() - count - 1] != current)
        return lca;
    }
    lca = current;
    count++;
  }
  // This line is never reached!
  return lca;
}

/******************************************************************************/

void TreeTools::constrainedMidPointRooting(Tree& tree)
{
  // is the tree rooted?
  if (!tree.isRooted())
    throw Exception("The tree has to be rooted on the branch of interest to determine the midpoint position of the root");
  // is the tree multifurcating?
  if (tree.isMultifurcating())
    throw Exception("The tree is multifurcated, which is not allowed.");

  double length = 0.;
  vector<int> sonsIds = tree.getSonsId(tree.getRootId());
  // Length of the branch containing the root:
  length = tree.getDistanceToFather(sonsIds.at(0)) + tree.getDistanceToFather(sonsIds.at(1));
  // The fraction of the original branch allowing to split its length and to place the root:
  double x = bestRootPosition_(tree, sonsIds.at(0), sonsIds.at(1), length);
  // The new branch lengths are then computed:
  tree.setDistanceToFather(sonsIds.at(0), length * x);
  tree.setDistanceToFather(sonsIds.at(1), length * (1 - x));
}

/******************************************************************************/

double TreeTools::bestRootPosition_(Tree& tree, int nodeId1, int nodeId2, double length)
{
  double x;
  Moments_ m1, m2;
  double A, B; // C;
  // The variance is expressed as a degree 2 polynomial : variance(x) = A * x * x + B * x + C
  // The fraction x is then obtained by differentiating this equation.
  m1 = statFromNode_(tree, nodeId1);
  m2 = statFromNode_(tree, nodeId2);
  A = 4 * m1.N * (m2.N * length) * length;
  B = 4 * length * (m2.N * m1.sum - m1.N * m2.sum - length * m1.N * m2.N);
//   C = (m1.N + m2.N) * (m1.squaredSum + m2.squaredSum) + m1.N * length * m2.N * length +
//     2 * m1.N * length * m2.sum - 2 * m2.N * length * m1.sum -
//     (m1.sum + m2.sum) * (m1.sum + m2.sum);

  if (A < 1e-20)
    x = 0.5;
  else
    x = -B / (2 * A);
  if (x < 0)
    x = 0;
  else if (x > 1)
    x = 1;

  return x;
}

/******************************************************************************/

TreeTools::Moments_ TreeTools::statFromNode_(Tree& tree, int rootId)
{
  // This function recursively calculates both the sum of the branch lengths and the sum of the squared branch lengths down the node whose ID is rootId.
  // If below a particular node there are N leaves, the branch between this node and its father is taken into account N times in the calculation.
  Moments_ m;
  static Moments_ mtmp;

  if (tree.isLeaf(rootId))
  {
    m.N = 1;
    m.sum = 0.;
    m.squaredSum = 0.;
  }
  else
  {
    vector<int> sonsId = tree.getSonsId(rootId);
    for (size_t i = 0; i < sonsId.size(); i++)
    {
      mtmp = statFromNode_(tree, sonsId.at(i));
      double bLength = tree.getDistanceToFather(sonsId.at(i));
      m.N += mtmp.N;
      m.sum += mtmp.sum + bLength * mtmp.N;
      m.squaredSum += mtmp.squaredSum + 2 * bLength * mtmp.sum + mtmp.N * bLength * bLength;
    }
  }

  return m;
}

/******************************************************************************/

Tree* TreeTools::MRPMultilabel(const vector<Tree*>& vecTr)
{
    // matrix representation
    VectorSiteContainer* sites = TreeTools::MRPEncode(vecTr);
    
    // starting bioNJ tree
    const DNA* alphabet = dynamic_cast<const DNA*>(sites->getAlphabet());
    JCnuc* jc = new JCnuc(alphabet);
    ConstantDistribution* constRate = new ConstantDistribution(1.);
    DistanceEstimation distFunc(jc, constRate, sites, 0, true);
    BioNJ bionjTreeBuilder(false, false);
    bionjTreeBuilder.setDistanceMatrix(*(distFunc.getMatrix()));
    bionjTreeBuilder.computeTree();
    if (ApplicationTools::message)
        ApplicationTools::message->endLine();
    TreeTemplate<Node>* startTree = new TreeTemplate<Node>(*bionjTreeBuilder.getTree());
    
    // MP optimization
    DRTreeParsimonyScore* MPScore = new DRTreeParsimonyScore(*startTree, *sites, false);
    MPScore = OptimizationTools::optimizeTreeNNI(MPScore, 0);
    delete startTree;
    Tree* retTree = new TreeTemplate<Node>(MPScore->getTree());
    delete MPScore;
    
    return retTree;
}


