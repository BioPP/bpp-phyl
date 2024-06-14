// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "PhyloStatistics.h"
#include "Tree/Node.h"
#include "Tree/TreeTemplate.h"

using namespace bpp;
using namespace std;

void PhyloStatistics::setTree(const Tree& tree)
{
  const Tree* treeP = &tree;
  const TreeTemplate<Node>* ttreeP = dynamic_cast<const TreeTemplate<Node>*>(treeP);
  bool copy = false;
  if (!ttreeP)
  {
    ttreeP = new TreeTemplate<Node>(tree);
    copy = true;
  }

  // Perform computations:
  numberOfLeaves_ = 0;
  numberOfAncestors_ = 0;
  branchLengths_.clear();
  nodeHeights_.clear();
  nodeDepths_.clear();
  nodeNumberOfSons_.clear();
  nodeIds_.clear();
  double h;
  size_t d;
  computeForSubtree_(ttreeP->getRootNode(), h, d);

  if (copy)
    delete ttreeP;
}

void PhyloStatistics::computeForSubtree_(const Node* node, double& height, size_t& depth)
{
  if (node->isLeaf())
    numberOfLeaves_++;
  else
    numberOfAncestors_++;

  nodeNumberOfSons_.push_back(node->getNumberOfSons());

  height = 0;
  depth = 0;
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    const Node* son = (*node)[static_cast<int>(i)];
    double dist = 0;
    if (son->hasDistanceToFather())
      dist = son->getDistanceToFather();
    else
      dist = 0;
    double subHeight;
    size_t subDepth;
    computeForSubtree_(son, subHeight, subDepth);
    subHeight += dist;
    subDepth++;
    if (subHeight > height)
      height = subHeight;
    if (subDepth  > depth)
      depth  = subDepth;
  }

  if (node->hasDistanceToFather())
    branchLengths_.push_back(node->getDistanceToFather());
  else
    branchLengths_.push_back(log(0)); // -Inf if no branch length.

  nodeHeights_.push_back(height);
  nodeDepths_.push_back(depth);
}
