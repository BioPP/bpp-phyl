//
// File: PHyloStatistics.cpp
// Created by: Julien Dutheil
// Created on: Sat Aug 08 07:29 2009
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

#include "PhyloStatistics.h"
#include "Tree/TreeTemplate.h"
#include "Tree/Node.h"

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

  //Perform computations:
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

  if (copy) delete ttreeP;
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
    if (son->hasDistanceToFather()) dist = son->getDistanceToFather();
    else dist = 0;
    double subHeight;
    size_t subDepth;
    computeForSubtree_(son, subHeight, subDepth);
    subHeight += dist;
    subDepth++;
    if (subHeight > height) height = subHeight;
    if (subDepth  > depth ) depth  = subDepth ;
  }

  if (node->hasDistanceToFather())
    branchLengths_.push_back(node->getDistanceToFather());
  else 
    branchLengths_.push_back(log(0)); //-Inf if no branch length.

  nodeHeights_.push_back(height);
  nodeDepths_.push_back(depth);
}


