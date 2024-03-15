// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/App/ApplicationTools.h>

#include "../Tree/Node.h"
#include "AbstractAgglomerativeDistanceMethod.h"

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

void AbstractAgglomerativeDistanceMethod::setDistanceMatrix(const DistanceMatrix& matrix)
{
  if (matrix.size() <= 3)
    throw Exception("AbstractAgglomerativeDistanceMethod::setDistanceMatrix(): matrix must be at least of dimension 3.");
  matrix_ = matrix;
  currentNodes_.clear();
  tree_.reset(nullptr);
}

void AbstractAgglomerativeDistanceMethod::computeTree()
{
  // Initialization:
  for (size_t i = 0; i < matrix_.size(); ++i)
  {
    currentNodes_[i] = getLeafNode(static_cast<int>(i), matrix_.getName(i));
  }
  int idNextNode = static_cast<int>(matrix_.size());
  vector<double> newDist(matrix_.size());

  // Build tree:
  while (currentNodes_.size() > (rootTree_ ? 2 : 3))
  {
    if (verbose_)
      ApplicationTools::displayGauge(matrix_.size() - currentNodes_.size(), matrix_.size() - (rootTree_ ? 2 : 3) - 1);
    vector<size_t> bestPair = getBestPair();
    vector<double> distances = computeBranchLengthsForPair(bestPair);
    Node* best1 = currentNodes_[bestPair[0]];
    Node* best2 = currentNodes_[bestPair[1]];
    // Distances may be used by getParentNodes (PGMA for instance).
    best1->setDistanceToFather(distances[0]);
    best2->setDistanceToFather(distances[1]);
    Node* parent = getParentNode(idNextNode, best1, best2);
    idNextNode++;
    for (map<size_t, Node*>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++)
    {
      size_t id = i->first;
      if (id != bestPair[0] && id != bestPair[1])
      {
        assert (id < newDist.size()); // DEBUG
        newDist[id] = computeDistancesFromPair(bestPair, distances, id);
      }
      else
      {
        newDist[id] = 0;
      }
    }
    // Actualize currentNodes_:
    currentNodes_[bestPair[0]] = parent;
    currentNodes_.erase(bestPair[1]);
    for (map<size_t, Node*>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++)
    {
      size_t id = i->first;
      matrix_(bestPair[0], id) = matrix_(id, bestPair[0]) = newDist[id];
    }
  }
  finalStep(idNextNode);
}

Node* AbstractAgglomerativeDistanceMethod::getLeafNode(int id, const std::string& name)
{
  return new Node(id, name);
}

Node* AbstractAgglomerativeDistanceMethod::getParentNode(int id, Node* son1, Node* son2)
{
  Node* parent = new Node(id);
  parent->addSon(son1);
  parent->addSon(son2);
  return parent;
}
