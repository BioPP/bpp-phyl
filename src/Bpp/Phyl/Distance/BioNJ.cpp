// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/App/ApplicationTools.h>

#include "../Tree/Tree.h"
#include "BioNJ.h"

using namespace bpp;

// From the STL:
#include <cmath>
#include <iostream>

using namespace std;

double BioNJ::computeDistancesFromPair(const vector<size_t>& pair, const vector<double>& branchLengths, size_t pos)
{
  return positiveLengths_ ?
         std::max(lambda_ * (matrix_(pair[0], pos) - branchLengths[0]) + (1 - lambda_) * (matrix_(pair[1], pos) - branchLengths[1]), 0.)
         :          lambda_* (matrix_(pair[0], pos) - branchLengths[0]) + (1 - lambda_) * (matrix_(pair[1], pos) - branchLengths[1]);
}

void BioNJ::computeTree()
{
  // Initialization:
  for (size_t i = 0; i < matrix_.size(); i++)
  {
    currentNodes_[i] = getLeafNode(static_cast<int>(i), matrix_.getName(i));
  }
  int idNextNode = static_cast<int>(matrix_.size());
  vector<double> newDist(matrix_.size());
  vector<double> newVar(matrix_.size());

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
    Node* parent = getParentNode(idNextNode++, best1, best2);
    // compute lambda
    lambda_ = 0;
    if (variance_(bestPair[0], bestPair[1]) == 0)
      lambda_ = .5;
    else
    {
      for (map<size_t, Node*>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++)
      {
        size_t id = i->first;
        if (id != bestPair[0] && id != bestPair[1])
          lambda_ += (variance_(bestPair[1], id) - variance_(bestPair[0], id));
      }
      double div = 2 * static_cast<double>(currentNodes_.size() - 2) * variance_(bestPair[0], bestPair[1]);
      lambda_ /= div;
      lambda_ += .5;
    }
    if (lambda_ < 0.)
      lambda_ = 0.;
    if (lambda_ > 1.)
      lambda_ = 1.;

    for (map<size_t, Node*>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++)
    {
      size_t id = i->first;
      if (id != bestPair[0] && id != bestPair[1])
      {
        newDist[id] = computeDistancesFromPair(bestPair, distances, id);
        newVar[id] = lambda_ * variance_(bestPair[0], id) + (1 - lambda_) * variance_(bestPair[1], id) - lambda_ * (1 - lambda_) * variance_(bestPair[0], bestPair[1]);
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
      matrix_(  bestPair[0], id) =    matrix_(id, bestPair[0]) = newDist[id];
      variance_(bestPair[0], id) =  variance_(id, bestPair[0]) = newVar[id];
    }
  }
  finalStep(idNextNode);
}
