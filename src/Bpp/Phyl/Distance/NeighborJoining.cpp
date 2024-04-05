// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../Tree/Tree.h"
#include "NeighborJoining.h"

using namespace bpp;

#include <cmath>
#include <iostream>

using namespace std;

std::vector<size_t> NeighborJoining::getBestPair()
{
  for (std::map<size_t, Node*>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++)
  {
    size_t id = i->first;
    sumDist_[id] = 0;
    for (map<size_t, Node*>::iterator j = currentNodes_.begin(); j != currentNodes_.end(); j++)
    {
      size_t jd = j->first;
      sumDist_[id] += matrix_(id, jd);
    }
  }

  vector<size_t> bestPair(2);
  double critMax = std::log(0.);
  for (map<size_t, Node*>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++)
  {
    size_t id = i->first;
    map<size_t, Node*>::iterator j = i;
    j++;
    for ( ; j != currentNodes_.end(); j++)
    {
      size_t jd = j->first;
      double crit = sumDist_[id] + sumDist_[jd] - static_cast<double>(currentNodes_.size() - 2) * matrix_(id, jd);
      // cout << "\t" << id << "\t" << jd << "\t" << crit << endl;
      if (crit > critMax)
      {
        critMax = crit;
        bestPair[0] = id;
        bestPair[1] = jd;
      }
    }
  }

  if (critMax == std::log(0.))
  {
    throw Exception("Unexpected error: no maximum criterium found.");
  }
  return bestPair;
}

std::vector<double> NeighborJoining::computeBranchLengthsForPair(const std::vector<size_t>& pair)
{
  double ratio = (sumDist_[pair[0]] - sumDist_[pair[1]]) / static_cast<double>(currentNodes_.size() - 2);
  vector<double> d(2);
  if (positiveLengths_)
  {
    d[0] = std::max(.5 * (matrix_(pair[0], pair[1]) + ratio), 0.);
    d[1] = std::max(.5 * (matrix_(pair[0], pair[1]) - ratio), 0.);
  }
  else
  {
    d[0] = .5 * (matrix_(pair[0], pair[1]) + ratio);
    d[1] = .5 * (matrix_(pair[0], pair[1]) - ratio);
  }
  return d;
}

double NeighborJoining::computeDistancesFromPair(const std::vector<size_t>& pair, const std::vector<double>& branchLengths, size_t pos)
{
  return
    positiveLengths_ ?
    std::max(.5 * (matrix_(pair[0], pos) - branchLengths[0] + matrix_(pair[1], pos) - branchLengths[1]), 0.)
    :          .5 * (matrix_(pair[0], pos) - branchLengths[0] + matrix_(pair[1], pos) - branchLengths[1]);
}

void NeighborJoining::finalStep(int idRoot)
{
  Node* root = new Node(idRoot);
  map<size_t, Node* >::iterator it = currentNodes_.begin();
  size_t i1 = it->first;
  Node* n1       = it->second;
  it++;
  size_t i2 = it->first;
  Node* n2       = it->second;
  if (currentNodes_.size() == 2)
  {
    // Rooted
    double d = matrix_(i1, i2) / 2;
    root->addSon(n1);
    root->addSon(n2);
    n1->setDistanceToFather(d);
    n2->setDistanceToFather(d);
  }
  else
  {
    // Unrooted
    it++;
    size_t i3 = it->first;
    Node* n3       = it->second;
    double d1 = positiveLengths_ ?
        std::max(matrix_(i1, i2) + matrix_(i1, i3) - matrix_(i2, i3), 0.)
                :          matrix_(i1, i2) + matrix_(i1, i3) - matrix_(i2, i3);
    double d2 = positiveLengths_ ?
        std::max(matrix_(i2, i1) + matrix_(i2, i3) - matrix_(i1, i3), 0.)
                :          matrix_(i2, i1) + matrix_(i2, i3) - matrix_(i1, i3);
    double d3 = positiveLengths_ ?
        std::max(matrix_(i3, i1) + matrix_(i3, i2) - matrix_(i1, i2), 0.)
                :          matrix_(i3, i1) + matrix_(i3, i2) - matrix_(i1, i2);
    root->addSon(n1);
    root->addSon(n2);
    root->addSon(n3);
    n1->setDistanceToFather(d1 / 2.);
    n2->setDistanceToFather(d2 / 2.);
    n3->setDistanceToFather(d3 / 2.);
  }
  tree_.reset(new TreeTemplate<Node>(root));
}
