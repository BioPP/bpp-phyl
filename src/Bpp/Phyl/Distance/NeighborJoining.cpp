//
// File: NeighborJoining.cpp
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Thu jun 23 10:39 2005
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

#include "NeighborJoining.h"
#include "../Tree/Tree.h"

using namespace bpp;

#include <cmath>
#include <iostream>

using namespace std;

std::vector<size_t> NeighborJoining::getBestPair() throw (Exception)
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
  tree_ = new TreeTemplate<Node>(root);
}

