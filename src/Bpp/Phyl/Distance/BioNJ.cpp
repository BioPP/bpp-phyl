//
// File: BioNJ.h
// Created by: Vincent Ranwez
// Created on: Tue Apr 11 14:23 2006
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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

#include "BioNJ.h"
#include "../Tree/Tree.h"

#include <Bpp/App/ApplicationTools.h>

using namespace bpp;

// From the STL:
#include <cmath>
#include <iostream>

using namespace std;

double BioNJ::computeDistancesFromPair(const vector<size_t>& pair, const vector<double>& branchLengths, size_t pos)
{
  return positiveLengths_ ?
         std::max(lambda_ * (matrix_(pair[0], pos) - branchLengths[0]) + (1 - lambda_) * (matrix_(pair[1], pos) - branchLengths[1]), 0.)
         :          lambda_ * (matrix_(pair[0], pos) - branchLengths[0]) + (1 - lambda_) * (matrix_(pair[1], pos) - branchLengths[1]);
}

void BioNJ::computeTree() throw (Exception)
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

