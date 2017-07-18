//
// File: AbstractAgglomerativeDistanceMethod.cpp
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Wed jun 22 10:00 2005
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

#include "AbstractAgglomerativeDistanceMethod.h"
#include "../Tree/Node.h"

#include <Bpp/App/ApplicationTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

void AbstractAgglomerativeDistanceMethod::setDistanceMatrix(const DistanceMatrix& matrix) throw (Exception)
{
  if (matrix.size() <= 3)
    throw Exception("AbstractAgglomerativeDistanceMethod::setDistanceMatrix(): matrix must be at least of dimension 3.");
  matrix_ = matrix;
  currentNodes_.clear();
  if (tree_) delete tree_;
}
    
void AbstractAgglomerativeDistanceMethod::computeTree() throw (Exception)
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
    for (map<size_t, Node *>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++)
    {
      size_t id = i->first;
      if (id != bestPair[0] && id != bestPair[1])
      {
        assert (id < newDist.size()); //DEBUG
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
    for (map<size_t, Node *>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++)
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

