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
#include "Tree.h"

using namespace bpp;

#include <cmath>
#include <iostream>

using namespace std;

vector<unsigned int> NeighborJoining::getBestPair() throw (Exception)
{
  for(map<unsigned int, Node *>::iterator i = _currentNodes.begin(); i != _currentNodes.end(); i++) {
    unsigned int id = i -> first;
    _sumDist[id] = 0;
    for(map<unsigned int, Node *>::iterator j = _currentNodes.begin(); j != _currentNodes.end(); j++) {
      unsigned int jd = j -> first;
      _sumDist[id] += _matrix(id, jd);
    }
  }

  vector<unsigned int> bestPair(2);
  double critMax = std::log(0.);
  for(map<unsigned int, Node *>::iterator i = _currentNodes.begin(); i != _currentNodes.end(); i++) {
    unsigned int id = i -> first;
    map<unsigned int, Node *>::iterator j = i;
    j++;
    for(; j != _currentNodes.end(); j++) {
      unsigned int jd = j -> first;
      double crit = _sumDist[id] + _sumDist[jd] - (_currentNodes.size() - 2) * _matrix(id, jd);
      //cout << "\t" << id << "\t" << jd << "\t" << crit << endl;
      if(crit > critMax) {
        critMax = crit;
        bestPair[0] = id;
        bestPair[1] = jd;
      }
    }
  }

  if(critMax == std::log(0.)) {
    throw Exception("Unexpected error: no maximum criterium found.");
  }
  return bestPair;  
}

vector<double> NeighborJoining::computeBranchLengthsForPair(const vector<unsigned int> & pair)
{
  double ratio = (_sumDist[pair[0]] - _sumDist[pair[1]]) / (_currentNodes.size() - 2);
  vector<double> d(2);
  if(_positiveLengths) {
    d[0] = std::max(.5 * (_matrix(pair[0], pair[1]) + ratio), 0.); 
    d[1] = std::max(.5 * (_matrix(pair[0], pair[1]) - ratio), 0.); 
  } else {
    d[0] = .5 * (_matrix(pair[0], pair[1]) + ratio); 
    d[1] = .5 * (_matrix(pair[0], pair[1]) - ratio); 
  }
  return d;
}

double NeighborJoining::computeDistancesFromPair(const vector<unsigned int> & pair, const vector<double> & branchLengths, unsigned int pos)
{
  return 
    _positiveLengths ?
      std::max(.5 * (_matrix(pair[0], pos) - branchLengths[0] + _matrix(pair[1], pos) - branchLengths[1]), 0.)
    :          .5 * (_matrix(pair[0], pos) - branchLengths[0] + _matrix(pair[1], pos) - branchLengths[1]); 
}

void NeighborJoining::finalStep(int idRoot)
{
  Node * root = new Node(idRoot);
  map<unsigned int, Node* >::iterator it = _currentNodes.begin();
  unsigned int i1 = it->first;
  Node * n1       = it->second;
  it++;
  unsigned int i2 = it->first;
  Node * n2       = it->second;
  if(_currentNodes.size() == 2)
  {
    //Rooted
    double d = _matrix(i1, i2) / 2;
    root->addSon(*n1);
    root->addSon(*n2);
    n1->setDistanceToFather(d);
    n2->setDistanceToFather(d);
  }
  else
  {
    //Unrooted
    it++;
    unsigned int i3 = it->first;
    Node * n3       = it->second;
    double d1 = _positiveLengths ?
        std::max(_matrix(i1, i2) + _matrix(i1, i3) - _matrix(i2, i3), 0.)
      :          _matrix(i1, i2) + _matrix(i1, i3) - _matrix(i2, i3);
    double d2 = _positiveLengths ?
        std::max(_matrix(i2, i1) + _matrix(i2, i3) - _matrix(i1, i3), 0.)
      :          _matrix(i2, i1) + _matrix(i2, i3) - _matrix(i1, i3);
    double d3 = _positiveLengths ?
        std::max(_matrix(i3, i1) + _matrix(i3, i2) - _matrix(i1, i2), 0.)
      :          _matrix(i3, i1) + _matrix(i3, i2) - _matrix(i1, i2);
    root->addSon(*n1);
    root->addSon(*n2);
    root->addSon(*n3);
    n1->setDistanceToFather(d1/2.);
    n2->setDistanceToFather(d2/2.);
    n3->setDistanceToFather(d3/2.);
  }
  _tree = new TreeTemplate<Node>(*root);
}

