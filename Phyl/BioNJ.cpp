//
// File: BioNJ.h
// Created by: Vincent Ranwez
// Created on: Tue Apr 11 14:23 2006
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004, 2005, 2006)

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
#include "Tree.h"

//From Utils:
#include <Utils/ApplicationTools.h>

using namespace bpp;

// From the STL:
#include <cmath>
#include <iostream>

using namespace std;

double BioNJ::computeDistancesFromPair(const vector<unsigned int> & pair, const vector<double> & branchLengths, unsigned int pos)
{
  return _positiveLengths ?
    std::max(_lambda * (_matrix(pair[0], pos) - branchLengths[0]) + (1-_lambda)*(_matrix(pair[1], pos) - branchLengths[1]), 0.)
  :          _lambda * (_matrix(pair[0], pos) - branchLengths[0]) + (1-_lambda)*(_matrix(pair[1], pos) - branchLengths[1]); 
}

void BioNJ::computeTree(bool rooted) throw (Exception)
{
  // Initialization:
  for(unsigned int i = 0; i < _matrix.size(); i++)
  {
    _currentNodes[i] = getLeafNode(i, _matrix.getName(i));
  }
  unsigned int idNextNode = _matrix.size();
  vector<double> newDist(_matrix.size());
  vector<double> newVar(_matrix.size());
  
  // Build tree:
  while(_currentNodes.size() > (rooted ? 2 : 3))
  {
    if(_verbose)
      ApplicationTools::displayGauge(_matrix.size() - _currentNodes.size(), _matrix.size() - (rooted ? 2 : 3) - 1);
    vector<unsigned int> bestPair = getBestPair();
    vector<double> distances = computeBranchLengthsForPair(bestPair);
    Node * best1 = _currentNodes[bestPair[0]];
    Node * best2 = _currentNodes[bestPair[1]];
    // Distances may be used by getParentNodes (PGMA for instance).
    best1->setDistanceToFather(distances[0]);
    best2->setDistanceToFather(distances[1]);
    Node * parent = getParentNode(idNextNode++, best1, best2);
    // compute lambda
    _lambda=0;
    if(_variance[bestPair[0]][bestPair[1]] == 0)
      _lambda=.5;
    else
    {
      for(map<unsigned int, Node *>::iterator i = _currentNodes.begin(); i != _currentNodes.end(); i++)
      {
        unsigned int id = i -> first;
        if(id != bestPair[0] && id != bestPair[1]) 
          _lambda += (_variance[bestPair[1]][id] - _variance[bestPair[0]][id]);
      }
      double div =  2*(_currentNodes.size() - 2)*_variance[bestPair[0]][bestPair[1]];
      _lambda /= div;
      _lambda += .5;
    }
    if(_lambda <0.) _lambda =0.;
    if(_lambda > 1.) _lambda=1.;
    
    for(map<unsigned int, Node *>::iterator i = _currentNodes.begin(); i != _currentNodes.end(); i++)
    {
      unsigned int id = i -> first;
      if(id != bestPair[0] && id != bestPair[1])
      {
        newDist[id] = computeDistancesFromPair(bestPair, distances, id);
        newVar[id] = _lambda*_variance[bestPair[0]][id]+ (1-_lambda)*_variance[bestPair[1]][id]- _lambda*(1-_lambda)*_variance[bestPair[0]][bestPair[1]];
      }
      else
      {
        newDist[id] = 0;
      }
    }
    // Actualize _currentNodes:
    _currentNodes[bestPair[0]] = parent;
    _currentNodes.erase(bestPair[1]);
    for(map<unsigned int, Node *>::iterator i = _currentNodes.begin(); i != _currentNodes.end(); i++)
    {
      unsigned int id = i -> first;
      _matrix(  bestPair[0], id) =    _matrix(id, bestPair[0]) = newDist[id];
      _variance(bestPair[0], id) =  _variance(id, bestPair[0]) = newVar[id];
    }
    
  }
  finalStep(idNextNode);
}

