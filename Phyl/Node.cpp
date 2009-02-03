//
// File: Node.cpp
// Created by: Julien Dutheil
// Created on: Thu Mar 13 12:03:18 2003
//

/*
Copyright or © or Copr. Julien Dutheil, (November 16, 2004)

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

#include "Node.h"
#include "TreeTools.h"

//From Utils:
#include <Utils/Exceptions.h>
#include <Utils/TextTools.h>

using namespace bpp;

//from the STL:
#include <algorithm>
#include <iostream>

using namespace std;

/** Copy constructor: *********************************************************/
  
Node::Node(const Node & node):
  _id(node._id), _name(NULL), _sons(node._sons), _father(node._father),
  _distanceToFather(NULL), _nodeProperties(), _branchProperties()
{
  _name             = node.hasName() ? new string(* node._name) : NULL;
  _distanceToFather = node.hasDistanceToFather() ? new double(* node._distanceToFather) : NULL;
  for(map<string, Clonable *>::iterator i = node._nodeProperties.begin(); i != node._nodeProperties.end(); i++)
    _nodeProperties[i->first] = i->second->clone();
  for(map<string, Clonable *>::iterator i = node._branchProperties.begin(); i != node._branchProperties.end(); i++)
    _branchProperties[i->first] = i->second->clone();
}

/** Assignation operator: *****************************************************/

Node & Node::operator=(const Node & node)
{
  _id               = node._id;
  if(_name) delete _name;
  _name             = node.hasName() ? new string(* node._name) : NULL;
  _father           = node._father;
  if(_distanceToFather) delete _distanceToFather;
  _distanceToFather = node.hasDistanceToFather() ? new double(* node._distanceToFather) : NULL;
  _sons             = node._sons;
  for(map<string, Clonable *>::iterator i = node._nodeProperties.begin(); i != node._nodeProperties.end(); i++)
  {
    Clonable * p = _nodeProperties[i->first];
    if(p) delete p;
    _nodeProperties[i->first] = i->second->clone();
  }
  for(map<string, Clonable *>::iterator i = node._branchProperties.begin(); i != node._branchProperties.end(); i++)
  {
    Clonable * p = _branchProperties[i->first];
    if(p) delete p;
    _branchProperties[i->first] = i->second->clone();
  }
  return * this;
}
      
/** Sons: *********************************************************************/
      
void Node::swap(unsigned int branch1, unsigned int branch2) throw (IndexOutOfBoundsException)
{
    Node* node1 = getSon(branch1);
    Node* node2 = getSon(branch2);
    removeSon(*node1);
    removeSon(*node2);
    addSon(branch1, *node2);
    addSon(branch2, *node1);
}

vector<const Node *> Node::getNeighbors() const
{
  vector<const Node *> neighbors;
  if(hasFather()) neighbors.push_back(_father);
  for(unsigned int i = 0; i < _sons.size(); i++) neighbors.push_back(_sons[i]);
  return neighbors;
}
    
vector<Node *> Node::getNeighbors()
{
  vector<Node *> neighbors;
  if(hasFather()) neighbors.push_back(_father);
  for(unsigned int i = 0; i < _sons.size(); i++) neighbors.push_back(_sons[i]);
  return neighbors;
}

unsigned int Node::getSonPosition(const Node & son) const throw (NodeNotFoundException)
{
  for(unsigned int i = 0; i < _sons.size(); i++) {
    if(_sons[i] == &son) return i;
  }
  throw NodeNotFoundException("Son not found", TextTools::toString(son.getId()));
}

double Node::getBootstrapValue() const throw (PropertyNotFoundException)
{
  if(hasBranchProperty(TreeTools::BOOTSTRAP))
    return dynamic_cast<const Number<double> *>(getBranchProperty(TreeTools::BOOTSTRAP))->getValue();
  else
    throw PropertyNotFoundException("", TreeTools::BOOTSTRAP, this);
}

/******************************************************************************/

