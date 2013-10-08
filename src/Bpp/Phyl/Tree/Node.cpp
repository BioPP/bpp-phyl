//
// File: Node.cpp
// Created by: Julien Dutheil
// Created on: Thu Mar 13 12:03:18 2003
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

#include "Node.h"
#include "TreeTools.h"

#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>

using namespace bpp;

//from the STL:
#include <algorithm>
#include <iostream>

using namespace std;

/** Copy constructor: *********************************************************/
  
Node::Node(const Node& node):
  id_(node.id_), name_(0),
  sons_(), father_(0),
  //, sons_(node.sons_), father_(node.father_),
  distanceToFather_(0), nodeProperties_(), branchProperties_()
{
  name_             = node.hasName() ? new string(* node.name_) : 0;
  distanceToFather_ = node.hasDistanceToFather() ? new double(* node.distanceToFather_) : 0;
  for (map<string, Clonable *>::iterator i = node.nodeProperties_.begin(); i != node.nodeProperties_.end(); i++)
    nodeProperties_[i->first] = i->second->clone();
  for (map<string, Clonable *>::iterator i = node.branchProperties_.begin(); i != node.branchProperties_.end(); i++)
    branchProperties_[i->first] = i->second->clone();
}

/** Assignation operator: *****************************************************/

Node& Node::operator=(const Node & node)
{
  id_               = node.id_;
  if(name_) delete name_;
  name_             = node.hasName() ? new string(* node.name_) : 0;
  //father_           = node.father_;
  if(distanceToFather_) delete distanceToFather_;
  distanceToFather_ = node.hasDistanceToFather() ? new double(* node.distanceToFather_) : 0;
  //sons_             = node.sons_;
  for(map<string, Clonable *>::iterator i = node.nodeProperties_.begin(); i != node.nodeProperties_.end(); i++)
  {
    Clonable * p = nodeProperties_[i->first];
    if(p) delete p;
    nodeProperties_[i->first] = i->second->clone();
  }
  for(map<string, Clonable *>::iterator i = node.branchProperties_.begin(); i != node.branchProperties_.end(); i++)
  {
    Clonable * p = branchProperties_[i->first];
    if(p) delete p;
    branchProperties_[i->first] = i->second->clone();
  }
  return * this;
}
      
/** Sons: *********************************************************************/
      
void Node::swap(size_t branch1, size_t branch2) throw (IndexOutOfBoundsException)
{
  if (branch1 > branch2)
  {
    size_t tmp = branch1;
    branch1 = branch2;
    branch2 = tmp;
  }
  Node* node1 = getSon(branch1);
  Node* node2 = getSon(branch2);
  removeSon(node1);
  removeSon(node2);
  addSon(branch1, node2);
  addSon(branch2, node1);
}

vector<const Node *> Node::getNeighbors() const
{
  vector<const Node *> neighbors;
  if(hasFather()) neighbors.push_back(father_);
  for(size_t i = 0; i < sons_.size(); i++) neighbors.push_back(sons_[i]);
  return neighbors;
}
    
vector<Node *> Node::getNeighbors()
{
  vector<Node *> neighbors;
  if(hasFather()) neighbors.push_back(father_);
  for(size_t i = 0; i < sons_.size(); i++) neighbors.push_back(sons_[i]);
  return neighbors;
}

size_t Node::getSonPosition(const Node* son) const throw (NodeNotFoundException, NullPointerException)
{
  if (!son)
    throw NullPointerException("Node::getSonPosition(). Empty node given as input.");
  for(size_t i = 0; i < sons_.size(); i++)
  {
    if(sons_[i] == son) return i;
  }
  throw NodeNotFoundException("Son not found", TextTools::toString(son->getId()));
}

bool Node::hasBootstrapValue() const
{
  return hasBranchProperty(TreeTools::BOOTSTRAP);
}

double Node::getBootstrapValue() const throw (PropertyNotFoundException)
{
  if(hasBranchProperty(TreeTools::BOOTSTRAP))
    return dynamic_cast<const Number<double> *>(getBranchProperty(TreeTools::BOOTSTRAP))->getValue();
  else
    throw PropertyNotFoundException("", TreeTools::BOOTSTRAP, this);
}

/******************************************************************************/

