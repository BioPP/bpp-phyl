// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>

#include "Node.h"
#include "TreeTools.h"

using namespace bpp;

// from the STL:
#include <algorithm>
#include <iostream>

using namespace std;

/** Copy constructor: *********************************************************/

Node::Node(const Node& node) :
  id_(node.id_), name_(0),
  sons_(), father_(0),
  // , sons_(node.sons_), father_(node.father_),
  distanceToFather_(0), nodeProperties_(), branchProperties_()
{
  name_             = node.hasName() ? new string(*node.name_) : 0;
  distanceToFather_ = node.hasDistanceToFather() ? new double(*node.distanceToFather_) : 0;
  for (map<string, Clonable*>::iterator i = node.nodeProperties_.begin(); i != node.nodeProperties_.end(); i++)
  {
    nodeProperties_[i->first] = i->second->clone();
  }
  for (map<string, Clonable*>::iterator i = node.branchProperties_.begin(); i != node.branchProperties_.end(); i++)
  {
    branchProperties_[i->first] = i->second->clone();
  }
}

/** Assignation operator: *****************************************************/

Node& Node::operator=(const Node& node)
{
  id_               = node.id_;
  if (name_)
    delete name_;
  name_             = node.hasName() ? new string(*node.name_) : 0;
  // father_           = node.father_;
  if (distanceToFather_)
    delete distanceToFather_;
  distanceToFather_ = node.hasDistanceToFather() ? new double(*node.distanceToFather_) : 0;
  // sons_             = node.sons_;
  for (map<string, Clonable*>::iterator i = node.nodeProperties_.begin(); i != node.nodeProperties_.end(); i++)
  {
    Clonable* p = nodeProperties_[i->first];
    if (p)
      delete p;
    nodeProperties_[i->first] = i->second->clone();
  }
  for (map<string, Clonable*>::iterator i = node.branchProperties_.begin(); i != node.branchProperties_.end(); i++)
  {
    Clonable* p = branchProperties_[i->first];
    if (p)
      delete p;
    branchProperties_[i->first] = i->second->clone();
  }
  return *this;
}

/** Sons: *********************************************************************/

void Node::swap(size_t branch1, size_t branch2)
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

vector<const Node*> Node::getNeighbors() const
{
  vector<const Node*> neighbors;
  if (hasFather())
    neighbors.push_back(father_);
  for (size_t i = 0; i < sons_.size(); i++)
  {
    neighbors.push_back(sons_[i]);
  }
  return neighbors;
}

vector<Node*> Node::getNeighbors()
{
  vector<Node*> neighbors;
  if (hasFather())
    neighbors.push_back(father_);
  for (size_t i = 0; i < sons_.size(); i++)
  {
    neighbors.push_back(sons_[i]);
  }
  return neighbors;
}

size_t Node::getSonPosition(const Node* son) const
{
  if (!son)
    throw NullPointerException("Node::getSonPosition(). Empty node given as input.");
  for (size_t i = 0; i < sons_.size(); i++)
  {
    if (sons_[i] == son)
      return i;
  }
  throw NodeNotFoundException("Son not found", TextTools::toString(son->getId()));
}

bool Node::hasBootstrapValue() const
{
  return hasBranchProperty(TreeTools::BOOTSTRAP);
}

double Node::getBootstrapValue() const
{
  if (hasBranchProperty(TreeTools::BOOTSTRAP))
    return dynamic_cast<const Number<double>*>(getBranchProperty(TreeTools::BOOTSTRAP))->getValue();
  else
    throw PropertyNotFoundException("", TreeTools::BOOTSTRAP, this);
}

/******************************************************************************/
