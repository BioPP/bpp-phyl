// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <iostream>

#include "PhyloNode.h"

using namespace bpp;
using namespace std;

const NodeEvent NodeEvent::speciationEvent(NodeEvent::NodeType::Speciation);
const NodeEvent NodeEvent::mixtureEvent(NodeEvent::NodeType::Mixture);

/** Copy constructor: *********************************************************/

PhyloNode::PhyloNode(const PhyloNode& node) :
  name_(node.name_),
  properties_()
{
  for (map<string, Clonable*>::iterator i = node.properties_.begin(); i != node.properties_.end(); i++)
  {
    properties_[i->first] = i->second->clone();
  }
}

/** Assignation operator: *****************************************************/

PhyloNode& PhyloNode::operator=(const PhyloNode& node)
{
  name_  = node.name_;
  for (map<string, Clonable*>::iterator i = node.properties_.begin(); i != node.properties_.end(); i++)
  {
    Clonable* p = properties_[i->first];
    if (p)
      delete p;
    properties_[i->first] = i->second->clone();
  }

  return *this;
}

/******************************************************************************/
