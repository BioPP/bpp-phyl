// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>

#include "AwareNode.h"

using namespace bpp;
using namespace std;

/** Copy constructor: *********************************************************/

AwareNode::AwareNode(const AwareNode& node) :
  id_(node.id_),
  sons_(), father_(),
  distanceToFather_(node.distanceToFather_)
{}

/** Assignation operator: *****************************************************/

AwareNode& AwareNode::operator=(const AwareNode& node)
{
  id_               = node.id_;
  father_           = 0;
  distanceToFather_ = node.distanceToFather_;
  sons_.clear();

  return *this;
}

/** Sons: *********************************************************************/

void AwareNode::swap(size_t branch1, size_t branch2)
{
  if (branch1 > branch2)
  {
    size_t tmp = branch1;
    branch1 = branch2;
    branch2 = tmp;
  }
  AwareNode* node1 = getSon(branch1);
  AwareNode* node2 = getSon(branch2);
  removeSon(node1);
  removeSon(node2);
  addSon(branch1, node2);
  addSon(branch2, node1);
}

size_t AwareNode::getSonPosition(const AwareNode* son) const
{
  if (!son)
    throw NullPointerException("Node::getSonPosition(). Empty node given as input.");
  for (size_t i = 0; i < sons_.size(); i++)
  {
    if (sons_[i] == son)
      return i;
  }
  throw Exception("Son not found: " + TextTools::toString(son->getId()));
}

/******************************************************************************/
