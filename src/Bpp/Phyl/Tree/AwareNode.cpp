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

#include "AwareNode.h"
#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>

using namespace bpp;
using namespace std;

/** Copy constructor: *********************************************************/
  
AwareNode::AwareNode(const AwareNode& node):
  id_(node.id_),
  sons_(), father_(),
  distanceToFather_(node.distanceToFather_)
{
}

AwareNode::AwareNode(const PhyloNode& node):
  id_(),
  sons_(), father_(),
  distanceToFather_(0)
{
}

/** Assignation operator: *****************************************************/

AwareNode& AwareNode::operator=(const AwareNode & node)
{
  id_               = node.id_;
  father_           = 0;
  distanceToFather_ = node.distanceToFather_;
  sons_.clear();

  return * this;
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
  for(size_t i = 0; i < sons_.size(); i++)
  {
    if(sons_[i] == son) return i;
  }
  throw Exception("Son not found: " + TextTools::toString(son->getId()));
}

/******************************************************************************/

