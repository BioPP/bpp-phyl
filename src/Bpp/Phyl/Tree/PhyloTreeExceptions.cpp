// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Text/TextTools.h>

#include "PhyloBranch.h"
#include "PhyloNode.h"
#include "PhyloTree.h"
#include "PhyloTreeExceptions.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

PhyloNodePException::PhyloNodePException(const std::string& text, const PhyloTree& tree, const std::shared_ptr<PhyloNode> node) :
  PhyloNodeException(text, tree.getNodeIndex(node)), node_(node.get())
{}

/******************************************************************************/

PhyloNodePException::PhyloNodePException(const std::string& text, const PhyloNode* node) :
  PhyloNodeException(text, 0), node_(node)
{}

/******************************************************************************/

PhyloNodeNotFoundException::PhyloNodeNotFoundException(const std::string& text, const std::string& id) :
  Exception("NodeNotFoundException: " + text + "(" + id + ")"),
  id_(id) {}

PhyloNodeNotFoundException::PhyloNodeNotFoundException(const std::string& text, unsigned int id) :
  Exception("NodeNotFoundException: " + text + "(" + TextTools::toString(id) + ")"),
  id_(TextTools::toString(id)) {}

/******************************************************************************/

PhyloBranchPException::PhyloBranchPException(const std::string& text, const PhyloTree& tree, const std::shared_ptr<PhyloBranch> branch) :
  PhyloBranchException(text, tree.getEdgeIndex(branch)), branch_(branch.get())
{}

/******************************************************************************/

PhyloBranchPException::PhyloBranchPException(const std::string& text, const PhyloBranch* branch) :
  PhyloBranchException(text, 0), branch_(branch)
{}

/******************************************************************************/

PhyloBranchNotFoundException::PhyloBranchNotFoundException(const std::string& text, const std::string& id) :
  Exception("BranchNotFoundException: " + text + "(" + id + ")"),
  id_(id) {}

/******************************************************************************/

PhyloBranchNotFoundException::PhyloBranchNotFoundException(const std::string& text, unsigned int id) :
  Exception("BranchNotFoundException: " + text + "(" + TextTools::toString(id) + ")"),
  id_(TextTools::toString(id)) {}

/******************************************************************************/

PhyloTreeException::PhyloTreeException(const std::string& text, const PhyloTree* tree) :
  Exception("PhyloTreeException: " + text + (tree != 0 ? "(" + tree->getName() + ")" : "")),
  tree_(tree) {}

/******************************************************************************/

UnrootedPhyloTreeException::UnrootedPhyloTreeException(const std::string& text, const PhyloTree* tree) :
  PhyloTreeException("UnrootedPhyloTreeException: " + text, tree) {}

/******************************************************************************/
