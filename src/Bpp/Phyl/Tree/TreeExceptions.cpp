// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Text/TextTools.h>

#include "Node.h"
#include "Tree.h"
#include "TreeExceptions.h"

using namespace bpp;

/******************************************************************************/

NodePException::NodePException(const std::string& text, const Node* node) :
  NodeException(text, node->getId()), node_(node)
{}

/******************************************************************************/

NodeNotFoundException::NodeNotFoundException(const std::string& text, const std::string& id) :
  Exception("NodeNotFoundException: " + text + "(" + id + ")"),
  id_(id) {}

NodeNotFoundException::NodeNotFoundException(const std::string& text, int id) :
  Exception("NodeNotFoundException: " + text + "(" + TextTools::toString(id) + ")"),
  id_(TextTools::toString(id)) {}

/******************************************************************************/

TreeException::TreeException(const std::string& text, const Tree* tree) :
  Exception("TreeException: " + text + (tree != 0 ? "(" + tree->getName() + ")" : "")),
  tree_(tree) {}

/******************************************************************************/

UnrootedTreeException::UnrootedTreeException(const std::string& text, const Tree* tree) :
  TreeException("UnrootedTreeException: " + text, tree) {}

/******************************************************************************/
