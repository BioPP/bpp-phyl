// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_TREE_TREEEXCEPTIONS_H
#define BPP_PHYL_TREE_TREEEXCEPTIONS_H

#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>


// From the STL:
#include <string>

namespace bpp
{
class Node;
class Tree;

/**
 * @brief General exception thrown when something is wrong with a particular node.
 */
class NodeException :
  public Exception
{
protected:
  int nodeId_;

public:
  /**
   * @brief Build a new NodePException.
   * @param text A message to be passed to the exception hierarchy.
   * @param nodeId The id of the node that threw the exception.
   */
  NodeException(const std::string& text, int nodeId) :
    Exception("NodeException: " + text + "(id:" + TextTools::toString(nodeId) + ")"),
    nodeId_(nodeId) {}

  virtual ~NodeException() {}

public:
  /**
   * @brief Get the id of node that threw the exception.
   *
   * @return The id of the faulty node.
   */
  virtual int getNodeId() const { return nodeId_; }
};


/**
 * @brief General exception thrown when something is wrong with a particular node.
 */
class NodePException :
  public NodeException
{
private:
  const Node* node_;

public:
  /**
   * @brief Build a new NodePException.
   * @param text A message to be passed to the exception hierarchy.
   * @param node A const pointer toward the node that threw the exception.
   */
  NodePException(const std::string& text, const Node* node = 0);

  /**
   * @brief Build a new NodePException.
   * @param text A message to be passed to the exception hierarchy.
   * @param nodeId The id of the node that threw the exception.
   */
  NodePException(const std::string& text, int nodeId) :
    NodeException(text, nodeId), node_(0) {}

  NodePException(const NodePException& nex) :
    NodeException(nex),
    node_(nex.node_)
  {}

  NodePException& operator=(const NodePException& nex)
  {
    NodeException::operator=(nex);
    node_ = nex.node_;
    return *this;
  }

  virtual ~NodePException() {}

public:
  /**
   * @brief Get the node that threw the exception.
   *
   * @return A pointer toward the faulty node.
   */
  virtual const Node* getNode() const { return node_; }
  /**
   * @brief Get the id of node that threw the exception.
   *
   * @return The id of the faulty node.
   */
  virtual int getNodeId() const { return nodeId_; }
};

/**
 * @brief General exception thrown if a property could not be found.
 */
class PropertyNotFoundException :
  public NodePException
{
private:
  std::string propertyName_;

public:
  /**
   * @brief Build a new PropertyNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param propertyName The name of the property.
   * @param node A const pointer toward the node that threw the exception.
   */
  PropertyNotFoundException(const std::string& text, const std::string& propertyName, const Node* node = 0) :
    NodePException("Property not found: " + propertyName + ". " + text, node),
    propertyName_(propertyName) {}

  /**
   * @brief Build a new PropertyNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param propertyName The name of the property.
   * @param nodeId The id of the node that threw the exception.
   */
  PropertyNotFoundException(const std::string& text, const std::string& propertyName, int nodeId) :
    NodePException("Property not found: " + propertyName + ". " + text, nodeId),
    propertyName_(propertyName) {}

  virtual ~PropertyNotFoundException() {}

public:
  /**
   * @brief Get the name of the property that could not be found.
   *
   * @return The name of the missing property.
   */
  virtual const std::string& getPropertyName() const { return propertyName_; }
};

/**
 * @brief Exception thrown when something is wrong with a particular node.
 */
class NodeNotFoundException :
  public Exception
{
private:
  std::string id_;

public:
  /**
   * @brief Build a new NodeNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param id   A string describing the node.
   */
  NodeNotFoundException(const std::string& text, const std::string& id);

  /**
   * @brief Build a new NodeNotFoundException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param id   A node identifier.
   */
  NodeNotFoundException(const std::string& text, int id);

  virtual ~NodeNotFoundException() {}

public:
  /**
   * @brief Get the node id that threw the exception.
   *
   * @return The id of the node.
   */
  virtual std::string getId() const { return id_; }
};

/**
 * @brief General exception thrown when something wrong happened in a tree.
 */
class TreeException :
  public Exception
{
private:
  const Tree* tree_;

public:
  /**
   * @brief Build a new TreeException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param tree A const pointer toward the tree that threw the exception.
   */
  TreeException(const std::string& text, const Tree* tree = 0);

  TreeException(const TreeException& tex) :
    Exception(tex),
    tree_(tex.tree_)
  {}

  TreeException& operator=(const TreeException& tex)
  {
    Exception::operator=(tex);
    tree_ = tex.tree_;
    return *this;
  }

  virtual ~TreeException() {}

public:
  /**
   * @brief Get the tree that threw the exception.
   *
   * @return The faulty tree
   */
  virtual const Tree* getTree() const { return tree_; }
};

/**
 * @brief Exception thrown when a tree is expected to be rooted.
 */
class UnrootedTreeException :
  public TreeException
{
public:
  /**
   * @brief Build a new UnrootedTreeException.
   *
   * @param text A message to be passed to the exception hierarchy.
   * @param tree A const pointer toward the tree that threw the exception.
   */
  UnrootedTreeException(const std::string& text, const Tree* tree = 0);

  virtual ~UnrootedTreeException() {}
};
} // end of namespace bpp.
#endif // BPP_PHYL_TREE_TREEEXCEPTIONS_H
