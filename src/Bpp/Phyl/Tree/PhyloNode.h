// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_TREE_PHYLONODE_H
#define BPP_PHYL_TREE_PHYLONODE_H

#include <Bpp/Clonable.h>
#include <Bpp/Utils/MapTools.h>
#include <memory>

#include "PhyloTreeExceptions.h"

namespace bpp
{
/**
 * @brief For each PhyloNode of the tree, a description of the event at the node.
 *
 * It will determine the choice of the computing operator on this node.
 *
 */
class NodeEvent : public Clonable
{
  enum class NodeType
  {
    Speciation = 0,
    Mixture = 1,
    Hybridization = 2
  };

  NodeType nodeType_;

public:
  NodeEvent(NodeType type) : nodeType_(type) {}

  // NodeEvent(const NodeEvent& event) :
  //   nodeType_(event.nodeType_) {}

  NodeEvent* clone() const { return new NodeEvent(*this);}

  bool isSpeciation() const { return nodeType_ == NodeType::Speciation;}

  bool isMixture() const {return nodeType_ == NodeType::Mixture;}

  bool isHybridization() const { return nodeType_ == NodeType::Hybridization;}

public:
  static const NodeEvent speciationEvent;
  static const NodeEvent mixtureEvent;
  static const NodeEvent hybridizationEvent;

  std::string toString() const
  {
    if (isSpeciation())
      return "speciation";
    else if (isMixture())
      return "mixture";
    else
      return "hybridization";
  }
};


class PhyloNode
{
private:
  // a name, if specified
  std::string name_;

  // Node properties
  mutable std::map<std::string, Clonable*> properties_;

public:
  /**
   * @brief Build a new void Node object.
   */
  PhyloNode() :
    name_(""),
    properties_()
  {}

  /**
   * @brief Build a new Node with specified name.
   *
   *
   */
  PhyloNode(const std::string& name) :
    name_(name),
    properties_()
  {}

  /**
   * @brief Copy constructor.
   *
   * @param node The node to copy.
   */
  PhyloNode(const PhyloNode& node);

  /**
   * @brief Assignation operator.
   *
   * @param node the node to copy.
   * @return A reference toward this node.
   */

  PhyloNode& operator=(const PhyloNode& node);

  PhyloNode* clone() const { return new PhyloNode(*this); }

  /**
   * @brief destructor.
   *
   */
  virtual ~PhyloNode()
  {
    deleteProperties();
  }

public:
  /**
   * @name Name:
   *
   * @{
   */

  /**
   * @brief Get the name associated to this node, if there is one,
   * otherwise throw a NodeException.
   *
   * @return The name associated to this node.
   */
  std::string getName() const
  {
    if (!hasName()) throw PhyloNodePException("Node::getName: no name associated to this node.", this);
    return name_;
  }

  /**
   * @brief Give a name or update the name associated to the node.
   *
   * @param name The name to give to the node.
   */
  void setName(const std::string& name)
  {
    name_ = name;
  }

  /**
   * @brief Delete the name associated to this node (do nothing if there is no name).
   */
  void deleteName()
  {
    name_ = "";
  }

  /**
   * @brief Tell is this node has a name.
   *
   * @return True if name != 0.
   */
  bool hasName() const { return name_ != ""; }

  /** @} */

  /**
   * @name Node properties:
   *
   * @{
   */

  /**
   * @brief Set/add a node property.
   *
   * If no property with the same name is found, the new property will be added to the list.
   * Conversely, the property will be deleted and replaced by the new one.
   * If you want to keep a copy of the old property, consider using the removeProperty function before.
   *
   * @param name The name of the property to set.
   * @param property The property object (will be cloned).
   */
  void setProperty(const std::string& name, const Clonable& property)
  {
    if (hasProperty(name))
      delete properties_[name];
    properties_[name] = property.clone();
  }

  Clonable* getProperty(const std::string& name)
  {
    if (hasProperty(name))
      return properties_[name];
    else
      throw PhyloNodePropertyNotFoundException("", name, this);
  }

  const Clonable* getProperty(const std::string& name) const
  {
    if (hasProperty(name))
      return const_cast<const Clonable*>(properties_[name]);
    else
      throw PhyloNodePropertyNotFoundException("", name, this);
  }

  Clonable* removeProperty(const std::string& name)
  {
    if (hasProperty(name))
    {
      Clonable* removed = properties_[name];
      properties_.erase(name);
      return removed;
    }
    else
      throw PhyloNodePropertyNotFoundException("", name, this);
  }

  void deleteProperty(const std::string& name)
  {
    if (hasProperty(name))
    {
      delete properties_[name];
      properties_.erase(name);
    }
    else
      throw PhyloNodePropertyNotFoundException("", name, this);
  }

  /**
   * @brief Remove all node properties.
   *
   * Attached objects will not be deleted.
   */
  void removeProperties()
  {
    properties_.clear();
  }

  /**
   * @brief Delete all node properties.
   */
  void deleteProperties()
  {
    for (std::map<std::string, Clonable*>::iterator i = properties_.begin(); i != properties_.end(); i++)
    {
      delete i->second;
    }
    properties_.clear();
  }

  bool hasProperty(const std::string& name) const { return properties_.find(name) != properties_.end(); }

  std::vector<std::string> getPropertyNames() const { return MapTools::getKeys(properties_); }
  /** @} */
}; // end of class node
} // end of namespace bpp.

#else
namespace bpp { class PhyloNode; }
#endif // BPP_PHYL_TREE_PHYLONODE_H
