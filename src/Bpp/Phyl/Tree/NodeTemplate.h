// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_TREE_NODETEMPLATE_H
#define BPP_PHYL_TREE_NODETEMPLATE_H


#include "Node.h"

namespace bpp
{
/**
 * @brief The NodeTemplate class.
 *
 * This class inherits from the Node class.
 * Its is a generic way to store any information to a node.
 * A NodeTemplate only add the setInfos and getInfos methods,
 * which set and retrieve a NodeInfo, whose type is given as a template
 * of the class.
 * This class is mainly for computation conveniency, one may define a NodeInfo
 * class with several results attached.
 * An example is provided in the PGMA class.
 * Another way is to use a map<Node *, NodeInfos>, with the limitation of the
 * map access.
 * One may also wish to use the property system of the Node class, but
 * properties are stored as in a map<string, Clonable *>, with the drawbacks
 * of the slow map access and the systematic use of dynamic_cast<NodeInfo *> to
 * convert from Clonable *.
 *
 * This class redefines all constructors and access methods (get*) with return
 * types as NodeTemplate and not Node (using covariant returns).
 *
 * @see Node, TreeTemplate
 */
template<class NodeInfos>
class NodeTemplate :
  public Node
{
  friend class TreeTemplateTools;

private:
  NodeInfos infos_;

public:
  /**
   * @brief Build a new void NodeTemplate object.
   */
  NodeTemplate() : Node(), infos_() {}

  /**
   * @brief Build a new NodeTemplate with specified id.
   */
  NodeTemplate(int id) : Node(id), infos_() {}

  /**
   * @brief Build a new NodeTemplate with specified name.
   */
  NodeTemplate(const std::string& name) : Node(name), infos_() {}

  /**
   * @brief Build a new NodeTemplate with specified id and name.
   */
  NodeTemplate(int id, const std::string& name) : Node(id, name), infos_() {}

protected:
  /**
   * @brief Copy constructor.
   *
   * @param node The node to copy.
   */
  NodeTemplate(const Node& node) : Node(node), infos_() {}

  /**
   * @brief Copy constructor.
   *
   * @param node The node to copy.
   */
  NodeTemplate(const NodeTemplate<NodeInfos>& node) :
    Node(node), infos_(node.infos_)
  {}

  /**
   * @brief Assignation operator.
   *
   * @param node the node to copy.
   * @return A reference toward this node.
   */
  NodeTemplate<NodeInfos>& operator=(const NodeTemplate<NodeInfos>& node)
  {
    Node::operator=(node);
    infos_ = node.infos_;
    return *this;
  }

  NodeTemplate<NodeInfos>* clone() const { return new NodeTemplate<NodeInfos>(*this); }

public:
  virtual ~NodeTemplate() {}

public:
  const NodeTemplate<NodeInfos>* getFather() const { return dynamic_cast<const NodeTemplate<NodeInfos>*>(father_); }

  NodeTemplate<NodeInfos>* getFather() { return dynamic_cast<NodeTemplate<NodeInfos>*>(father_); }

  NodeTemplate<NodeInfos>* removeFather() { NodeTemplate<NodeInfos>* f = dynamic_cast<NodeTemplate<NodeInfos>*>(father_); father_ = 0; return f; }

  const NodeTemplate<NodeInfos>* getSon(size_t i) const { return dynamic_cast<NodeTemplate<NodeInfos>*>(sons_[i]); }

  NodeTemplate<NodeInfos>* getSon(size_t i) { return dynamic_cast<NodeTemplate<NodeInfos>*>(sons_[i]); }

  std::vector<const NodeTemplate<NodeInfos>*> getNeighbors() const
  {
    std::vector<const Node*> neighbors = Node::getNeighbors();
    std::vector<const NodeTemplate<NodeInfos>*> neighbors2(neighbors.size());
    for (size_t i = 0; i < neighbors.size(); i++)
    {
      neighbors2[i] = dynamic_cast<const NodeTemplate<NodeInfos>*>(neighbors[i]);
    }
    return neighbors2;
  }

  std::vector<NodeTemplate<NodeInfos>*> getNeighbors()
  {
    std::vector<Node*> neighbors = Node::getNeighbors();
    std::vector<NodeTemplate<NodeInfos>*> neighbors2(neighbors.size());
    for (size_t i = 0; i < neighbors.size(); i++)
    {
      neighbors2[i] = dynamic_cast<NodeTemplate<NodeInfos>*>(neighbors[i]);
    }
    return neighbors2;
  }

  NodeTemplate<NodeInfos>* operator[](int i) { return dynamic_cast<NodeTemplate<NodeInfos>*>((i < 0) ? father_ : sons_[i]); }

  const NodeTemplate<NodeInfos>* operator[](int i) const { return dynamic_cast<const NodeTemplate<NodeInfos>*>((i < 0) ? father_ : sons_[i]); }


  // Specific methods:

  /**
   * @return A reference toward the information object associated to this node.
   */
  virtual const NodeInfos& getInfos() const { return infos_; }

  /**
   * @return A reference toward the information object associated to this node.
   */
  virtual NodeInfos& getInfos() { return infos_; }

  /**
   * @brief Set the information to be associated to this node.
   *
   * @param infos An information object.
   */
  virtual void setInfos(const NodeInfos& infos) { infos_ = infos; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_TREE_NODETEMPLATE_H
