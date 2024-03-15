// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_PARSIMONY_TREEPARSIMONYDATA_H
#define BPP_PHYL_PARSIMONY_TREEPARSIMONYDATA_H

#include <Bpp/Clonable.h>

#include "../Tree/Node.h"
#include "../Tree/TreeTemplate.h"

namespace bpp
{
/**
 * @brief TreeParsimonyScore node data structure.
 *
 * Stores inner computation for a given node.
 *
 * @see TreeParsimonyData
 */
class TreeParsimonyNodeDataInterface :
  public virtual Clonable
{
public:
  TreeParsimonyNodeDataInterface() {}
  virtual ~TreeParsimonyNodeDataInterface() {}

  TreeParsimonyNodeDataInterface* clone() const = 0;

public:
  /**
   * @brief Get the node associated to this data structure.
   *
   * @return The node associated to this structure.
   */
  virtual const Node* getNode() const = 0;

  /**
   * @brief Set the node associated to this data
   *
   * @param node A pointer toward the node to be associated to this data.
   */
  virtual void setNode(const Node* node) = 0;
};

/**
 * @brief TreeParsimonyScore data structure.
 *
 * Stores all the inner computations:
 * - subtree scores and ancestral states for each node,
 * - correspondance between sites in the dataset and array indices.
 *
 * @see TreeParsimonyNodeData
 */
class TreeParsimonyDataInterface :
  public virtual Clonable
{
public:
  TreeParsimonyDataInterface() {}
  virtual ~TreeParsimonyDataInterface() {}

  TreeParsimonyDataInterface* clone() const override = 0;

public:
  virtual const TreeTemplate<Node>& tree() const = 0;
  virtual std::shared_ptr<const TreeTemplate<Node>> getTree() const = 0;
  virtual size_t getArrayPosition(int parentId, int sonId, size_t currentPosition) const = 0;
  virtual size_t getRootArrayPosition(size_t site) const = 0;
  virtual TreeParsimonyNodeDataInterface& nodeData(int nodeId) = 0;
  virtual const TreeParsimonyNodeDataInterface& nodeData(int nodeId) const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_PARSIMONY_TREEPARSIMONYDATA_H
