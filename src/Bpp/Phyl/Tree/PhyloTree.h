// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_TREE_PHYLOTREE_H
#define BPP_PHYL_TREE_PHYLOTREE_H

#include <Bpp/Graph/AssociationTreeGraphImplObserver.h>

#include "PhyloBranch.h"
#include "PhyloNode.h"
#include "TreeTemplate.h"

namespace bpp
{
/**
 * @brief Defines a Phylogenetic Tree based on a TreeGraph & its
 * associationObserver.
 *
 * @author Thomas Bigot
 */
class ParametrizablePhyloTree;

class PhyloTree :
  public AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranch>
{
private:
  std::string name_;

public:
  PhyloTree(bool rooted = false);

  PhyloTree(const PhyloTree* tree);

  PhyloTree(const ParametrizablePhyloTree& tree);

  template<class T, class U>
  PhyloTree(AssociationTreeGlobalGraphObserver<T, U> tree) :
    AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranch>(tree),
    name_("")
  {}

  PhyloTree* clone() const
  {
    return new PhyloTree(*this);
  }

  /**
   * @brief Tree name.
   *
   * @{
   */
  std::string getName() const
  {
    return name_;
  }

  void setName(const std::string& name)
  {
    name_ = name;
  }

  /** @} */

  std::vector<std::string> getAllLeavesNames() const;

  /*
   *@brief Get PhyloNode with given name, or null shared_ptr if the
   * name is not present in the PhyloTree.
   *
   */
  std::shared_ptr<PhyloNode> getPhyloNode(const std::string& name) const;

  Vdouble getBranchLengths() const;

  void resetNodesId();

  void setBranchLengths(double l);

  /**
   * @brief Multiply all branch lengths by a given factor.
   *
   * @param factor The factor to multiply all branch lengths with.
   */
  void scaleTree(double factor);

  /**
   * @brief Prune a tree to a given set of leaf names
   *
   * @param leaves  the vector of leaf names restricting the tree
   */

  void pruneTree(std::vector<std::string> leaves);

  /**
   * @brief Multiply all branch lengths under a Node by a given factor.
   *
   * @param node The node defining the subtree.
   * @param factor The factor to multiply all branch lengths with.
   */
  void scaleTree(std::shared_ptr<PhyloNode> node, double factor);

  /**
   * @brief Add the lengths of branches of another phylotree to this
   * one. Just branch ids are considered, whatever the topology of
   * the trees.
   *
   * @param phylotree The added PhyloTree
   */
  PhyloTree& operator+=(const PhyloTree& phylotree);

  /**
   * @brief Substracts the lengths of branches of another phylotree to this
   * one. Just branch ids are considered, whatever the topology of
   * the trees.
   *
   * @param phylotree The added PhyloTree
   */
  PhyloTree& operator-=(const PhyloTree& phylotree);

  /**
   * @brief Divides the lengths of branches of this phylotree by the
   * ones of another phylotree. Just branch ids are considered,
   * whatever the topology of the trees.
   *
   * @param phylotree The dividant PhyloTree
   */
  PhyloTree& operator/=(const PhyloTree& phylotree);

  /**
   * @brief Multiplies the lengths of branches of this phylotree by the
   * ones of another phylotree. Just branch ids are considered,
   * whatever the topology of the trees.
   *
   * @param phylotree The dividant PhyloTree
   */
  PhyloTree& operator*=(const PhyloTree& phylotree);

  /**
   * @brief Concatenate the subtree under a Node (in a
   * TreeTemplate<Node>) to this PhyloTree, under the given phylonode.
   */
  void addSubTree(std::shared_ptr<PhyloNode> phyloNode, const Node& node);
};
}

#endif // BPP_PHYL_TREE_PHYLOTREE_H
