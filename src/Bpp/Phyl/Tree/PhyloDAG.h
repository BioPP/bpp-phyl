// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_TREE_PHYLODAG_H
#define BPP_PHYL_TREE_PHYLODAG_H

#include <Bpp/Graph/AssociationDAGraphImplObserver.h>

#include "PhyloBranch.h"
#include "PhyloNode.h"

namespace bpp
{
/**
 * @brief Defines a Phylogenetic DAG based on a DAGraph & its
 * associationObserver.
 *
 */

class ParametrizablePhyloDAG;

class PhyloDAG :
  public AssociationDAGlobalGraphObserver<PhyloNode, PhyloBranch>
{
private:
  std::string name_;

public:
  PhyloDAG();

  PhyloDAG(const PhyloDAG* dag);

  PhyloDAG(const ParametrizablePhyloDAG& dag);

  template<class T, class U>
  PhyloDAG(AssociationDAGlobalGraphObserver<T, U> dag) :
    AssociationDAGlobalGraphObserver<PhyloNode, PhyloBranch>(dag),
    name_("")
  {}

  PhyloDAG* clone() const
  {
    return new PhyloDAG(*this);
  }

  /**
   * @brief DAG name.
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
   * name is not present in the PhyloDAG.
   *
   */
  std::shared_ptr<PhyloNode> getPhyloNode(const std::string& name) const;

  Vdouble getBranchLengths() const;

  /*
   *@brief set indexes to nodes & edges.
   * For each node, one (arbitrary) edge has the same index.
   */
  
  void resetNodesId();

  void setBranchLengths(double l);

  /**
   * @brief Multiply all branch lengths by a given factor.
   *
   * @param factor The factor to multiply all branch lengths with.
   */
  void scaleDAG(double factor);

  /**
   * @brief Prune a tree to a given set of leaf names
   *
   * @param leaves  the vector of leaf names restricting the tree
   */

  void pruneDAG(std::vector<std::string> leaves);

  /**
   * @brief Multiply all branch lengths under a Node by a given factor.
   *
   * @param node The node defining the subtree.
   * @param factor The factor to multiply all branch lengths with.
   */
  void scaleDAG(std::shared_ptr<PhyloNode> node, double factor);

  /**
   * @brief Add the lengths of branches of another phylotree to this
   * one. Just branch ids are considered, whatever the topology of
   * the trees.
   *
   * @param phylotree The added PhyloDAG
   */
  PhyloDAG& operator+=(const PhyloDAG& phylotree);

  /**
   * @brief Subtracts the lengths of branches of another phylotree to this
   * one. Just branch ids are considered, whatever the topology of
   * the trees.
   *
   * @param phylotree The added PhyloDAG
   */
  PhyloDAG& operator-=(const PhyloDAG& phylotree);

  /**
   * @brief Divides the lengths of branches of this phylotree by the
   * ones of another phylotree. Just branch ids are considered,
   * whatever the topology of the trees.
   *
   * @param phylotree The dividant PhyloDAG
   */
  PhyloDAG& operator/=(const PhyloDAG& phylotree);

  /**
   * @brief Multiplies the lengths of branches of this phylotree by the
   * ones of another phylotree. Just branch ids are considered,
   * whatever the topology of the trees.
   *
   * @param phylotree The dividant PhyloDAG
   */
  PhyloDAG& operator*=(const PhyloDAG& phylotree);

  // /**
  //  * @brief Concatenate the subtree under a Node (in a
  //  * DAGTemplate<Node>) to this PhyloTree, under the given phylonode.
  //  */
  // void addSubTree(std::shared_ptr<PhyloNode> phyloNode, const Node& node);
};
}

#endif // BPP_PHYL_TREE_PHYLOTREE_H
