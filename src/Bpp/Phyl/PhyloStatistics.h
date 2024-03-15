// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_PHYLOSTATISTICS_H
#define BPP_PHYL_PHYLOSTATISTICS_H

#include <Bpp/Clonable.h>

#include "Tree/Tree.h"

// From the STL:
#include <vector>

namespace bpp
{
/**
 *  @brief Compute several quantities on a tree simulateously, optimizing the recursions on the tree.
 *
 *  This class uses a TreeTemplate. If the input tree is not a TreeTemplate, then a copy is performed
 *  before any computation.
 *
 *  @see TreeTools, TreeTemplateTools.
 */
class PhyloStatistics :
  public virtual Clonable
{
private:
  size_t numberOfLeaves_;
  size_t numberOfAncestors_;
  std::vector<double> branchLengths_;
  std::vector<double> nodeHeights_;
  std::vector<size_t> nodeDepths_;
  std::vector<size_t> nodeNumberOfSons_;
  std::vector<int> nodeIds_;

public:
  PhyloStatistics() :
    numberOfLeaves_(0), numberOfAncestors_(0),
    branchLengths_(), nodeHeights_(), nodeDepths_(), nodeNumberOfSons_(), nodeIds_()
  {}
  virtual ~PhyloStatistics() {}

  PhyloStatistics* clone() const { return new PhyloStatistics(*this); }

  /**
   * @brief Compute statistics for a given input tree.
   *
   * @param tree The tree for which the statistics should be computed.
   */
  void setTree(const Tree& tree);

  size_t getNumberOfLeaves() const { return numberOfLeaves_; }
  size_t getNumberOfAncestors() const { return numberOfAncestors_; }
  const std::vector<double>& getBranchLengths() const { return branchLengths_; }
  const std::vector<double>& getNodeHeights() const { return nodeHeights_; }
  const std::vector<size_t>& getNodeDepths() const { return nodeDepths_; }
  const std::vector<size_t>& getNodeNumberOfSons() const { return nodeNumberOfSons_; }
  const std::vector<int>& getNodeIds() const { return nodeIds_; }

private:
  void computeForSubtree_(const Node* node, double& height, size_t& depth);
};
} // end of namespace bpp.
#endif // BPP_PHYL_PHYLOSTATISTICS_H
