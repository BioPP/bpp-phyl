// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_PARSIMONY_ABSTRACTTREEPARSIMONYDATA_H
#define BPP_PHYL_PARSIMONY_ABSTRACTTREEPARSIMONYDATA_H


#include "TreeParsimonyData.h"

namespace bpp
{
/**
 * @brief Partial implementation of the TreeParsimonyData interface.
 *
 * This data structure provides a simple compression, by performing and storing computations
 * only one time per identical sites.
 *
 * The compression is achieved by the TreeParsimonyScore object.
 * The correspondance between sites in the dataset and the arrays in the structures is given
 * by the rootPatternLinks_ array: the array indice for site @f$i@f$ if given by:
 * @code
 * rootPatternLinks_[i]
 * @endcode
 *
 * Finally, the rootWeights_ array gives for each array position, the number of sites with this
 * pattern.
 * The global parsimony score is then given by the sum of all scores for each array position,
 * weighted by the corresponding number of sites.
 */
class AbstractTreeParsimonyData :
  public virtual TreeParsimonyDataInterface
{
protected:
  std::vector<size_t> rootPatternLinks_;
  std::vector<unsigned int> rootWeights_;
  std::shared_ptr<const TreeTemplate<Node>> tree_;

public:
  AbstractTreeParsimonyData(std::shared_ptr<const TreeTemplate<Node>> tree) :
    rootPatternLinks_(),
    rootWeights_(),
    tree_(tree)
  {}

  AbstractTreeParsimonyData(const AbstractTreeParsimonyData& atpd) :
    rootPatternLinks_(atpd.rootPatternLinks_),
    rootWeights_(atpd.rootWeights_),
    tree_(atpd.tree_)
  {}

  AbstractTreeParsimonyData& operator=(const AbstractTreeParsimonyData& atpd)
  {
    rootPatternLinks_ = atpd.rootPatternLinks_;
    rootWeights_      = atpd.rootWeights_;
    tree_             = atpd.tree_;
    return *this;
  }


  virtual ~AbstractTreeParsimonyData() {}

public:
  size_t getRootArrayPosition(size_t site) const override
  {
    return rootPatternLinks_[site];
  }

  unsigned int getWeight(size_t pos) const
  {
    return rootWeights_[pos];
  }

  const TreeTemplate<Node>& tree() const override { return *tree_; }
  
  std::shared_ptr<const TreeTemplate<Node>> getTree() const override { return tree_; }

protected:
  void setTree(std::shared_ptr<const TreeTemplate<Node>> tree) { tree_ = tree; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_PARSIMONY_ABSTRACTTREEPARSIMONYDATA_H
