// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_ABSTRACTTREELIKELIHOODDATA_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_ABSTRACTTREELIKELIHOODDATA_H


#include "TreeLikelihoodData.h"

// From the STL:
#include <vector>
#include <map>

namespace bpp
{
/**
 * @brief Partial implementation of the TreeLikelihoodData interface.
 *
 * This data structure provides a simple compression, by performing and storing computations
 * only one time per identical sites.
 *
 * The compression is achieved by the TreeLikelihood object.
 * The correspondence between sites in the dataset and the arrays in the structures is given
 * by the rootPatternLinks_ array: the array indice for site @f$i@f$ if given by:
 * @code
 * rootPatternLinks_[i]
 * @endcode
 *
 * Finally, the rootWeights_ array gives for each array position, the number of sites with this
 * pattern.
 * The global likelihood is then given by the product of all likelihoods for each array position,
 * weighted by the corresponding number of sites.
 */
class AbstractTreeLikelihoodData :
  public TreeLikelihoodData
{
protected:
  /**
   * @brief Links between sites and patterns.
   *
   * The size of this vector is equal to the number of sites in the container,
   * each element corresponds to a site in the container and points to the
   * corresponding column in the likelihood array of the root node.
   * If the container contains no repeated site, there will be a strict
   * equivalence between each site and the likelihood array of the root node.
   * However, if this is not the case, some pointers may point toward the same
   * element in the likelihood array.
   */
  std::vector<size_t> rootPatternLinks_;

  /**
   * @brief The frequency of each site.
   */
  std::vector<unsigned int> rootWeights_;

  std::shared_ptr<const TreeTemplate<Node>> tree_;

  std::shared_ptr<const Alphabet> alphabet_;

public:
  AbstractTreeLikelihoodData(std::shared_ptr< const TreeTemplate<Node>> tree) :
    rootPatternLinks_(), rootWeights_(), tree_(tree), alphabet_(0) {}

  AbstractTreeLikelihoodData(const AbstractTreeLikelihoodData& atd) :
    rootPatternLinks_(atd.rootPatternLinks_),
    rootWeights_(atd.rootWeights_),
    tree_(atd.tree_),
    alphabet_(atd.alphabet_)
  {}

  AbstractTreeLikelihoodData& operator=(const AbstractTreeLikelihoodData& atd)
  {
    rootPatternLinks_ = atd.rootPatternLinks_;
    rootWeights_      = atd.rootWeights_;
    tree_             = atd.tree_;
    alphabet_         = atd.alphabet_;
    return *this;
  }


  virtual ~AbstractTreeLikelihoodData() {}

public:
  std::vector<size_t>& getRootArrayPositions() { return rootPatternLinks_; }
  const std::vector<size_t>& getRootArrayPositions() const { return rootPatternLinks_; }
  size_t getRootArrayPosition(const size_t site) const
  {
    return rootPatternLinks_[site];
  }
  unsigned int getWeight(size_t pos) const
  {
    return rootWeights_[pos];
  }
  const std::vector<unsigned int>& getWeights() const
  {
    return rootWeights_;
  }

  std::shared_ptr<const Alphabet> getAlphabet() const { return alphabet_; }

  std::shared_ptr< const TreeTemplate<Node>> getTree() const { return tree_; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_ABSTRACTTREELIKELIHOODDATA_H
