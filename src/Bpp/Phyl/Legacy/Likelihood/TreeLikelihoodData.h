// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_TREELIKELIHOODDATA_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_TREELIKELIHOODDATA_H


#include "../../Tree/Node.h"
#include "../../Tree/TreeTemplate.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignmentData.h>

namespace bpp
{
/**
 * @brief TreeLikelihood partial data structure.
 *
 * Stores inner computation for a given node.
 *
 * @see TreeLikelihoodData
 */
class TreeLikelihoodNodeData :
  public virtual Clonable
{
public:
  TreeLikelihoodNodeData() {}
  virtual ~TreeLikelihoodNodeData() {}

  TreeLikelihoodNodeData* clone() const = 0;

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
   * A pointer toward this node will be created and associated to this data.
   *
   * @param node The node to be associated to this data.
   */
  virtual void setNode(const Node* node) = 0;
};

/**
 * @brief TreeLikelihood data structure.
 *
 * Stores all the inner computations:
 * - conditionnal likelihoods for each node,
 * - correspondance between sites in the dataset and array indices.
 *
 * @see TreeLikelihoodNodeData
 */
class TreeLikelihoodData :
  public virtual Clonable
{
public:
  TreeLikelihoodData() {}
  virtual ~TreeLikelihoodData() {}

  TreeLikelihoodData* clone() const = 0;

public:
  virtual std::shared_ptr<const Alphabet> getAlphabet() const = 0;
  virtual std::shared_ptr< const TreeTemplate<Node> > getTree() const = 0;
  virtual size_t getArrayPosition(int parentId, int sonId, size_t currentPosition) const = 0;
  virtual size_t getRootArrayPosition(size_t site) const = 0;
  virtual TreeLikelihoodNodeData& getNodeData(int nodeId) = 0;
  virtual const TreeLikelihoodNodeData& getNodeData(int nodeId) const = 0;

  /**
   * @return The number of non redundant patterns.
   */
  virtual size_t getNumberOfDistinctSites() const = 0;

  /**
   * @return The total number of sites.
   */
  virtual size_t getNumberOfSites() const = 0;

  /**
   * @return Get the number of states used in the model.
   */
  virtual size_t getNumberOfStates() const = 0;

  /**
   * @return The frequency of a given pattern.
   */
  virtual unsigned int getWeight(size_t pos) const = 0;

  /**
   * @return Frequencies for each pattern.
   */
  virtual const std::vector<unsigned int>& getWeights() const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_TREELIKELIHOODDATA_H
