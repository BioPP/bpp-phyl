// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_PARSIMONY_DRTREEPARSIMONYDATA_H
#define BPP_PHYL_PARSIMONY_DRTREEPARSIMONYDATA_H


#include "../Model/StateMap.h"
#include "AbstractTreeParsimonyData.h"

// From SeqLib
#include <Bpp/Seq/Container/SiteContainer.h>

// From the STL:
#include <bitset>

namespace bpp
{
typedef std::bitset<21> Bitset; // 20AA + gaps, codon not lalowed so far :s

/**
 * @brief Parsimony data structure for a node.
 *
 * This class is for use with the DRTreeParsimonyData class.
 *
 * Store for each neighbor node
 * - a vector of bitsets,
 * - a vector of score for the corresponding subtree.
 *
 * @see DRTreeParsimonyData
 */
class DRTreeParsimonyNodeData :
  public TreeParsimonyNodeDataInterface
{
private:
  mutable std::map<int, std::vector<Bitset>> nodeBitsets_;
  mutable std::map<int, std::vector<unsigned int>> nodeScores_;
  const Node* node_;

public:
  DRTreeParsimonyNodeData() :
    nodeBitsets_(),
    nodeScores_(),
    node_(0)
  {}

  DRTreeParsimonyNodeData(const DRTreeParsimonyNodeData& tpnd) :
    nodeBitsets_(tpnd.nodeBitsets_),
    nodeScores_(tpnd.nodeScores_),
    node_(tpnd.node_)
  {}

  DRTreeParsimonyNodeData& operator=(const DRTreeParsimonyNodeData& tpnd)
  {
    nodeBitsets_ = tpnd.nodeBitsets_;
    nodeScores_  = tpnd.nodeScores_;
    node_        = tpnd.node_;
    return *this;
  }

  DRTreeParsimonyNodeData* clone() const { return new DRTreeParsimonyNodeData(*this); }

public:
  const Node* getNode() const { return node_; }

  void setNode(const Node* node) { node_ = node; }

  std::vector<Bitset>& getBitsetsArrayForNeighbor(int neighborId)
  {
    return nodeBitsets_[neighborId];
  }
  const std::vector<Bitset>& getBitsetsArrayForNeighbor(int neighborId) const
  {
    return nodeBitsets_[neighborId];
  }
  std::vector<unsigned int>& getScoresArrayForNeighbor(int neighborId)
  {
    return nodeScores_[neighborId];
  }
  const std::vector<unsigned int>& getScoresArrayForNeighbor(int neighborId) const
  {
    return nodeScores_[neighborId];
  }

  bool isNeighbor(int neighborId) const
  {
    return nodeBitsets_.find(neighborId) != nodeBitsets_.end();
  }

  void eraseNeighborArrays()
  {
    nodeBitsets_.erase(nodeBitsets_.begin(), nodeBitsets_.end());
    nodeScores_.erase(nodeScores_.begin(), nodeScores_.end());
  }
};

/**
 * @brief Parsimony data structure for a leaf.
 *
 * This class is for use with the DRTreeParsimonyData class.
 *
 * Store the vector of bitsets associated to a leaf.
 *
 * @see DRTreeParsimonyData
 */
class DRTreeParsimonyLeafData :
  public TreeParsimonyNodeDataInterface
{
private:
  mutable std::vector<Bitset> leafBitsets_;
  const Node* leaf_;

public:
  DRTreeParsimonyLeafData() :
    leafBitsets_(),
    leaf_(0)
  {}

  DRTreeParsimonyLeafData(const DRTreeParsimonyLeafData& tpld) :
    leafBitsets_(tpld.leafBitsets_),
    leaf_(tpld.leaf_)
  {}

  DRTreeParsimonyLeafData& operator=(const DRTreeParsimonyLeafData& tpld)
  {
    leafBitsets_ = tpld.leafBitsets_;
    leaf_        = tpld.leaf_;
    return *this;
  }


  DRTreeParsimonyLeafData* clone() const { return new DRTreeParsimonyLeafData(*this); }

public:
  const Node* getNode() const { return leaf_; }
  void setNode(const Node* node) { leaf_ = node; }

  std::vector<Bitset>& getBitsetsArray()
  {
    return leafBitsets_;
  }
  const std::vector<Bitset>& getBitsetsArray() const
  {
    return leafBitsets_;
  }
};

/**
 * @brief Parsimony data structure for double-recursive (DR) algorithm.
 *
 * States are coded using bitsets for faster computing (@see AbstractTreeParsimonyData).
 * For each inner node in the tree, we store a DRTreeParsimonyNodeData object in nodeData_.
 * For each leaf node in the tree, we store a DRTreeParsimonyLeafData object in leafData_.
 *
 * The dataset is first compressed, removing all identical sites.
 * The resulting dataset is stored in shrunkData_.
 * The corresponding positions are stored in rootPatternLinks_, inherited from AbstractTreeParsimonyData.
 */
class DRTreeParsimonyData :
  public AbstractTreeParsimonyData
{
private:
  mutable std::map<int, DRTreeParsimonyNodeData> nodeData_;
  mutable std::map<int, DRTreeParsimonyLeafData> leafData_;
  mutable std::vector<Bitset> rootBitsets_;
  mutable std::vector<unsigned int> rootScores_;
  std::unique_ptr<SiteContainerInterface> shrunkData_;
  size_t nbSites_;
  size_t nbStates_;
  size_t nbDistinctSites_;

public:
  DRTreeParsimonyData(std::shared_ptr<const TreeTemplate<Node>> tree) :
    AbstractTreeParsimonyData(tree),
    nodeData_(),
    leafData_(),
    rootBitsets_(),
    rootScores_(),
    shrunkData_(nullptr),
    nbSites_(0),
    nbStates_(0),
    nbDistinctSites_(0)
  {}

  DRTreeParsimonyData(const DRTreeParsimonyData& data);

  DRTreeParsimonyData& operator=(const DRTreeParsimonyData& data);

  virtual ~DRTreeParsimonyData() {}

  DRTreeParsimonyData* clone() const override { return new DRTreeParsimonyData(*this); }

public:
  /**
   * @brief Set the tree associated to the data.
   *
   * All node data will be actualized accordingly by calling the setNode() method on the corresponding nodes.
   * @warning: the old tree and the new tree must be two clones! And particularly, they have to share the
   * same topology and nodes id.
   *
   * @param tree The tree to be associated to this data.
   */
  void setTree(std::shared_ptr<const TreeTemplate<Node>> tree)
  {
    AbstractTreeParsimonyData::setTree(tree);
    for (auto& it : nodeData_)
    {
      int id = it.second.getNode()->getId();
      it.second.setNode(tree_->getNode(id));
    }
    for (auto& it : leafData_)
    {
      int id = it.second.getNode()->getId();
      it.second.setNode(tree_->getNode(id));
    }
  }

  DRTreeParsimonyNodeData& nodeData(int nodeId) override
  {
    return nodeData_[nodeId];
  }
  const DRTreeParsimonyNodeData& nodeData(int nodeId) const override
  {
    return nodeData_[nodeId];
  }

  DRTreeParsimonyLeafData& leafData(int nodeId)
  {
    return leafData_[nodeId];
  }
  const DRTreeParsimonyLeafData& leafData(int nodeId) const
  {
    return leafData_[nodeId];
  }

  std::vector<Bitset>& getBitsetsArray(int nodeId, int neighborId)
  {
    return nodeData_[nodeId].getBitsetsArrayForNeighbor(neighborId);
  }
  const std::vector<Bitset>& getBitsetsArray(int nodeId, int neighborId) const
  {
    return nodeData_[nodeId].getBitsetsArrayForNeighbor(neighborId);
  }

  std::vector<unsigned int>& getScoresArray(int nodeId, int neighborId)
  {
    return nodeData_[nodeId].getScoresArrayForNeighbor(neighborId);
  }
  const std::vector<unsigned int>& getScoresArray(int nodeId, int neighborId) const
  {
    return nodeData_[nodeId].getScoresArrayForNeighbor(neighborId);
  }

  size_t getArrayPosition(int parentId, int sonId, size_t currentPosition) const override
  {
    return currentPosition;
  }

  std::vector<Bitset>& getRootBitsets() { return rootBitsets_; }
  const std::vector<Bitset>& getRootBitsets() const { return rootBitsets_; }
  const Bitset& getRootBitset(size_t i) const { return rootBitsets_[i]; }

  std::vector<unsigned int>& getRootScores() { return rootScores_; }
  const std::vector<unsigned int>& getRootScores() const { return rootScores_; }
  unsigned int getRootScore(size_t i) const { return rootScores_[i]; }

  size_t getNumberOfDistinctSites() const { return nbDistinctSites_; }
  size_t getNumberOfSites() const { return nbSites_; }
  size_t getNumberOfStates() const { return nbStates_; }

  void init(std::shared_ptr<const SiteContainerInterface> sites,
      std::shared_ptr<const StateMapInterface> stateMap);
  void reInit();

protected:
  void init_(const Node* node,
  std::shared_ptr<const SiteContainerInterface> sites,
  std::shared_ptr<const StateMapInterface> stateMap);

  void reInit_(const Node* node);
};
} // end of namespace bpp.
#endif // BPP_PHYL_PARSIMONY_DRTREEPARSIMONYDATA_H
