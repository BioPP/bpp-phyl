// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_PARSIMONY_DRTREEPARSIMONYSCORE_H
#define BPP_PHYL_PARSIMONY_DRTREEPARSIMONYSCORE_H


#include "../Tree/NNISearchable.h"
#include "../Tree/TreeTools.h"
#include "AbstractTreeParsimonyScore.h"
#include "DRTreeParsimonyData.h"

namespace bpp
{
/**
 * @brief Double recursive implementation of interface TreeParsimonyScore.
 *
 * Uses a DRTreeParsimonyData object for data storage.
 */
class DRTreeParsimonyScore :
  public AbstractTreeParsimonyScore,
  public virtual NNISearchable
{
private:
  std::unique_ptr<DRTreeParsimonyData> parsimonyData_;
  size_t nbDistinctSites_;

public:
  DRTreeParsimonyScore(
      std::shared_ptr<TreeTemplate<Node>> tree,
      std::shared_ptr<const SiteContainerInterface> data,
      bool verbose = true,
      bool includeGaps = false);

  DRTreeParsimonyScore(
      std::shared_ptr<TreeTemplate<Node>> tree,
      std::shared_ptr<const SiteContainerInterface> data,
      std::shared_ptr<const StateMapInterface> statesMap,
      bool verbose = true);

  DRTreeParsimonyScore(const DRTreeParsimonyScore& tp);

  DRTreeParsimonyScore& operator=(const DRTreeParsimonyScore& tp);

  virtual ~DRTreeParsimonyScore();

  DRTreeParsimonyScore* clone() const override { return new DRTreeParsimonyScore(*this); }

private:
  void init_(std::shared_ptr<const SiteContainerInterface> data, bool verbose);

protected:
  /**
   * @brief Compute all scores.
   *
   * Call the computeScoresPreorder and computeScoresPostorder methods, and then initialize rootBitsets_ and rootScores_.
   */
  virtual void computeScores();
  /**
   * @brief Compute scores (preorder algorithm).
   */
  virtual void computeScoresPreorder(const Node*);
  /**
   * @brief Compute scores (postorder algorithm).
   */
  virtual void computeScoresPostorder(const Node*);

public:
  unsigned int getScore() const override;
  unsigned int getScoreForSite(size_t site) const override;

  /**
   * @brief Compute bitsets and scores for each site for a node, in postorder.
   *
   * @param pData    The node data to use.
   * @param rBitsets The bitset array where to store the resulting bitsets.
   * @param rScores  The score array where to write the resulting scores.
   */
  static void computeScoresPostorderForNode(
      const DRTreeParsimonyNodeData& pData,
      std::vector<Bitset>& rBitsets,
      std::vector<unsigned int>& rScores);

  /**
   * @brief Compute bitsets and scores for each site for a node, in preorder.
   *
   * @param pData    The node data to use.
   * @param source   The node where we are coming from.
   * @param rBitsets The bitset array where to store the resulting bitsets.
   * @param rScores  The score array where to write the resulting scores.
   */
  static void computeScoresPreorderForNode(
      const DRTreeParsimonyNodeData& pData,
      const Node* source,
      std::vector<Bitset>& rBitsets,
      std::vector<unsigned int>& rScores);

  /**
   * @brief Compute bitsets and scores for each site for a node, in all directions.
   *
   * @param pData    The node data to use.
   * @param rBitsets The bitset array where to store the resulting bitsets.
   * @param rScores  The score array where to write the resulting scores.
   */
  static void computeScoresForNode(
      const DRTreeParsimonyNodeData& pData, std::vector<Bitset>& rBitsets,
      std::vector<unsigned int>& rScores);

  /**
   * @brief Compute bitsets and scores from an array of arrays.
   *
   * This method is the more general score computation.
   * Depending on what is passed as input, it may computes scores of a subtree
   * or the whole tree.
   *
   * @param iBitsets The vector of bitset arrays to use.
   * @param iScores  The vector of score arrays to use.
   * @param oBitsets The bitset array where to store the resulting bitsets.
   * @param oScores  The score array where to write the resulting scores.
   */
  static void computeScoresFromArrays(
      const std::vector<const std::vector<Bitset>*>& iBitsets,
      const std::vector<const std::vector<unsigned int>*>& iScores,
      std::vector<Bitset>& oBitsets,
      std::vector<unsigned int>& oScores);

  /**
   * @name Thee NNISearchable interface.
   *
   * @{
   */
  double getTopologyValue() const override { return getScore(); }

  double testNNI(int nodeId) const override;

  void doNNI(int nodeId) override;

  const Tree& topology() const override { return tree(); }

  void topologyChangeTested(const TopologyChangeEvent& event) override
  {
    parsimonyData_->reInit();
    computeScores();
  }

  void topologyChangeSuccessful(const TopologyChangeEvent& event) override {}
  /**@} */

  /**
   * @name Parsimony solution
   */
  static std::string PARSIMONY_SOLUTION_STATE;

  /**
   * @brief Sets the state of a node in a mapping 
   * @param node               The node to get the state of
   * @param state              The state that needs to be assigned to the node
   */
  void setNodeState(Node* node, size_t state);

  /**
   * @brief Extracts the state of a node in a mapping 
   * @param node              The node to get the state of
   * @return                  Node state is int
   */
  size_t getNodeState(const Node* node);

  /**
   * @brief Compute a maximum parsimony solution in DELTRAN manner.
   */  
  void computeSolution();
  /**@} */
};
} // end of namespace bpp.
#endif // BPP_PHYL_PARSIMONY_DRTREEPARSIMONYSCORE_H
