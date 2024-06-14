// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_DRASRTREELIKELIHOODDATA_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_DRASRTREELIKELIHOODDATA_H

#include <Bpp/Numeric/VectorTools.h>

#include "../../Model/SubstitutionModel.h"
#include "../../SitePatterns.h"
#include "AbstractTreeLikelihoodData.h"

// From the STL:
#include <map>

namespace bpp
{
/**
 * @brief Likelihood data structure for a node.
 *
 * This class is for use with the DRASRTreeParsimonyData class.
 *
 * Store all conditional likelihoods:
 * <pre>
 * x[i][c][s]
 *   |---------> Site i
 *      |------> Rate class c
 *         |---> Ancestral state s
 * </pre>
 * We call this the <i>likelihood array</i> for each node.
 * In the same way, we store first and second order derivatives.
 *
 * @see DRASRTreeLikelihoodData
 */
class DRASRTreeLikelihoodNodeData :
  public virtual TreeLikelihoodNodeData
{
private:
  mutable VVVdouble nodeLikelihoods_;
  mutable VVVdouble nodeDLikelihoods_;
  mutable VVVdouble nodeD2Likelihoods_;
  const Node* node_;

public:
  DRASRTreeLikelihoodNodeData() : nodeLikelihoods_(), nodeDLikelihoods_(), nodeD2Likelihoods_(), node_(0) {}

  DRASRTreeLikelihoodNodeData(const DRASRTreeLikelihoodNodeData& data) :
    nodeLikelihoods_(data.nodeLikelihoods_),
    nodeDLikelihoods_(data.nodeDLikelihoods_),
    nodeD2Likelihoods_(data.nodeD2Likelihoods_),
    node_(data.node_)
  {}

  DRASRTreeLikelihoodNodeData& operator=(const DRASRTreeLikelihoodNodeData& data)
  {
    nodeLikelihoods_   = data.nodeLikelihoods_;
    nodeDLikelihoods_  = data.nodeDLikelihoods_;
    nodeD2Likelihoods_ = data.nodeD2Likelihoods_;
    node_              = data.node_;
    return *this;
  }

  DRASRTreeLikelihoodNodeData* clone() const
  {
    return new DRASRTreeLikelihoodNodeData(*this);
  }

public:
  const Node* getNode() const { return node_; }
  void setNode(const Node* node) { node_ = node; }

  VVVdouble& getLikelihoodArray() { return nodeLikelihoods_; }
  const VVVdouble& getLikelihoodArray() const { return nodeLikelihoods_; }

  VVVdouble& getDLikelihoodArray() { return nodeDLikelihoods_; }
  const VVVdouble& getDLikelihoodArray() const { return nodeDLikelihoods_; }

  VVVdouble& getD2LikelihoodArray() { return nodeD2Likelihoods_; }
  const VVVdouble& getD2LikelihoodArray() const { return nodeD2Likelihoods_; }
};

/**
 * @brief discrete Rate Across Sites, (simple) Recursive likelihood data structure.
 */
class DRASRTreeLikelihoodData :
  public virtual AbstractTreeLikelihoodData
{
private:
  /**
   * @brief This contains all likelihood values used for computation.
   */
  mutable std::map<int, DRASRTreeLikelihoodNodeData> nodeData_;

  /**
   * @brief This map defines the pattern network.
   *
   * Let n1 be the id of a node in the tree, and n11 and n12 the ids of its sons.
   * Providing the likelihood array is known for nodes n11 and n12,
   * the likelihood array for node n1 and site <i>i</i> (_likelihood[n1][i]) must be computed
   * using arrays patternLinks_[n1][n11][i] and patternLinks_[n1][n12][i].
   * This network is initialized once for all in the constructor of this class.
   *
   * The double map contains the position of the site to use (second dimension)
   * of the likelihoods array.
   */
  mutable std::map<int, std::map<int, std::vector<size_t>>> patternLinks_;
  std::shared_ptr<AlignmentDataInterface> shrunkData_;
  size_t nbSites_;
  size_t nbStates_;
  size_t nbClasses_;
  size_t nbDistinctSites_;
  bool usePatterns_;

public:
  DRASRTreeLikelihoodData(std::shared_ptr< const TreeTemplate<Node>> tree, size_t nbClasses, bool usePatterns = true) :
    AbstractTreeLikelihoodData(tree),
    nodeData_(), patternLinks_(), shrunkData_(), nbSites_(0), nbStates_(0),
    nbClasses_(nbClasses), nbDistinctSites_(0), usePatterns_(usePatterns)
  {}

  DRASRTreeLikelihoodData(const DRASRTreeLikelihoodData& data) :
    AbstractTreeLikelihoodData(data),
    nodeData_(data.nodeData_),
    patternLinks_(data.patternLinks_),
    shrunkData_(),
    nbSites_(data.nbSites_), nbStates_(data.nbStates_),
    nbClasses_(data.nbClasses_), nbDistinctSites_(data.nbDistinctSites_),
    usePatterns_(data.usePatterns_)
  {
    if (data.shrunkData_)
      shrunkData_ = std::unique_ptr<AlignmentDataInterface>(data.shrunkData_->clone());
  }

  DRASRTreeLikelihoodData& operator=(const DRASRTreeLikelihoodData& data)
  {
    AbstractTreeLikelihoodData::operator=(data);
    nodeData_          = data.nodeData_;
    patternLinks_      = data.patternLinks_;
    nbSites_           = data.nbSites_;
    nbStates_          = data.nbStates_;
    nbClasses_         = data.nbClasses_;
    nbDistinctSites_   = data.nbDistinctSites_;
    if (data.shrunkData_)
      shrunkData_ = std::unique_ptr<AlignmentDataInterface>(data.shrunkData_->clone());
    else
      shrunkData_      = 0;
    usePatterns_       = data.usePatterns_;
    return *this;
  }

  virtual ~DRASRTreeLikelihoodData() {}

  DRASRTreeLikelihoodData* clone() const { return new DRASRTreeLikelihoodData(*this); }

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
  void setTree(std::shared_ptr< const TreeTemplate<Node>> tree)
  {
    tree_ = tree;
    for (auto& it : nodeData_)
    {
      int id = it.second.getNode()->getId();
      it.second.setNode(tree_->getNode(id));
    }
  }

  DRASRTreeLikelihoodNodeData& getNodeData(int nodeId)
  {
    return nodeData_[nodeId];
  }
  const DRASRTreeLikelihoodNodeData& getNodeData(int nodeId) const
  {
    return nodeData_[nodeId];
  }
  size_t getArrayPosition(int parentId, int sonId, size_t currentPosition) const
  {
    return patternLinks_[parentId][sonId][currentPosition];
  }
  size_t getRootArrayPosition(size_t currentPosition) const
  {
    return rootPatternLinks_[currentPosition];
  }
  const std::vector<size_t>& getArrayPositions(int parentId, int sonId) const
  {
    return patternLinks_[parentId][sonId];
  }
  std::vector<size_t>& getArrayPositions(int parentId, int sonId)
  {
    return patternLinks_[parentId][sonId];
  }
  size_t getArrayPosition(int parentId, int sonId, size_t currentPosition)
  {
    return patternLinks_[parentId][sonId][currentPosition];
  }

  VVVdouble& getLikelihoodArray(int nodeId)
  {
    return nodeData_[nodeId].getLikelihoodArray();
  }

  VVVdouble& getDLikelihoodArray(int nodeId)
  {
    return nodeData_[nodeId].getDLikelihoodArray();
  }

  VVVdouble& getD2LikelihoodArray(int nodeId)
  {
    return nodeData_[nodeId].getD2LikelihoodArray();
  }

  size_t getNumberOfDistinctSites() const { return nbDistinctSites_; }
  size_t getNumberOfSites() const { return nbSites_; }
  size_t getNumberOfStates() const { return nbStates_; }
  size_t getNumberOfClasses() const { return nbClasses_; }

  void initLikelihoods(
      const AlignmentDataInterface& sites,
      const TransitionModelInterface& model);

protected:
  /**
   * @brief This method initializes the leaves according to a sequence file.
   * likelihood is set to 1 for the state corresponding to the sequence site,
   * otherwise it is set to 0.
   *
   * All likelihood arrays at each nodes are initialized according to alphabet
   * size and sequences length, and filled with 1.
   *
   * NB: This method is recursive.
   *
   * @param node      The node defining the subtree to analyse.
   * @param sequences The data to be used for initialization.
   * @param model     The model to use.
   */
  virtual void initLikelihoods(
      const Node* node,
      const AlignmentDataInterface& sequences,
      const TransitionModelInterface& model);

  /**
   * @brief This method initializes the leaves according to a sequence file.
   *
   * likelihood is set to 1 for the state corresponding to the sequence site,
   * otherwise it is set to 0.
   *
   * All likelihood arrays at each nodes are initialized according to alphabet
   * size and sequences length, and filled with 1.
   *
   * NB: This method is recursive.
   *
   * @param node      The node defining the subtree to analyse.
   * @param sequences The data to be used for initialization.
   * @param model     The model to use.
   * @return The shrunk sub-dataset + indices for the subtree defined by <i>node</i>.
   */
  virtual std::unique_ptr<SitePatterns> initLikelihoodsWithPatterns(
      const Node* node,
      const AlignmentDataInterface& sequences,
      const TransitionModelInterface& model);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_DRASRTREELIKELIHOODDATA_H
