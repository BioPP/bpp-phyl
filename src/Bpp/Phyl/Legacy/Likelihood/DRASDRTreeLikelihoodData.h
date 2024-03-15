// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_DRASDRTREELIKELIHOODDATA_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_DRASDRTREELIKELIHOODDATA_H


#include "../../Model/SubstitutionModel.h"
#include "../../SitePatterns.h"
#include "../../PatternTools.h"
#include "AbstractTreeLikelihoodData.h"

// From SeqLib:
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>

// From the STL:
#include <map>

namespace bpp
{
/**
 * @brief Likelihood data structure for a leaf.
 *
 * This class is for use with the DRASDRTreeLikelihoodData class.
 *
 * Store the likelihoods arrays associated to a leaf.
 *
 * @see DRASDRTreeLikelihoodData
 */
class DRASDRTreeLikelihoodLeafData :
  public virtual TreeLikelihoodNodeData
{
private:
  mutable VVdouble leafLikelihood_;
  const Node* leaf_;

public:
  DRASDRTreeLikelihoodLeafData() : leafLikelihood_(), leaf_(0) {}

  DRASDRTreeLikelihoodLeafData(const DRASDRTreeLikelihoodLeafData& data) :
    leafLikelihood_(data.leafLikelihood_), leaf_(data.leaf_) {}

  DRASDRTreeLikelihoodLeafData& operator=(const DRASDRTreeLikelihoodLeafData& data)
  {
    leafLikelihood_ = data.leafLikelihood_;
    leaf_           = data.leaf_;
    return *this;
  }

  DRASDRTreeLikelihoodLeafData* clone() const
  {
    return new DRASDRTreeLikelihoodLeafData(*this);
  }

public:
  const Node* getNode() const { return leaf_; }
  void setNode(const Node* node) { leaf_ = node; }

  VVdouble& getLikelihoodArray()  { return leafLikelihood_;  }
};

/**
 * @brief Likelihood data structure for a node.
 *
 * This class is for use with the DRASDRTreeLikelihoodData class.
 *
 * Store for each neighbor node an array with conditionnal likelihoods.
 *
 * @see DRASDRTreeLikelihoodData
 */
class DRASDRTreeLikelihoodNodeData :
  public virtual TreeLikelihoodNodeData
{
private:
  /**
   * @brief This contains all likelihood values used for computation.
   *
   * <pre>
   * x[b][i][c][s]
   *   |------------> Neighbor node of n (id)
   *     |---------> Site i
   *         |------> Rate class c
   *            |---> Ancestral state s
   * </pre>
   * We call this the <i>likelihood array</i> for each node.
   */

  mutable std::map<int, VVVdouble> nodeLikelihoods_;
  /**
   * @brief This contains all likelihood first order derivatives values used for computation.
   *
   * <pre>
   * x[i]
   *   |---------> Site i
   * </pre>
   * We call this the <i>dLikelihood array</i> for each node.
   */
  mutable Vdouble nodeDLikelihoods_;

  /**
   * @brief This contains all likelihood second order derivatives values used for computation.
   *
   * <pre>
   * x[i]
   |---------> Site i
   * </pre>
   * We call this the <i>d2Likelihood array</i> for each node.
   */
  mutable Vdouble nodeD2Likelihoods_;

  const Node* node_;

public:
  DRASDRTreeLikelihoodNodeData() : nodeLikelihoods_(), nodeDLikelihoods_(), nodeD2Likelihoods_(), node_(0) {}

  DRASDRTreeLikelihoodNodeData(const DRASDRTreeLikelihoodNodeData& data) :
    nodeLikelihoods_(data.nodeLikelihoods_),
    nodeDLikelihoods_(data.nodeDLikelihoods_),
    nodeD2Likelihoods_(data.nodeD2Likelihoods_),
    node_(data.node_)
  {}

  DRASDRTreeLikelihoodNodeData& operator=(const DRASDRTreeLikelihoodNodeData& data)
  {
    nodeLikelihoods_   = data.nodeLikelihoods_;
    nodeDLikelihoods_  = data.nodeDLikelihoods_;
    nodeD2Likelihoods_ = data.nodeD2Likelihoods_;
    node_              = data.node_;
    return *this;
  }

  virtual ~DRASDRTreeLikelihoodNodeData() {}

  DRASDRTreeLikelihoodNodeData* clone() const
  {
    return new DRASDRTreeLikelihoodNodeData(*this);
  }

public:
  const Node* getNode() const { return node_; }

  void setNode(const Node* node) { node_ = node; }

  const std::map<int, VVVdouble>& getLikelihoodArrays() const { return nodeLikelihoods_; }

  std::map<int, VVVdouble>& getLikelihoodArrays() { return nodeLikelihoods_; }

  VVVdouble& getLikelihoodArrayForNeighbor(int neighborId)
  {
    return nodeLikelihoods_[neighborId];
  }

  const VVVdouble& getLikelihoodArrayForNeighbor(int neighborId) const
  {
    return nodeLikelihoods_[neighborId];
  }

  Vdouble& getDLikelihoodArray() { return nodeDLikelihoods_;  }

  const Vdouble& getDLikelihoodArray() const {  return nodeDLikelihoods_;  }

  Vdouble& getD2LikelihoodArray()  {  return nodeD2Likelihoods_; }

  const Vdouble& getD2LikelihoodArrayForNeighbor() const { return nodeD2Likelihoods_; }

  bool isNeighbor(int neighborId) const
  {
    return nodeLikelihoods_.find(neighborId) != nodeLikelihoods_.end();
  }

  void eraseNeighborArrays()
  {
    nodeLikelihoods_.erase(nodeLikelihoods_.begin(), nodeLikelihoods_.end());
    nodeDLikelihoods_.erase(nodeDLikelihoods_.begin(), nodeDLikelihoods_.end());
    nodeD2Likelihoods_.erase(nodeD2Likelihoods_.begin(), nodeD2Likelihoods_.end());
  }
};

/**
 * @brief Likelihood data structure for rate across sites models, using a double-recursive algorithm.
 */
class DRASDRTreeLikelihoodData :
  public virtual AbstractTreeLikelihoodData
{
private:
  mutable std::map<int, DRASDRTreeLikelihoodNodeData> nodeData_;
  mutable std::map<int, DRASDRTreeLikelihoodLeafData> leafData_;
  mutable VVVdouble rootLikelihoods_;
  mutable VVdouble rootLikelihoodsS_;
  mutable Vdouble rootLikelihoodsSR_;

  std::shared_ptr<AlignmentDataInterface> shrunkData_;
  size_t nbSites_;
  size_t nbStates_;
  size_t nbClasses_;
  size_t nbDistinctSites_;

public:
  DRASDRTreeLikelihoodData(std::shared_ptr< const TreeTemplate<Node> > tree, size_t nbClasses) :
    AbstractTreeLikelihoodData(tree),
    nodeData_(), leafData_(), rootLikelihoods_(), rootLikelihoodsS_(), rootLikelihoodsSR_(),
    shrunkData_(nullptr), nbSites_(0), nbStates_(0), nbClasses_(nbClasses), nbDistinctSites_(0)
  {}

  DRASDRTreeLikelihoodData(const DRASDRTreeLikelihoodData& data) :
    AbstractTreeLikelihoodData(data),
    nodeData_(data.nodeData_), leafData_(data.leafData_),
    rootLikelihoods_(data.rootLikelihoods_),
    rootLikelihoodsS_(data.rootLikelihoodsS_),
    rootLikelihoodsSR_(data.rootLikelihoodsSR_),
    shrunkData_(nullptr),
    nbSites_(data.nbSites_), nbStates_(data.nbStates_),
    nbClasses_(data.nbClasses_), nbDistinctSites_(data.nbDistinctSites_)
  {
    if (data.shrunkData_)
      shrunkData_ = std::shared_ptr<AlignmentDataInterface>(data.shrunkData_->clone());
  }

  DRASDRTreeLikelihoodData& operator=(const DRASDRTreeLikelihoodData& data)
  {
    AbstractTreeLikelihoodData::operator=(data);
    nodeData_          = data.nodeData_;
    leafData_          = data.leafData_;
    rootLikelihoods_   = data.rootLikelihoods_;
    rootLikelihoodsS_  = data.rootLikelihoodsS_;
    rootLikelihoodsSR_ = data.rootLikelihoodsSR_;
    nbSites_           = data.nbSites_;
    nbStates_          = data.nbStates_;
    nbClasses_         = data.nbClasses_;
    nbDistinctSites_   = data.nbDistinctSites_;
    if (data.shrunkData_)
      shrunkData_ = std::shared_ptr<AlignmentDataInterface>(data.shrunkData_->clone());
    else
      shrunkData_      = nullptr;
    return *this;
  }

  virtual ~DRASDRTreeLikelihoodData() {}

  DRASDRTreeLikelihoodData* clone() const { return new DRASDRTreeLikelihoodData(*this); }

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
  void setTree(std::shared_ptr< const TreeTemplate<Node> > tree)
  {
    tree_ = tree;
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

  DRASDRTreeLikelihoodNodeData& getNodeData(int nodeId)
  {
    return nodeData_[nodeId];
  }

  const DRASDRTreeLikelihoodNodeData& getNodeData(int nodeId) const
  {
    return nodeData_[nodeId];
  }

  DRASDRTreeLikelihoodLeafData& getLeafData(int nodeId)
  {
    return leafData_[nodeId];
  }

  const DRASDRTreeLikelihoodLeafData& getLeafData(int nodeId) const
  {
    return leafData_[nodeId];
  }

  size_t getArrayPosition(int parentId, int sonId, size_t currentPosition) const
  {
    return currentPosition;
  }

  const std::map<int, VVVdouble>& getLikelihoodArrays(int nodeId) const
  {
    return nodeData_[nodeId].getLikelihoodArrays();
  }

  std::map<int, VVVdouble>& getLikelihoodArrays(int nodeId)
  {
    return nodeData_[nodeId].getLikelihoodArrays();
  }

  VVVdouble& getLikelihoodArray(int parentId, int neighborId)
  {
    return nodeData_[parentId].getLikelihoodArrayForNeighbor(neighborId);
  }

  const VVVdouble& getLikelihoodArray(int parentId, int neighborId) const
  {
    return nodeData_[parentId].getLikelihoodArrayForNeighbor(neighborId);
  }

  Vdouble& getDLikelihoodArray(int nodeId)
  {
    return nodeData_[nodeId].getDLikelihoodArray();
  }

  const Vdouble& getDLikelihoodArray(int nodeId) const
  {
    return nodeData_[nodeId].getDLikelihoodArray();
  }

  Vdouble& getD2LikelihoodArray(int nodeId)
  {
    return nodeData_[nodeId].getD2LikelihoodArray();
  }

  const Vdouble& getD2LikelihoodArray(int nodeId) const
  {
    return nodeData_[nodeId].getD2LikelihoodArray();
  }

  VVdouble& getLeafLikelihoods(int nodeId)
  {
    return leafData_[nodeId].getLikelihoodArray();
  }

  const VVdouble& getLeafLikelihoods(int nodeId) const
  {
    return leafData_[nodeId].getLikelihoodArray();
  }

  VVVdouble& getRootLikelihoodArray() { return rootLikelihoods_; }
  const VVVdouble& getRootLikelihoodArray() const { return rootLikelihoods_; }

  VVdouble& getRootSiteLikelihoodArray() { return rootLikelihoodsS_; }
  const VVdouble& getRootSiteLikelihoodArray() const { return rootLikelihoodsS_; }

  Vdouble& getRootRateSiteLikelihoodArray() { return rootLikelihoodsSR_; }
  const Vdouble& getRootRateSiteLikelihoodArray() const { return rootLikelihoodsSR_; }

  size_t getNumberOfDistinctSites() const { return nbDistinctSites_; }

  size_t getNumberOfSites() const { return nbSites_; }

  size_t getNumberOfStates() const { return nbStates_; }

  size_t getNumberOfClasses() const { return nbClasses_; }

  std::shared_ptr<const AlignmentDataInterface> getShrunkData() const { return shrunkData_; }
  
  const AlignmentDataInterface& shrunkData() const { return *shrunkData_; }

  /**
   * @brief Resize and initialize all likelihood arrays according to the given data set and substitution model.
   *
   * @param sites The sequences to use as data.
   * @param model The substitution model to use.
   * @throw Exception if an error occures.
   */
  void initLikelihoods(const AlignmentDataInterface& sites, const TransitionModelInterface& model);

  /**
   * @brief Rebuild likelihood arrays at inner nodes.
   *
   * This method is to be called when the topology of the tree has changed.
   * Node arrays relationship are rebuilt according to the new topology of the tree.
   * The leaves likelihood remain unchanged, so as for the first and second order derivatives.
   */
  void reInit();

  void reInit(const Node* node);

protected:
  /**
   * @brief This method initializes the leaves according to a sequence container.
   *
   * Here the container shrunkData_ is used.
   * Likelihood is set to 1 for the state corresponding to the sequence site,
   * otherwise it is set to 0.
   *
   * All likelihood arrays at each nodes are initialized according to alphabet
   * size and sequences length, and filled with 1.
   *
   * NB: This method is recursive.
   *
   * @param node  The node defining the subtree to analyse.
   * @param sites The sequence container to use.
   * @param model The model, used for initializing leaves' likelihoods.
   */
  void initLikelihoods(const Node* node, const AlignmentDataInterface& sites, const TransitionModelInterface& model);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_DRASDRTREELIKELIHOODDATA_H
