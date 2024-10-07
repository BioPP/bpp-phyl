// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_TREE_AWARENODE_H
#define BPP_PHYL_LEGACY_TREE_AWARENODE_H

#include <Bpp/BppString.h>
#include <Bpp/Clonable.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/Utils/MapTools.h>

#include "../../Tree/TreeExceptions.h"

// From the STL:
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <cstddef>
#include <memory>

#include <Bpp/Graph/AssociationTreeGraphImplObserver.h>

namespace bpp
{
/**
 * @brief A node class aware of its neighbours.
 *
 *
 * This class is made the more general as possible, while keeping it very simple. It contains:</p>
 * - An identity tag, to identity it in the tree;
 * - A pointer toward the father node;
 * - A std::vector of pointer toward son nodes;
 * - The distance from the father node:
 *
 * @see PhyloTree
 */

class AwareNode // :
    // public std::enable_shared_from_this<AwareNode>
{
protected:
  unsigned int id_;
  std::vector<AwareNode* > sons_;
  AwareNode* father_;
  double distanceToFather_;

public:
  /**
   * @brief Build a new void Node object.
   */
  AwareNode() :
    id_(0),
    sons_(),
    father_(0),
    distanceToFather_(0)
  {}

  /**
   * @brief Build a new Node with specified id.
   */
  AwareNode(unsigned int id) :
    id_(id),
    sons_(),
    father_(0),
    distanceToFather_(0)
  {}

  /**
   * @brief Copy constructor.
   *
   * @warning THIS OPERATOR COPIES ALL FIELDS, EXCEPTED FATHER AND
   * SON NODE POINTERS.
   *
   * Need the call of updateTree after whole Tree is built.
   *
   * @param node The node to copy.
   */
  AwareNode(const AwareNode& node);

  /**
   * @brief Assignation operator.
   *
   * @warning THIS OPERATOR COPIES ALL FIELDS, EXCEPTED FATHER AND
   * SON NODE POINTERS.
   *
   * Need the call of updateTree after whole Tree is built.
   *
   * @param node the node to copy.
   * @return A reference toward this node.
   */
  AwareNode& operator=(const AwareNode& node);

  AwareNode* clone() const { return new AwareNode(*this); }

  virtual ~AwareNode()
  {}

public:
  /**
   * @brief update information from TreeObserver
   *
   */
  template<class N, class E, class I>
  void updateTree(AssociationTreeGraphImplObserver<N, E, I>* tree, unsigned int index);


  /**
   * @name Identity
   *
   * @{
   */

  /**
   * @brief Get the node's id.
   *
   * @return The identity tag of this node.
   */
  virtual unsigned int getId() const { return id_; }

  /**
   * @brief Set this node's id.
   *
   * @param id The new identity tag.
   */
  virtual void setId(unsigned int id) { id_ = id; }

  virtual std::vector<unsigned int> getSonsId() const
  {
    std::vector<unsigned int> sonsId(sons_.size());
    for (size_t i = 0; i < sons_.size(); i++)
    {
      sonsId[i] = sons_[i]->getId();
    }
    return sonsId;
  }

  /** @} */

  /**
   * @name Distances:
   *
   * @{
   */

  /**
   * @brief Get the distance to the father node is there is one,
   * otherwise throw a NodeException.
   *
   * @return The distance to the father node.
   */
  virtual double getDistanceToFather() const
  {
    return distanceToFather_;
  }

  /**
   * @brief Set or update the distance toward the father node.
   *
   * Warning: a distance to the father node may be set even if no father node is specified.
   * This is used by several tree reconstruction methods.
   * It may also be useful for manipulating subtrees.
   *
   * @param distance The new distance to the father node.
   */
  virtual void setDistanceToFather(double distance)
  {
    distanceToFather_ = distance;
  }

  /** @} */

  /**
   * @name Father:
   *
   * @{
   */

  /**
   * @brief Get the father of this node is there is one.
   *
   * @return A pointer toward the father node, 0 if there is not.
   */
  virtual const AwareNode* getFather() const { return father_; }

  /**
   * @brief Get the father of this node is there is one.
   *
   * @return A pointer toward the father node, 0 if there is not.
   */
  virtual AwareNode* getFather() { return father_; }

  //  virtual int getFatherId() const { return father_->getId(); }

  /**
   * @brief Set the father node of this node.
   *
   * @param node The father node.
   */
  virtual void setFather(AwareNode* node)
  {
    father_ = node;
  }

  /**
   * @brief Remove the father of this node.
   */
  virtual void removeFather()
  {
    father_ = 0;
  }

  /**
   * @brief Tell if this node has a father node.
   */
  virtual bool hasFather() const { return father_ != 0; }

  /** @} */

  /**
   * @name Sons:
   *
   * @{
   */
  virtual size_t getNumberOfSons() const { return sons_.size(); }

  bool isLeaf() const
  {
    return getNumberOfSons() == 0;
  }

  virtual std::vector<AwareNode*>& getSons()
  {
    return sons_;
  }

  virtual const AwareNode* getSon(size_t pos) const
  {
    if (pos >= sons_.size()) throw IndexOutOfBoundsException("AwareNode::getSon().", pos, 0, sons_.size() - 1);
    return sons_[pos];
  }

  virtual AwareNode* getSon(size_t pos)
  {
    if (pos >= sons_.size()) throw IndexOutOfBoundsException("AwareNode::getSon().", pos, 0, sons_.size() - 1);
    return sons_[pos];
  }

  virtual void addSon(size_t pos, AwareNode* node)
  {
    if (!node)
      throw NullPointerException("AwareNode::addSon(). Empty node given as input.");
    if (find(sons_.begin(), sons_.end(), node) == sons_.end())
      sons_.insert(sons_.begin() + static_cast<ptrdiff_t>(pos), node);
    else // Otherwise node is already present.
      std::cerr << "DEVEL warning: AwareNode::addSon. Son node already registered! No pb here, but could be a bug in your implementation..." << std::endl;
  }

  virtual void addSon(AwareNode* node)
  {
    if (!node)
      throw NullPointerException("AwareNode::addSon(). Empty node given as input.");
    if (find(sons_.begin(), sons_.end(), node) == sons_.end())
      sons_.push_back(node);
    else // Otherwise node is already present.
      throw Exception("AwareNode::addSon. Trying to add a node which is already present.");
  }

  virtual void setSon(size_t pos, AwareNode* node)
  {
    if (!node)
      throw NullPointerException("AwareNode::setSon(). Empty node given as input.");
    if (pos >= sons_.size())
      throw IndexOutOfBoundsException("AwareNode::setSon(). Invalid node position.", pos, 0, sons_.size() - 1);
    std::vector<AwareNode*>::iterator search = find(sons_.begin(), sons_.end(), node);
    if (search == sons_.end() || search == sons_.begin() + static_cast<ptrdiff_t>(pos))
      sons_[pos] = node;
    else
      throw Exception("AwareNode::setSon. Trying to set a node which is already present.");
  }

  virtual void removeSon(size_t pos)
  {
    if (pos >= sons_.size())
      throw IndexOutOfBoundsException("AwareNode::removeSon(). Invalid node position.", pos, 0, sons_.size() - 1);
    sons_.erase(sons_.begin() + static_cast<ptrdiff_t>(pos));
  }

  virtual void removeSon(AwareNode* node)
  {
    if (!node)
      throw NullPointerException("AwareNode::removeSon(). Empty node given as input.");
    for (size_t i = 0; i < sons_.size(); i++)
    {
      if (sons_[i] == node)
      {
        sons_.erase(sons_.begin() + static_cast<ptrdiff_t>(i));
        return;
      }
    }
    throw Exception("AwareNode::removeSon, not found node : " + TextTools::toString(node->getId()));
  }

  virtual void removeSons()
  {
    sons_.clear();
  }

  virtual void swap(size_t branch1, size_t branch2);

  virtual size_t getSonPosition(const AwareNode* son) const;

  /** @} */

  /**
   * @name Operators:
   *
   * - a positive value send the corresponding son;
   * - a negative value send the father.
   *
   * @{
   */
  AwareNode* operator[](int i) { return (i < 0) ? father_ : sons_[static_cast<std::size_t>(i)]; }

  const AwareNode* operator[](int i) const { return (i < 0) ? father_ : sons_[static_cast<std::size_t>(i)]; }

  /** @} */

  virtual bool hasNoSon() const { return getNumberOfSons() == 0; }
};

template<class N, class E, class I>
void AwareNode::updateTree(AssociationTreeGraphImplObserver<N, E, I>* tree, unsigned int index)
{
  id_ = index;

  std::shared_ptr<N> thisN = tree->getNode(index);
  father_ = tree->hasFather(thisN) ? dynamic_cast<AwareNode*>(tree->getFatherOfNode(thisN).get()) : 0;

  sons_.clear();
  std::vector< std::shared_ptr<N>> vS = tree->getSons(thisN);

  for (typename std::vector<std::shared_ptr<N>>::iterator it = vS.begin(); it != vS.end(); it++)
  {
    sons_.push_back(dynamic_cast<AwareNode*>(it->get()));
  }
}
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_TREE_AWARENODE_H
