// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_TREE_NNISEARCHABLE_H
#define BPP_PHYL_TREE_NNISEARCHABLE_H


#include "Node.h"
#include "TopologySearch.h"
#include "TreeTemplate.h"

namespace bpp
{
/**
 * @brief Interface for Nearest Neighbor Interchanges algorithms.
 *
 * This interface defines the methods to work with NNI algorithms.
 *
 * NNISearchable objects are supposed to work with TreeTemplate objects.
 * NNI are defined as follow:
 * <pre>
 * ------------->
 *     +------- C
 *     |
 * D --+ X  +-- B
 *     |    |
 *     +----+ F
 *          |
 *          +-- A
 * </pre>
 * Where:
 * -F is the focus (parent) node,
 * -A and B are the son of F
 * -X is the parent of F and so on.
 * Two NNI's are possible for branch (XF):
 * - swaping B and C, which is the same as D and A
 * - swaping A and C, which is the same as D and B
 * Because of the rooted representation, we'll consider B \f$\leftrightarrow\f$ C and A \f$\leftrightarrow\f$ C,
 * which are simpler to perform.
 *
 * For unrooted tree, we have at the 'root' node:
 * <pre>
 * ------------->
 *   +------- D
 *   |
 * X +------- C
 *   |
 *   |    +-- B
 *   |    |
 *   +----+ F
 *        |
 *        +-- A
 * </pre>
 * In this case, we swap A or B with one of C or D.
 * Which one of C or D depends on the implementation, but will always be the same,
 * so that swapping A or B involve 2 distinct NNI.
 *
 */
class NNISearchable :
  public TopologyListener,
  public virtual Clonable
{
public:
  NNISearchable() {}
  virtual ~NNISearchable() {}

  virtual NNISearchable* clone() const = 0;

public:
  /**
   * @brief Send the score of a NNI movement, without performing it.
   *
   * This methods sends the score variation.
   * This variation must be negative if the new point is better,
   * i.e. the object is to be used with a minimizing optimization
   * (for consistence with Optimizer objects).
   *
   * @param nodeId The id of the node defining the NNI movement.
   * @return The score variation of the NNI.
   * @throw NodeException If the node does not define a valid NNI.
   */
  virtual double testNNI(int nodeId) const = 0;

  /**
   * @brief Perform a NNI movement.
   *
   * @param nodeId The id of the node defining the NNI movement.
   * @throw NodeException If the node does not define a valid NNI.
   */
  virtual void doNNI(int nodeId) = 0;

  /**
   * @brief Get the tree associated to this NNISearchable object.
   *
   * @return The tree associated to this instance.
   */
  virtual const Tree& topology() const = 0;

  /**
   * @brief Get the current score of this NNISearchable object.
   *
   * @return The current score of this instance.
   */
  virtual double getTopologyValue() const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_TREE_NNISEARCHABLE_H
