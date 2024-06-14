// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _LEGACY_ANCESTRAL_STATES_RECONSTRUCTION_H_
#define _LEGACY_ANCESTRAL_STATES_RECONSTRUCTION_H_

// From SeqLib:
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/Container/SiteContainer.h>

// From the STL:
#include <vector>
#include <map>

namespace bpp
{
class Node;

/**
 * @brief Interface for ancestral states reconstruction methods.
 */
class LegacyAncestralStateReconstruction
{
public:
  LegacyAncestralStateReconstruction() {}
  virtual ~LegacyAncestralStateReconstruction() {}

public:
  /**
   * @brief Get ancestral states for a given node as a vector of int.
   *
   * The size of the vector depends on the implementation.
   * This method is mainly for efficient internal use in other classes.
   * Consider using the getAncestralSequenceForNode() method for a more
   * general output.
   *
   * @param nodeId the id of the node at which the states must be reconstructed.
   * @return A vector of states indices.
   * @see getAncestralSequenceForNode
   */
  virtual std::vector<size_t> getAncestralStatesForNode(int nodeId) const = 0;

  /**
   * @brief Get all ancestral states for all nodes.
   *
   * Call the getAncestralSequenceForNode() method on each node in the tree.
   *
   * @return A map with nodes id as key, and a vector of states indices as value.
   * @see getAncestralSequenceForNode
   */
  virtual std::map<int, std::vector<size_t>> getAllAncestralStates() const = 0;

  /**
   * @brief Get the ancestral sequence for a given node.
   *
   * @param nodeId The id of the node at which the sequence must be reconstructed.
   * @return A sequence object.
   */
  virtual std::unique_ptr<Sequence> getAncestralSequenceForNode(int nodeId) const = 0;

  /**
   * @brief Get all the ancestral sequences for all nodes.
   *
   * @return A new SiteContainer object.
   */
  virtual std::unique_ptr<SiteContainerInterface> getAncestralSequences() const = 0;
};
} // end of namespace bpp.

#endif // _LEGACY_ANCESTRAL_STATES_RECONSTRUCTION_H_
