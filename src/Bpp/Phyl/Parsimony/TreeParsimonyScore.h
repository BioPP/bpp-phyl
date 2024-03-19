// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_PARSIMONY_TREEPARSIMONYSCORE_H
#define BPP_PHYL_PARSIMONY_TREEPARSIMONYSCORE_H

#include <Bpp/Clonable.h>

#include "../Tree/TreeTemplate.h"

// From the STL:
#include <vector>

namespace bpp
{
/**
 * @brief Compute a parsimony score.
 */
class TreeParsimonyScoreInterface :
  public virtual Clonable
{
public:
  TreeParsimonyScoreInterface() {}
  virtual ~TreeParsimonyScoreInterface() {}

  TreeParsimonyScoreInterface* clone() const = 0;

public:
  /**
   * @brief Get the score for the current tree, i.e. the total minimum number of changes in the tree.
   *
   * @return The minimum total number of changes in the tree.
   */
  virtual unsigned int getScore() const = 0;

  /**
   * @brief Get the score for a given site for the current tree, i.e. the total minimum number of changes in the tree for this site.
   *
   * @param site The corresponding site.
   * @return The minimum total number of changes in the tree for site 'site'.
   */
  virtual unsigned int getScoreForSite(size_t site) const = 0;

  /**
   * @brief Get the score for each site for the current tree, i.e. the total minimum number of changes in the tree for each site.
   *
   * @return The minimum total number of changes in the tree for each site.
   */
  virtual std::vector<unsigned int> getScorePerSite() const = 0;

  /**
   * @brief Get the tree for wich scores are computed.
   *
   * @return The tree associated to this object.
   */
  virtual const Tree& tree() const = 0;

  /**
   * @brief Get the state map associated to this instance.
   *
   * @return The underlying state map.
   */
  virtual const StateMapInterface& stateMap() const = 0;

  /**
   * @brief Get the state map associated to this instance.
   *
   * @return A shared pointer toward the underlying state map.
   */
  virtual std::shared_ptr<const StateMapInterface> getStateMap() const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_PARSIMONY_TREEPARSIMONYSCORE_H
