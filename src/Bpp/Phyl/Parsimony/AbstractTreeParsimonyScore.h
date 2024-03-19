// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_PARSIMONY_ABSTRACTTREEPARSIMONYSCORE_H
#define BPP_PHYL_PARSIMONY_ABSTRACTTREEPARSIMONYSCORE_H


#include "../Model/StateMap.h"
#include "TreeParsimonyScore.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainer.h>

namespace bpp
{
/**
 * @brief Partial implementation of the TreeParsimonyScore interface.
 */
class AbstractTreeParsimonyScore :
  public virtual TreeParsimonyScoreInterface
{
private:
  std::shared_ptr<TreeTemplate<Node>> treePtr_;
  std::shared_ptr<const SiteContainerInterface> data_;
  std::shared_ptr<const Alphabet> alphabet_;
  std::shared_ptr<const StateMapInterface> statesMap_;
  size_t nbStates_;

public:
  AbstractTreeParsimonyScore(
      std::shared_ptr<TreeTemplate<Node>> tree,
      std::shared_ptr<const SiteContainerInterface> data,
      bool verbose,
      bool includeGaps);

  AbstractTreeParsimonyScore(
      std::shared_ptr<TreeTemplate<Node>> tree,
      std::shared_ptr<const SiteContainerInterface> data,
      std::shared_ptr<const StateMapInterface> statesMap,
      bool verbose);

  virtual ~AbstractTreeParsimonyScore() {}

private:
  void init_(std::shared_ptr<const SiteContainerInterface> data, bool verbose);

public:
  const Tree& tree() const override { return *treePtr_; }
  virtual const TreeTemplate<Node>& treeTemplate() const { return *treePtr_; }
  virtual std::shared_ptr<const TreeTemplate<Node>> getTreeTemplate() const { return treePtr_; }
  std::vector<unsigned int> getScorePerSite() const override;
  const StateMapInterface& stateMap() const override { return *statesMap_; }
  std::shared_ptr<const StateMapInterface> getStateMap() const override { return statesMap_; }

protected:
  virtual Tree& tree_() { return *treePtr_; }
  virtual TreeTemplate<Node>& treeTemplate_() { return *treePtr_; }
  virtual std::shared_ptr<TreeTemplate<Node>> getTreeTemplate_() { return treePtr_; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_PARSIMONY_ABSTRACTTREEPARSIMONYSCORE_H
