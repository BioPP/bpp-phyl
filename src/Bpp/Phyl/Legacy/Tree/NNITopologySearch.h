// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_TREE_NNITOPOLOGYSEARCH_H
#define BPP_PHYL_LEGACY_TREE_NNITOPOLOGYSEARCH_H


#include "../../Tree/NNISearchable.h"
#include "../../Tree/TopologySearch.h"

namespace bpp
{
/**
 * @brief NNI topology search method.
 *
 * Several algorithm are implemented:
 * - Fast algorithm: loop over all nodes, check all NNIs and perform the corresponding change if it improve the score.
 *   When a NNI is done, reloop from the first node without checking the remaining ones.
 * - Better algorithm: loop over all nodes, check all NNIS.
 *   Then choose the NNI corresponding to the best improvement and perform it.
 *   Then re-loop over all nodes.
 * - PhyML algorithm (not fully tested, use with care): as the previous one, but perform all NNI improving the score at the same time.
 *   Leads to faster convergence.
 */
class NNITopologySearch :
  public virtual TopologySearch
{
public:
  const static std::string FAST;
  const static std::string BETTER;
  const static std::string PHYML;

private:
  std::shared_ptr<NNISearchable> searchableTree_;
  std::string algorithm_;
  unsigned int verbose_;
  std::vector<std::shared_ptr<TopologyListener>> topoListeners_;

public:
  NNITopologySearch(
      std::shared_ptr<NNISearchable> tree,
      const std::string& algorithm = FAST,
      unsigned int verbose = 2) :
    searchableTree_(tree), algorithm_(algorithm), verbose_(verbose), topoListeners_()
  {}

  NNITopologySearch(const NNITopologySearch& ts) :
    searchableTree_(ts.searchableTree_),
    algorithm_(ts.algorithm_),
    verbose_(ts.verbose_),
    topoListeners_(ts.topoListeners_)
  {
    // Hard-copy all listeners:
    for (unsigned int i = 0; i < topoListeners_.size(); i++)
    {
      topoListeners_[i].reset(ts.topoListeners_[i]->clone());
    }
  }

  NNITopologySearch& operator=(const NNITopologySearch& ts)
  {
    searchableTree_ = ts.searchableTree_;
    algorithm_      = ts.algorithm_;
    verbose_        = ts.verbose_;
    topoListeners_  = ts.topoListeners_;
    // Hard-copy all listeners:
    for (unsigned int i = 0; i < topoListeners_.size(); i++)
    {
      topoListeners_[i].reset(ts.topoListeners_[i]->clone());
    }
    return *this;
  }


  virtual ~NNITopologySearch()
  {}

public:
  void search();

  /**
   * @brief Add a listener to the list.
   *
   * All listeners will be notified in the order of the list.
   * The first listener to be notified is the NNISearchable object itself.
   *
   * The listener will be owned by this instance, and copied when needed.
   */
  void addTopologyListener(std::shared_ptr<TopologyListener> listener)
  {
    if (listener)
      topoListeners_.push_back(listener);
  }

public:
  /**
   * @brief Retrieve the tree.
   *
   * @return The tree associated to this instance.
   */
  const Tree& topology() const { return searchableTree_->topology(); }

  /**
   * @return The NNISearchable object associated to this instance.
   */
  std::shared_ptr<NNISearchable> getSearchableObject() { return searchableTree_; }

  /**
   * @return The NNISearchable object associated to this instance.
   */
  std::shared_ptr<const NNISearchable> getSearchableObject() const { return searchableTree_; }

protected:
  void searchFast();
  void searchBetter();
  void searchPhyML();

  /**
   * @brief Process a TopologyChangeEvent to all listeners.
   */
  void notifyAllPerformed(const TopologyChangeEvent& event);
  /**
   * @brief Process a TopologyChangeEvent to all listeners.
   */
  void notifyAllTested(const TopologyChangeEvent& event);
  /**
   * @brief Process a TopologyChangeEvent to all listeners.
   */
  void notifyAllSuccessful(const TopologyChangeEvent& event);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_TREE_NNITOPOLOGYSEARCH_H
