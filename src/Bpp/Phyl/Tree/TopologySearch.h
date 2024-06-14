// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_TREE_TOPOLOGYSEARCH_H
#define BPP_PHYL_TREE_TOPOLOGYSEARCH_H

#include <Bpp/Clonable.h>
#include <Bpp/Exceptions.h>


// From the STL:
#include <string>
#include <vector>

namespace bpp
{
/**
 * @brief Class for notifying new toplogy change events.
 */
class TopologyChangeEvent
{
protected:
  std::string message_;

public:
  TopologyChangeEvent() : message_("") {}
  TopologyChangeEvent(const std::string& message) : message_(message) {}
  virtual ~TopologyChangeEvent() {}

public:
  /**
   * @brief Get the message associated to this event.
   *
   * @return The message associated to this event.
   */
  virtual const std::string& getMessage() const { return message_; }
};

class TopologySearch;

/**
 * @brief Implement this interface to be notified when the topology of a tree
 * has changed during topology search.
 */
class TopologyListener :
  public virtual Clonable
{
public:
  TopologyListener() {}
  virtual ~TopologyListener() {}

  TopologyListener* clone() const = 0;

public:
  /**
   * @brief Notify a topology change event.
   *
   * This method is to be invoked after one or several NNI are performed.
   * It allows appropriate recomputations.
   *
   * In most case, this is the same as
   * topologyChangeTested() + topologyChangeSuccessful().
   *
   * @param event The topology change event.
   */
  virtual void topologyChangePerformed(const TopologyChangeEvent& event)
  {
    topologyChangeTested(event);
    topologyChangeSuccessful(event);
  }

  /**
   * @brief Notify a topology change event.
   *
   * @param event The topology change event.
   */
  virtual void topologyChangeTested(const TopologyChangeEvent& event) = 0;

  /**
   * @brief Tell that a topology change is definitive.
   *
   * This method is called after the topologyChangeTested() method.
   *
   * @param event The topology change event.
   */
  virtual void topologyChangeSuccessful(const TopologyChangeEvent& event) = 0;
};


/**
 * @brief Interface for topology search methods.
 */
class TopologySearch
{
public:
  TopologySearch() {}
  virtual ~TopologySearch() {}

public:
  /**
   * @brief Performs the search.
   */
  virtual void search() = 0;

  /**
   * @brief Add a topology listener to this class.
   *
   * TopologyListeners will be notified when the topology of the tree is modified.
   */
  virtual void addTopologyListener(std::shared_ptr<TopologyListener> listener) = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_TREE_TOPOLOGYSEARCH_H
