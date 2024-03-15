// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_GRAPHICS_TREEDRAWINGLISTENER_H
#define BPP_PHYL_GRAPHICS_TREEDRAWINGLISTENER_H

#include <Bpp/Clonable.h>

#include "TreeDrawing.h"

namespace bpp
{
/**
 * @brief Interface allowing to capture drawing events.
 *
 * Implementing this interface allows you to easily and efficiently tune a plot,
 * and/or add elements to it.
 */
class TreeDrawingListener :
  public virtual Clonable
{
public:
  TreeDrawingListener* clone() const = 0;

  virtual void beforeDrawTree(const DrawTreeEvent& event) = 0;
  virtual void afterDrawTree(const DrawTreeEvent& event) = 0;
  virtual void beforeDrawNode(const DrawNodeEvent& event) = 0;
  virtual void afterDrawNode(const DrawNodeEvent& event) = 0;
  virtual void beforeDrawBranch(const DrawBranchEvent& event) = 0;
  virtual void afterDrawBranch(const DrawBranchEvent& event) = 0;

  /**
   * @brief Tells if the listener is autonomous. If so, it
   * will never be hard-copied or deleted.
   */
  virtual bool isAutonomous() const = 0;

  virtual bool isEnabled() const = 0;
  virtual void enable(bool tf) = 0;
};


/**
 * @brief An empty implementation of the TreeDrawingListener interface.
 */
class TreeDrawingListenerAdapter :
  public virtual TreeDrawingListener
{
private:
  bool autonomous_;
  bool enabled_;

public:
  TreeDrawingListenerAdapter(bool autonomous) : autonomous_(autonomous), enabled_(true) {}

public:
  void beforeDrawTree(const DrawTreeEvent& event) {}
  void afterDrawTree(const DrawTreeEvent& event) {}
  void beforeDrawNode(const DrawNodeEvent& event) {}
  void afterDrawNode(const DrawNodeEvent& event) {}
  void beforeDrawBranch(const DrawBranchEvent& event) {}
  void afterDrawBranch(const DrawBranchEvent& event) {}
  bool isAutonomous() const { return autonomous_; }
  bool isEnabled() const { return enabled_; }
  void enable(bool tf) { enabled_ = tf; }
};


/**
 * @brief A TreeDrawingListener implementation that writes nodes id.
 */
class NodesIdTreeDrawingListener :
  public TreeDrawingListenerAdapter
{
private:
  const TreeDrawingSettings* settings_;

public:
  NodesIdTreeDrawingListener(const TreeDrawingSettings* settings, bool autonomous = false) :
    TreeDrawingListenerAdapter(autonomous),
    settings_(settings)
  {}

  NodesIdTreeDrawingListener(const NodesIdTreeDrawingListener& lntdl) :
    TreeDrawingListenerAdapter(lntdl),
    settings_(lntdl.settings_) {}

  NodesIdTreeDrawingListener& operator=(const NodesIdTreeDrawingListener& lntdl)
  {
    TreeDrawingListenerAdapter::operator=(lntdl);
    settings_ = lntdl.settings_;
    return *this;
  }

  NodesIdTreeDrawingListener* clone() const { return new NodesIdTreeDrawingListener(*this); }

public:
  void afterDrawNode(const DrawNodeEvent& event);
};


/**
 * @brief A TreeDrawingListener implementation that write leaf names.
 */
class LeafNamesTreeDrawingListener :
  public TreeDrawingListenerAdapter
{
private:
  const TreeDrawingSettings* settings_;

public:
  LeafNamesTreeDrawingListener(const TreeDrawingSettings* settings, bool autonomous = false) :
    TreeDrawingListenerAdapter(autonomous),
    settings_(settings)
  {}

  LeafNamesTreeDrawingListener(const LeafNamesTreeDrawingListener& lntdl) :
    TreeDrawingListenerAdapter(lntdl),
    settings_(lntdl.settings_)
  {}

  LeafNamesTreeDrawingListener& operator=(const LeafNamesTreeDrawingListener& lntdl)
  {
    TreeDrawingListenerAdapter::operator=(lntdl);
    settings_   = lntdl.settings_;
    return *this;
  }

  LeafNamesTreeDrawingListener* clone() const { return new LeafNamesTreeDrawingListener(*this); }

public:
  void afterDrawNode(const DrawNodeEvent& event);
};


/**
 * @brief A TreeDrawingListener implementation that write the branch lengths of inner nodes.
 *
 * Collapsed nodes are not labelled.
 *
 * This listener works with TreeDrawing classes, but is more efficient when used with a class that fires DrawINodeEvent events.
 */
class BranchLengthsTreeDrawingListener :
  public TreeDrawingListenerAdapter
{
private:
  const TreeDrawingSettings* settings_;

public:
  BranchLengthsTreeDrawingListener(const TreeDrawingSettings* settings, bool autonomous = false) :
    TreeDrawingListenerAdapter(autonomous),
    settings_(settings)
  {}

  BranchLengthsTreeDrawingListener(const BranchLengthsTreeDrawingListener& lntdl) :
    TreeDrawingListenerAdapter(lntdl),
    settings_(lntdl.settings_)
  {}

  BranchLengthsTreeDrawingListener& operator=(const BranchLengthsTreeDrawingListener& lntdl)
  {
    TreeDrawingListenerAdapter::operator=(lntdl);
    settings_   = lntdl.settings_;
    return *this;
  }

  BranchLengthsTreeDrawingListener* clone() const { return new BranchLengthsTreeDrawingListener(*this); }

public:
  void afterDrawBranch(const DrawBranchEvent& event);
};


/**
 * @brief A TreeDrawingListener implementation that write the bootstrap values of inner nodes.
 *
 * Collapsed nodes are not labelled.
 *
 * This listener works with TreeDrawing classes, but is more efficient when used with a class that fires DrawINodeEvent events.
 */
class BootstrapValuesTreeDrawingListener :
  public TreeDrawingListenerAdapter
{
private:
  const TreeDrawingSettings* settings_;

public:
  BootstrapValuesTreeDrawingListener(const TreeDrawingSettings* settings, bool autonomous = false) :
    TreeDrawingListenerAdapter(autonomous),
    settings_(settings)
  {}

  BootstrapValuesTreeDrawingListener(const BootstrapValuesTreeDrawingListener& lntdl) :
    TreeDrawingListenerAdapter(lntdl),
    settings_(lntdl.settings_) {}

  BootstrapValuesTreeDrawingListener& operator=(const BootstrapValuesTreeDrawingListener& lntdl)
  {
    TreeDrawingListenerAdapter::operator=(lntdl);
    settings_   = lntdl.settings_;
    return *this;
  }

  BootstrapValuesTreeDrawingListener* clone() const { return new BootstrapValuesTreeDrawingListener(*this); }

public:
  void afterDrawBranch(const DrawBranchEvent& event);
};


/**
 * @brief A TreeDrawingListener implementation that write the names of inner nodes.
 *
 * Collapsed nodes are not labelled.
 *
 * This listener works with TreeDrawing classes, but is more efficient when used with a class that fires DrawINodeEvent events.
 */
class LabelInnerNodesTreeDrawingListener :
  public TreeDrawingListenerAdapter
{
public:
  LabelInnerNodesTreeDrawingListener(bool autonomous = false) :
    TreeDrawingListenerAdapter(autonomous) {}

  LabelInnerNodesTreeDrawingListener* clone() const { return new LabelInnerNodesTreeDrawingListener(*this); }

public:
  void afterDrawNode(const DrawNodeEvent& event);
};


/**
 * @brief A TreeDrawingListener implementation that label the collapsed nodes.
 *
 * This listener works with TreeDrawing classes, but is more efficient when used with a class that fires DrawINodeEvent events.
 */
class LabelCollapsedNodesTreeDrawingListener :
  public TreeDrawingListenerAdapter
{
public:
  LabelCollapsedNodesTreeDrawingListener(bool autonomous = false) :
    TreeDrawingListenerAdapter(autonomous) {}

  LabelCollapsedNodesTreeDrawingListener* clone() const { return new LabelCollapsedNodesTreeDrawingListener(*this); }

public:
  void afterDrawNode(const DrawNodeEvent& event);
};
} // end of namespace bpp
#endif // BPP_PHYL_GRAPHICS_TREEDRAWINGLISTENER_H
