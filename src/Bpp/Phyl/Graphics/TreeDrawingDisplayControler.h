// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_GRAPHICS_TREEDRAWINGDISPLAYCONTROLER_H
#define BPP_PHYL_GRAPHICS_TREEDRAWINGDISPLAYCONTROLER_H


#include "TreeDrawingListener.h"

// From the STL:
#include <string>
#include <vector>
#include <map>
#include <algorithm>

namespace bpp
{
/**
 * @brief Easy tune of tree drawings display.
 *
 * This class maintains a set of autonomous TreeDrawing listeners that
 * are used for annotating a tree drawing.
 *
 * @author Julien Dutheil
 */
class TreeDrawingDisplayControler
{
private:
  std::map<std::string, TreeDrawingListener*> listeners_;
  std::vector<TreeDrawing*> registeredTreeDrawings_;

public:
  TreeDrawingDisplayControler() :
    listeners_(), registeredTreeDrawings_()
  {}

private:
  TreeDrawingDisplayControler(const TreeDrawingDisplayControler& tddc) :
    listeners_(), registeredTreeDrawings_(tddc.registeredTreeDrawings_)
  {
    for (std::map<std::string, TreeDrawingListener*>::const_iterator it = tddc.listeners_.begin();
         it != tddc.listeners_.end(); ++it)
    {
      listeners_[it->first] = dynamic_cast<TreeDrawingListener*>(it->second->clone());
    }
  }
  TreeDrawingDisplayControler& operator=(const TreeDrawingDisplayControler& tddc)
  {
    listeners_.clear();
    registeredTreeDrawings_ = tddc.registeredTreeDrawings_;
    for (std::map<std::string, TreeDrawingListener*>::const_iterator it = tddc.listeners_.begin();
         it != tddc.listeners_.end(); ++it)
    {
      listeners_[it->first] = dynamic_cast<TreeDrawingListener*>(it->second->clone());
    }
    return *this;
  }

public:
  virtual ~TreeDrawingDisplayControler();

public:
  /**
   * @brief Add a listener to the controler. The controler then owns the object, and will
   * copy or delete it when needed.
   */
  void addListener(const std::string& propertyName, TreeDrawingListener* listener);

  bool hasListenerFor(const std::string& propertyName) const
  {
    return listeners_.find(propertyName) != listeners_.end();
  }

  void enableListener(const std::string& propertyName, bool tf)
  {
    if (!hasListenerFor(propertyName))
      throw Exception("TreeDrawingDisplayControler::enableListener. No listener is registered for property " + propertyName + ".");
    listeners_[propertyName]->enable(tf);
  }

  bool isListenerEnabled(const std::string& propertyName) const
  {
    if (!hasListenerFor(propertyName))
      throw Exception("TreeDrawingDisplayControler::enableListener. No listener is registered for property " + propertyName + ".");
    return listeners_.find(propertyName)->second->isEnabled();
  }

  void registerTreeDrawing(TreeDrawing* td)
  {
    if (std::find(registeredTreeDrawings_.begin(), registeredTreeDrawings_.end(), td) != registeredTreeDrawings_.end())
      throw Exception("TreeDrawingDisplayControler::registerTreeDrawing. TreeDrawing is already associated to this controler.");
    for (std::map<std::string, TreeDrawingListener*>::iterator it = listeners_.begin();
         it != listeners_.end(); ++it)
    {
      td->addTreeDrawingListener(it->second);
    }
    registeredTreeDrawings_.push_back(td);
  }
};


/**
 * @brief Easy tune of tree drawings display, a basic implementation:
 *
 * This class maintains several "standard" drawing listener for:
 * - Plotting node id,
 * - Plotting leaves names,
 * - Plotting branch lengths,
 * - Plotting plotting bootstrap values.
 *
 * This controler takes as an argument a TreeDrawingSettings object that is used by
 * all listeners that require one.
 */
class BasicTreeDrawingDisplayControler :
  public TreeDrawingDisplayControler
{
public:
  static const std::string PROPERTY_NODE_IDS;
  static const std::string PROPERTY_LEAF_NAMES;
  static const std::string PROPERTY_BRANCH_LENGTHS;
  static const std::string PROPERTY_BOOTSTRAP_VALUES;

private:
  const TreeDrawingSettings* settings_;

public:
  BasicTreeDrawingDisplayControler(const TreeDrawingSettings* settings) :
    settings_(settings)
  {
    if (!settings)
      throw NullPointerException("BasicTreeDrawingDisplayControler::constructor. Trying to use NULL settings.");
    addListener(PROPERTY_NODE_IDS, new NodesIdTreeDrawingListener        (settings_, true));
    addListener(PROPERTY_LEAF_NAMES, new LeafNamesTreeDrawingListener      (settings_, true));
    addListener(PROPERTY_BRANCH_LENGTHS, new BranchLengthsTreeDrawingListener  (settings_, true));
    addListener(PROPERTY_BOOTSTRAP_VALUES, new BootstrapValuesTreeDrawingListener(settings_, true));
  }

private:
  BasicTreeDrawingDisplayControler(const BasicTreeDrawingDisplayControler&) : TreeDrawingDisplayControler(), settings_(0) {}
  BasicTreeDrawingDisplayControler& operator=(const BasicTreeDrawingDisplayControler&) { return *this; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_GRAPHICS_TREEDRAWINGDISPLAYCONTROLER_H
