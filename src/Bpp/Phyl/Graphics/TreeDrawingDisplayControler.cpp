// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "TreeDrawingDisplayControler.h"

using namespace bpp;
using namespace std;

TreeDrawingDisplayControler::~TreeDrawingDisplayControler()
{
  for (std::map<std::string, TreeDrawingListener*>::iterator it = listeners_.begin();
       it != listeners_.end(); ++it)
  {
    for (unsigned int j = 0; j < registeredTreeDrawings_.size(); ++j)
    {
      registeredTreeDrawings_[j]->removeTreeDrawingListener(it->second);
    }
    delete it->second;
  }
}

void TreeDrawingDisplayControler::addListener(const std::string& propertyName, TreeDrawingListener* listener)
{
  if (listeners_.find(propertyName) != listeners_.end())
    throw Exception("TreeDrawingDisplayControler::addListener. A listener is already added with this name: " + propertyName + ".");
  if (!listener)
    throw Exception("TreeDrawingDisplayControler::addListener. Trying to add a NULL listener!");
  if (!listener->isAutonomous())
    throw Exception("TreeDrawingDisplayControler::addListener. Trying to add a non-autonomous listener!");
  listeners_[propertyName] = listener;
  // Add listener to already registered tree drawings:
  for (unsigned int i = 0; i < registeredTreeDrawings_.size(); ++i)
  {
    registeredTreeDrawings_[i]->addTreeDrawingListener(listener);
  }
}

const string BasicTreeDrawingDisplayControler::PROPERTY_NODE_IDS         = "Node ids";
const string BasicTreeDrawingDisplayControler::PROPERTY_LEAF_NAMES       = "Leaf names";
const string BasicTreeDrawingDisplayControler::PROPERTY_BRANCH_LENGTHS   = "Branch lengths";
const string BasicTreeDrawingDisplayControler::PROPERTY_BOOTSTRAP_VALUES = "Bootstrap values";
