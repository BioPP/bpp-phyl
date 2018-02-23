//
// File: TreeDrawingDisplayControler.cpp
// Created by: Julien Dutheil
// Created on: Tue May 18 12:37 2010
//

/*
Copyright or © or Copr. Bio++ Development Team, (2010)

This software is a computer program whose purpose is to provide
graphic components to develop bioinformatics applications.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "TreeDrawingDisplayControler.h"

using namespace bpp;
using namespace std;

TreeDrawingDisplayControler::~TreeDrawingDisplayControler()
{
  for (std::map<std::string, TreeDrawingListener*>::iterator it = listeners_.begin();
    it != listeners_.end(); ++it)
  {
    for (unsigned int j = 0; j < registeredTreeDrawings_.size(); ++j)
      registeredTreeDrawings_[j]->removeTreeDrawingListener(it->second);
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
  //Add listener to already registered tree drawings:
  for (unsigned int i = 0; i < registeredTreeDrawings_.size(); ++i)
    registeredTreeDrawings_[i]->addTreeDrawingListener(listener);
}

const string BasicTreeDrawingDisplayControler::PROPERTY_NODE_IDS         = "Node ids";
const string BasicTreeDrawingDisplayControler::PROPERTY_LEAF_NAMES       = "Leaf names";
const string BasicTreeDrawingDisplayControler::PROPERTY_BRANCH_LENGTHS   = "Branch lengths";
const string BasicTreeDrawingDisplayControler::PROPERTY_BOOTSTRAP_VALUES = "Bootstrap values";

