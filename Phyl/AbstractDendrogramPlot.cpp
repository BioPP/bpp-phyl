//
// File: AbstractDendrogramPlot.cpp
// Created by: Julien Dutheil
// Created on: Fri Jul 17 11:23 2009
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#include "AbstractDendrogramPlot.h"

using namespace bpp;
using namespace std;

short AbstractDendrogramPlot::ORIENTATION_LEFT_TO_RIGHT = 1;
short AbstractDendrogramPlot::ORIENTATION_RIGHT_TO_LEFT = 2;
short AbstractDendrogramPlot::ORIENTATION_TOP_TO_BOTTOM = 3;
short AbstractDendrogramPlot::ORIENTATION_BOTTOM_TO_TOP = 4;

string AbstractDendrogramPlot::PROPERTY_BRLEN = "Branch lengths";
string AbstractDendrogramPlot::PROPERTY_BOOTSTRAP = "Bootstrap values";
string AbstractDendrogramPlot::PROPERTY_IDS = "Nodes ids";

void AbstractDendrogramPlot::plot(GraphicDevice& gDevice) const throw (Exception)
{
  drawDendrogram_(gDevice);
  if (getDisplaySettings().drawLeafNames)
    drawLeafNames_(gDevice);
}

bool AbstractDendrogramPlot::drawProperty(GraphicDevice& gDevice, const string& property) const
{
  if (property == PROPERTY_BRLEN)
  {
    drawBranchLengthValues_(gDevice);
    return true;
  }
  if (property == PROPERTY_BOOTSTRAP)
  {
    drawBootstrapValues_(gDevice);
    return true;
  }
  if (property == PROPERTY_IDS)
  {
    drawNodesId_(gDevice);
    return true;
  }
  return false;
}

void AbstractDendrogramPlot::drawNodesId_(GraphicDevice& gDevice) const
{
  if (!hasTree()) return;
  vector<const INode*> nodes = getTree_()->getNodes();
  for (unsigned int i = 0; i < nodes.size(); i++)
  {
    const INode* node = nodes[i];
    drawAtNode(gDevice, *node, TextTools::toString(node->getId()), 0, 0,
        horOrientation_ == ORIENTATION_LEFT_TO_RIGHT
        ? GraphicDevice::TEXT_HORIZONTAL_LEFT
        : GraphicDevice::TEXT_HORIZONTAL_RIGHT);
  }
}

void AbstractDendrogramPlot::drawLeafNames_(GraphicDevice& gDevice) const
{
  if (!hasTree()) return;
  Font fontBck = gDevice.getCurrentFont();
  gDevice.setCurrentFont(getDisplaySettings().fontLeafNames);
  vector<const INode *> leaves = getTree_()->getLeaves();
  for(unsigned int i = 0; i < leaves.size(); i++)
  {
    const INode* node = leaves[i];
    drawAtNode(gDevice, *node, TextTools::toString(node->getName()), 0, 0,
        horOrientation_ == ORIENTATION_LEFT_TO_RIGHT
        ? GraphicDevice::TEXT_HORIZONTAL_LEFT
        : GraphicDevice::TEXT_HORIZONTAL_RIGHT);
  }
  gDevice.setCurrentFont(fontBck);
}

void AbstractDendrogramPlot::drawBranchLengthValues_(GraphicDevice& gDevice) const
{
  if (!hasTree()) return;
  vector<const INode*> nodes = getTree_()->getNodes();
  for (unsigned int i = 0; i < nodes.size(); i++)
  {
    const INode* node = nodes[i];
    if (node->hasDistanceToFather())
    {
      drawAtBranch(gDevice, *node, TextTools::toString(node->getDistanceToFather()), 0, 0,GraphicDevice::TEXT_HORIZONTAL_CENTER);
    }
  }
}

void AbstractDendrogramPlot::drawBootstrapValues_(GraphicDevice& gDevice) const
{
  if (!hasTree()) return;
  vector<const INode *> nodes = getTree_()->getNodes();
  for (unsigned int i = 0; i < nodes.size(); i++)
  {
    const INode* node = nodes[i];
    if (node->hasBranchProperty(TreeTools::BOOTSTRAP))
    {
      const Clonable* b = node->getBranchProperty(TreeTools::BOOTSTRAP);
      drawAtNode(gDevice, *node, TextTools::toString(dynamic_cast<const Number<double> *>(b)->getValue()),
        horOrientation_ == ORIENTATION_LEFT_TO_RIGHT ? -5 : 5, 0,
        horOrientation_ == ORIENTATION_LEFT_TO_RIGHT
        ? GraphicDevice::TEXT_HORIZONTAL_LEFT
        : GraphicDevice::TEXT_HORIZONTAL_RIGHT);
    }
  }
}

