//
// File: AbstractTreeDrawing.cpp
// Created by: Julien Dutheil
// Created on: Sun Oct 16 10:06 2006
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "AbstractTreeDrawing.h"
#include "../Tree/TreeTemplateTools.h"
#include "../Tree/TreeTools.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/Number.h>

// From the STL:
#include <algorithm>

using namespace bpp;
using namespace std;

const TreeDrawingSettings AbstractTreeDrawing::DEFAULT_SETTINGS;
  
Point2D<double> AbstractTreeDrawing::getNodePosition(int nodeId) const
throw (NodeNotFoundException)
{
  vector<INode *> nodes = tree_->getNodes();
  for(unsigned int i = 0; i < nodes.size(); i++)
  {
    INode * node = nodes[i];
    if(node->getId() == nodeId)
    {
      return node->getInfos().getPosition();
    }
  }
  throw NodeNotFoundException("AbstractTreeDrawing::getNodePosition.", TextTools::toString(nodeId));
}

int AbstractTreeDrawing::getNodeAt(const Point2D<double>& position) const
throw (NodeNotFoundException)
{
  vector<INode *> nodes = tree_->getNodes();
  for(unsigned int i = 0; i < nodes.size(); i++)
  {
    INode * node = nodes[i];
    Point2D<double> nodePos = node->getInfos().getPosition();
    if(belongsTo(position, nodePos))
    {
      return node->getId();
    }
  }
  throw NodeNotFoundException("AbstractTreeDrawing::getNode.", "");
}

bool AbstractTreeDrawing::belongsTo(const Point2D<double>& p1, const Point2D<double>& p2) const
{
  return (p1.getX() >= p2.getX() - settings_->pointArea && p1.getX() <= p2.getX() + settings_->pointArea
       && p1.getY() >= p2.getY() - settings_->pointArea && p1.getY() <= p2.getY() + settings_->pointArea);
}

void AbstractTreeDrawing::drawAtNode(GraphicDevice& gDevice, const INode& node, const string& text, double xOffset, double yOffset, short hpos, short vpos, double angle) const
{
  gDevice.drawText(node.getInfos().getX() + xOffset * xUnit_, node.getInfos().getY() + yOffset * yUnit_, text, hpos, vpos, angle);
}

void AbstractTreeDrawing::drawAtBranch(GraphicDevice& gDevice, const INode& node, const string& text, double xOffset, double yOffset, short hpos, short vpos, double angle) const
{
  if (node.hasFather()) {
    gDevice.drawText((node.getInfos().getX() + node.getFather()->getInfos().getX()) / 2. + xOffset * xUnit_, node.getInfos().getY() + yOffset * yUnit_, text, hpos, vpos, angle);
  }
}  

