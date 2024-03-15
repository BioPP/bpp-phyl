// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Number.h>
#include <Bpp/Text/TextTools.h>

#include "../Tree/TreeTemplateTools.h"
#include "../Tree/TreeTools.h"
#include "AbstractTreeDrawing.h"

// From the STL:
#include <algorithm>

using namespace bpp;
using namespace std;

const TreeDrawingSettings AbstractTreeDrawing::DEFAULT_SETTINGS;

Point2D<double> AbstractTreeDrawing::getNodePosition(int nodeId) const
{
  vector<INode*> nodes = tree_->getNodes();
  for (unsigned int i = 0; i < nodes.size(); i++)
  {
    INode* node = nodes[i];
    if (node->getId() == nodeId)
    {
      return node->getInfos().getPosition();
    }
  }
  throw NodeNotFoundException("AbstractTreeDrawing::getNodePosition.", TextTools::toString(nodeId));
}

int AbstractTreeDrawing::getNodeAt(const Point2D<double>& position) const
{
  vector<INode*> nodes = tree_->getNodes();
  for (unsigned int i = 0; i < nodes.size(); i++)
  {
    INode* node = nodes[i];
    Point2D<double> nodePos = node->getInfos().getPosition();
    if (belongsTo(position, nodePos))
    {
      return node->getId();
    }
  }
  throw NodeNotFoundException("AbstractTreeDrawing::getNode.", "");
}

bool AbstractTreeDrawing::belongsTo(const Point2D<double>& p1, const Point2D<double>& p2) const
{
  return p1.getX() >= p2.getX() - settings_->pointArea && p1.getX() <= p2.getX() + settings_->pointArea
         && p1.getY() >= p2.getY() - settings_->pointArea && p1.getY() <= p2.getY() + settings_->pointArea;
}

void AbstractTreeDrawing::drawAtNode(GraphicDevice& gDevice, const INode& node, const string& text, double xOffset, double yOffset, short hpos, short vpos, double angle) const
{
  gDevice.drawText(node.getInfos().getX() + xOffset * xUnit_, node.getInfos().getY() + yOffset * yUnit_, text, hpos, vpos, angle);
}

void AbstractTreeDrawing::drawAtBranch(GraphicDevice& gDevice, const INode& node, const string& text, double xOffset, double yOffset, short hpos, short vpos, double angle) const
{
  if (node.hasFather())
  {
    gDevice.drawText((node.getInfos().getX() + node.getFather()->getInfos().getX()) / 2. + xOffset * xUnit_, node.getInfos().getY() + yOffset * yUnit_, text, hpos, vpos, angle);
  }
}
