// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../Tree/TreeTemplateTools.h"
#include "CladogramPlot.h"

// From the STL:
#include <memory>

using namespace bpp;
using namespace std;

CladogramDrawBranchEvent::CladogramDrawBranchEvent(const TreeDrawing* source, GraphicDevice* gd, const INode* node, double length, const Cursor& cursor, short orientation) :
  DrawIBranchEvent(source, gd, node, cursor), orientation_(), length_(length)
{
  orientation_ = (orientation == AbstractDendrogramPlot::ORIENTATION_LEFT_TO_RIGHT ? -1 : 1);
}

Cursor CladogramDrawBranchEvent::getBranchCursor(double position) const
{
  double offset = 0;
  if (getINode()->hasDistanceToFather())
  {
    offset = orientation_ * length_ * position;
  }
  return getCursor().getTranslation(offset, 0);
}

void CladogramPlot::setTree(const Tree* tree)
{
  AbstractDendrogramPlot::setTree(tree);
  if (hasTree())
    totalDepth_ = static_cast<double>(TreeTemplateTools::getDepth(*getTree_()->getRootNode()));
}

void CladogramPlot::drawDendrogram_(GraphicDevice& gDevice) const
{
  if (hasTree())
  {
    DrawTreeEvent treeEvent(this, &gDevice);
    fireBeforeTreeEvent_(treeEvent);
    unsigned int tipCounter = 0;
    double y;
    recursivePlot_(gDevice, *const_cast<INode*>(getTree_()->getRootNode()),
        getHorizontalOrientation() == ORIENTATION_LEFT_TO_RIGHT ? 0 : getWidth() * getXUnit(),
        y,
        getHorizontalOrientation() == ORIENTATION_LEFT_TO_RIGHT ? 1. : -1.,
        getVerticalOrientation() == ORIENTATION_TOP_TO_BOTTOM ? 1. : -1.,
        tipCounter);
    fireAfterTreeEvent_(treeEvent);
  }
}

void CladogramPlot::recursivePlot_(GraphicDevice& gDevice, INode& node, double x, double& y, double hDirection, double vDirection, unsigned int& tipCounter) const
{
  double depth = static_cast<double>(TreeTemplateTools::getDepth(node));
  double x2 = ((getHorizontalOrientation() == ORIENTATION_LEFT_TO_RIGHT ? totalDepth_ : 0) - depth) * getXUnit() * hDirection;
  unique_ptr<Cursor> cursor;
  unique_ptr<DrawINodeEvent> nodeEvent;
  unique_ptr<DrawIBranchEvent> branchEvent;
  short hpos = (getHorizontalOrientation() == ORIENTATION_LEFT_TO_RIGHT ? GraphicDevice::TEXT_HORIZONTAL_LEFT : GraphicDevice::TEXT_HORIZONTAL_RIGHT);
  if (node.isLeaf())
  {
    y = ((getVerticalOrientation() == ORIENTATION_TOP_TO_BOTTOM ? 0 : getHeight()) + static_cast<double>(tipCounter) * vDirection) * getYUnit();
    tipCounter++;
    cursor.reset(new Cursor(x2, y, 0, hpos));
    nodeEvent.reset(new DrawINodeEvent(this, &gDevice, &node, *cursor));
    fireBeforeNodeEvent_(*nodeEvent);
  }
  else if (node.getInfos().isCollapsed())
  {
    y = ((getVerticalOrientation() == ORIENTATION_TOP_TO_BOTTOM ? 0 : getHeight()) + static_cast<double>(tipCounter) * vDirection) * getYUnit();
    tipCounter++;
    cursor.reset(new Cursor(x2, y, 0, hpos));
    nodeEvent.reset(new DrawINodeEvent(this, &gDevice, &node, *cursor));
    fireBeforeNodeEvent_(*nodeEvent);
  }
  else
  {
    // Vertical line. Call the method on son nodes first:
    double miny = 1000000; // (unsigned int)(-log(0));
    double maxy = 0;
    for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
    {
      double yson;
      recursivePlot_(gDevice, *node.getSon(i), x2, yson, hDirection, vDirection, tipCounter);
      if (yson < miny)
        miny = yson;
      if (yson > maxy)
        maxy = yson;
    }
    y = (maxy + miny) / 2.;
    cursor.reset(new Cursor(x2, y, 0, hpos));
    nodeEvent.reset(new DrawINodeEvent(this, &gDevice, &node, *cursor));
    fireBeforeNodeEvent_(*nodeEvent);
    gDevice.drawLine(x2, miny, x2, maxy);
  }
  // Actualize node infos:
  node.getInfos().setX(x2);
  node.getInfos().setY(y);
  nodeEvent.reset(new DrawINodeEvent(this, &gDevice, &node, *cursor));
  fireAfterNodeEvent_(*nodeEvent);

  // Horizontal line
  branchEvent.reset(new CladogramDrawBranchEvent(this, &gDevice, &node, x2 - x, *cursor, getHorizontalOrientation()));
  fireBeforeBranchEvent_(*branchEvent);
  gDevice.drawLine(x, y, x2, y);
  fireAfterBranchEvent_(*branchEvent);
}
