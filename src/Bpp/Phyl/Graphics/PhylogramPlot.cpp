//
// File: PhylogramPlot.cpp
// Created by: Julien Dutheil
// Created on: Tue Oct 3 20:52 2006
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

#include "PhylogramPlot.h"
#include "../Tree/TreeTemplateTools.h"

using namespace bpp;

#include <string>

using namespace std;

PhylogramDrawBranchEvent::PhylogramDrawBranchEvent(const TreeDrawing* source, GraphicDevice* gd, const INode* node, const Cursor& cursor, short orientation) :
  DrawIBranchEvent(source, gd, node, cursor), orientation_()
{
  orientation_ = (orientation == AbstractDendrogramPlot::ORIENTATION_LEFT_TO_RIGHT ? -1 : 1);
}

Cursor PhylogramDrawBranchEvent::getBranchCursor(double position) const
{
  double offset = 0;
  if (getINode()->hasDistanceToFather())
  {
    double l = getINode()->getDistanceToFather();
    offset = orientation_ * l * position * getTreeDrawing()->getXUnit();
  }
  return getCursor().getTranslation(offset, 0);
}

void PhylogramPlot::setTree(const Tree* tree)
{
  AbstractDendrogramPlot::setTree(tree);
  if (tree)
  {
    getTree_()->setVoidBranchLengths(0.);
  }
}

void PhylogramPlot::drawDendrogram_(GraphicDevice& gDevice) const throw (Exception)
{
  if (hasTree())
  {
    DrawTreeEvent treeEvent(this, &gDevice);
    fireBeforeTreeEvent_(treeEvent);
    unsigned int* tipCounter = new unsigned int(0);
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

void PhylogramPlot::recursivePlot_(GraphicDevice& gDevice, INode& node, double x, double& y, double hDirection, double vDirection, unsigned int* tipCounter) const
{
  double x2;
  bool drawBranch = true;
  if (node.hasDistanceToFather())
  {
    double length = node.hasDistanceToFather() ? node.getDistanceToFather() : 0.;
    if (length < -10000000)
    {
      x2 = x;
      drawBranch = false;
    }
    else
    {
      x2 = x + hDirection * length * getXUnit();
    }
  }
  else
  {
    x2 = x;
    drawBranch = false;
  }
 
  auto_ptr<Cursor> cursor;
  auto_ptr<DrawINodeEvent> nodeEvent;
  auto_ptr<DrawIBranchEvent> branchEvent;
  short hpos = (getHorizontalOrientation() == ORIENTATION_LEFT_TO_RIGHT ? GraphicDevice::TEXT_HORIZONTAL_LEFT : GraphicDevice::TEXT_HORIZONTAL_RIGHT);
  if (node.isLeaf())
  {
    y = ((getVerticalOrientation() == ORIENTATION_TOP_TO_BOTTOM ? 0 : getHeight()) + static_cast<double>(*tipCounter) * vDirection) * getYUnit();
    (*tipCounter)++;
    cursor.reset(new Cursor(x2, y, 0, hpos));
    nodeEvent.reset(new DrawINodeEvent(this, &gDevice, &node, *cursor));
    fireBeforeNodeEvent_(*nodeEvent);
  }
  else if (node.getInfos().isCollapsed())
  {
    y = ((getVerticalOrientation() == ORIENTATION_TOP_TO_BOTTOM ? 0 : getHeight()) + static_cast<double>(*tipCounter) * vDirection) * getYUnit();
    (*tipCounter)++;
    cursor.reset(new Cursor(x2, y, 0, hpos));
    nodeEvent.reset(new DrawINodeEvent(this, &gDevice, &node, *cursor));
    fireBeforeNodeEvent_(*nodeEvent);
  }
  else
  {
    //Vertical line. Call the method on son nodes first:
    double miny = 1000000; //(unsigned int)(-log(0));
    double maxy = 0;
    for (unsigned int i = 0; i < node.getNumberOfSons(); i++)
    {
      double yson;
      recursivePlot_(gDevice, *node.getSon(i), x2, yson, hDirection, vDirection, tipCounter);
      if(yson < miny) miny = yson;
      if(yson > maxy) maxy = yson;
    }
    y = (maxy + miny) / 2.;
    cursor.reset(new Cursor(x2, y, 0, hpos));
    nodeEvent.reset(new DrawINodeEvent(this, &gDevice, &node, *cursor));
    fireBeforeNodeEvent_(*nodeEvent);
    gDevice.drawLine(x2, miny, x2, maxy);
  }
  
  //Actualize node infos:
  node.getInfos().setX(x2);
  node.getInfos().setY(y);
  nodeEvent.reset(new DrawINodeEvent(this, &gDevice, &node, *cursor));
  fireAfterNodeEvent_(*nodeEvent);

  if (drawBranch)
  {
    //Horizontal line
    branchEvent.reset(new PhylogramDrawBranchEvent(this, &gDevice, &node, *cursor, getHorizontalOrientation()));
    fireBeforeBranchEvent_(*branchEvent);
    gDevice.drawLine(x, y, x2, y);
    fireAfterBranchEvent_(*branchEvent);
  }
}


