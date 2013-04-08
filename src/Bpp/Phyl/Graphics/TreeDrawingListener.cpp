//
// File: TreeDrawingListener.cpp
// Created by: Julien Dutheil
// Created on: Tue May 18 10:33 2010
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

#include "TreeDrawingListener.h"
#include "AbstractTreeDrawing.h"

//From the STL:
#include <string>
#include <exception>
#include <typeinfo>

using namespace std;

using namespace bpp;

void NodesIdTreeDrawingListener::afterDrawNode(const DrawNodeEvent& event)
{
  GraphicDevice* gd = event.getGraphicDevice();
  Cursor cursor     = event.getCursor();
  Font fontBck      = gd->getCurrentFont();
  if (settings_)
    gd->setCurrentFont(settings_->fontNodesId);
  string name = "#" + TextTools::toString(event.getNodeId());
  gd->drawText(cursor.getX(), cursor.getY(), name, cursor.getHPos(), cursor.getVPos(), cursor.getAngle());
  gd->setCurrentFont(fontBck);
}

void LeafNamesTreeDrawingListener::afterDrawNode(const DrawNodeEvent& event)
{
  try
  {
    //Pointer-based event (efficient):
    const DrawINodeEvent& eventC = dynamic_cast<const DrawINodeEvent&>(event);
    if (eventC.getINode()->isLeaf())
    {
      GraphicDevice* gd = event.getGraphicDevice();
      Cursor cursor     = event.getCursor();
      Font fontBck      = gd->getCurrentFont();
      if (settings_)
        gd->setCurrentFont(settings_->fontLeafNames);
      string name = eventC.getINode()->getName();
      gd->drawText(cursor.getX(), cursor.getY(), name, cursor.getHPos(), cursor.getVPos(), cursor.getAngle());
      gd->setCurrentFont(fontBck);
    }
  }
  catch (bad_cast& e)
  {
    //Id-based event (less-efficient):
    const TreeDrawing* td = event.getTreeDrawing();
    if (td->getTree()->isLeaf(event.getNodeId()))
    {
      GraphicDevice* gd = event.getGraphicDevice();
      Cursor cursor     = event.getCursor();
      Font fontBck      = gd->getCurrentFont();
      if (settings_)
        gd->setCurrentFont(settings_->fontLeafNames);
      string name = td->getTree()->getNodeName(event.getNodeId());
      gd->drawText(cursor.getX(), cursor.getY(), name, cursor.getHPos(), cursor.getVPos(), cursor.getAngle());
      gd->setCurrentFont(fontBck);
    }
  }
}

void BranchLengthsTreeDrawingListener::afterDrawBranch(const DrawBranchEvent& event)
{
  try
  {
    //Pointer-based event (efficient):
    const DrawINodeEvent& eventC = dynamic_cast<const DrawINodeEvent&>(event);
    if (eventC.getINode()->hasDistanceToFather())
    {
      GraphicDevice* gd = event.getGraphicDevice();
      Cursor cursor     = event.getBranchCursor(0.5);
      Font fontBck      = gd->getCurrentFont();
      if (settings_)
        gd->setCurrentFont(settings_->fontBranchLengths);
      gd->drawText(cursor.getX(), cursor.getY(),
          TextTools::toString(eventC.getINode()->getDistanceToFather()),
          GraphicDevice::TEXT_HORIZONTAL_CENTER, GraphicDevice::TEXT_VERTICAL_BOTTOM, cursor.getAngle());
      gd->setCurrentFont(fontBck);
    }
  }
  catch (std::bad_cast& e)
  {
    //Id-based event (less-efficient):
    const TreeDrawing* td = event.getTreeDrawing();
    if (td->getTree()->hasDistanceToFather(event.getNodeId()))
    {
      GraphicDevice* gd = event.getGraphicDevice();
      Cursor cursor     = event.getBranchCursor(0.5);
      Font fontBck      = gd->getCurrentFont();
      if (settings_)
        gd->setCurrentFont(settings_->fontLeafNames);
      gd->drawText(cursor.getX(), cursor.getY(),
          TextTools::toString(td->getTree()->getDistanceToFather(event.getNodeId())),
          GraphicDevice::TEXT_HORIZONTAL_CENTER, GraphicDevice::TEXT_VERTICAL_BOTTOM, cursor.getAngle());
      gd->setCurrentFont(fontBck);
    }
  }
}

void BootstrapValuesTreeDrawingListener::afterDrawBranch(const DrawBranchEvent& event)
{
  try
  {
    //Pointer-based event (efficient):
    const DrawINodeEvent& eventC = dynamic_cast<const DrawINodeEvent&>(event);
    if (eventC.getINode()->hasBranchProperty(TreeTools::BOOTSTRAP))
    {
      const Clonable* b = eventC.getINode()->getBranchProperty(TreeTools::BOOTSTRAP);
      GraphicDevice* gd = event.getGraphicDevice();
      Cursor cursor     = event.getCursor();
      Font fontBck      = gd->getCurrentFont();
      if (settings_)
        gd->setCurrentFont(settings_->fontBranchLengths);
      gd->drawText(cursor.getX(), cursor.getY(),
          TextTools::toString(dynamic_cast<const Number<double> *>(b)->getValue()),
          cursor.getHPos(), GraphicDevice::TEXT_VERTICAL_CENTER, cursor.getAngle());
      gd->setCurrentFont(fontBck);
    }
  }
  catch (std::bad_cast& e)
  {
    //Id-based event (less-efficient):
    const TreeDrawing* td = event.getTreeDrawing();
    if (td->getTree()->hasBranchProperty(event.getNodeId(), TreeTools::BOOTSTRAP))
    {
      const Clonable* b = td->getTree()->getBranchProperty(event.getNodeId(), TreeTools::BOOTSTRAP);
      GraphicDevice* gd = event.getGraphicDevice();
      Cursor cursor     = event.getCursor();
      Font fontBck      = gd->getCurrentFont();
      if (settings_)
        gd->setCurrentFont(settings_->fontLeafNames);
      gd->drawText(cursor.getX(), cursor.getY(),
          TextTools::toString(dynamic_cast<const Number<double> *>(b)->getValue()),
          cursor.getHPos(), GraphicDevice::TEXT_VERTICAL_CENTER, cursor.getAngle());
      gd->setCurrentFont(fontBck);
    }
  }
}

void LabelInnerNodesTreeDrawingListener::afterDrawNode(const DrawNodeEvent& event)
{
  try
  {
    //Pointer-based event (efficient):
    const DrawINodeEvent& eventC = dynamic_cast<const DrawINodeEvent&>(event);
    if (!eventC.getINode()->getInfos().isCollapsed())
    {
      GraphicDevice* gd = event.getGraphicDevice();
      Cursor cursor     = event.getCursor();
      if (eventC.getINode()->hasName())
      {
        string name = eventC.getINode()->getName();
        gd->drawText(cursor.getX(), cursor.getY(), name, cursor.getHPos(), cursor.getVPos(), cursor.getAngle());
      }
    }
  }
  catch(std::bad_cast& e)
  {
    //Id-based event (less-efficient):
    const TreeDrawing* td = event.getTreeDrawing();
    if (! td->isNodeCollapsed(event.getNodeId()))
    {
      GraphicDevice* gd = event.getGraphicDevice();
      Cursor cursor     = event.getCursor();
      if (td->getTree()->hasNodeName(event.getNodeId()))
      {
        string name = td->getTree()->getNodeName(event.getNodeId());
        gd->drawText(cursor.getX(), cursor.getY(), name, cursor.getHPos(), cursor.getVPos(), cursor.getAngle());
      }
    }
  }
}

void LabelCollapsedNodesTreeDrawingListener::afterDrawNode(const DrawNodeEvent& event)
{
  try
  {
    //Pointer-based event (efficient):
    const DrawINodeEvent& eventC = dynamic_cast<const DrawINodeEvent&>(event);
    if (eventC.getINode()->getInfos().isCollapsed())
    {
      GraphicDevice* gd = event.getGraphicDevice();
      Cursor cursor     = event.getCursor();
      size_t size       = TreeTemplateTools::getNumberOfLeaves(*eventC.getINode());
      string text = "";
      if (eventC.getINode()->hasName())
        text += eventC.getINode()->getName() + " ";
      text += "(" + TextTools::toString(size) + " leaves)";
      gd->drawText(cursor.getX(), cursor.getY(), text, cursor.getHPos(), cursor.getVPos(), cursor.getAngle());
    }
  }
  catch(std::bad_cast& e)
  {
    //Id-based event (less-efficient):
    const TreeDrawing* td = event.getTreeDrawing();
    if (td->isNodeCollapsed(event.getNodeId()))
    {
      GraphicDevice* gd = event.getGraphicDevice();
      Cursor cursor     = event.getCursor();
      size_t size = TreeTools::getNumberOfLeaves(*td->getTree(), event.getNodeId());
      string text = "";
      if (td->getTree()->hasNodeName(event.getNodeId()))
        text += td->getTree()->getNodeName(event.getNodeId()) + " ";
      text += "(" + TextTools::toString(size) + " leaves)";
      gd->drawText(cursor.getX(), cursor.getY(), text, cursor.getHPos(), cursor.getVPos(), cursor.getAngle());
    }
  }
}

