// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractTreeDrawing.h"
#include "TreeDrawingListener.h"

// From the STL:
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
    // Pointer-based event (efficient):
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
    // Id-based event (less-efficient):
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
    // Pointer-based event (efficient):
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
    // Id-based event (less-efficient):
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
    // Pointer-based event (efficient):
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
                   TextTools::toString(dynamic_cast<const Number<double>*>(b)->getValue()),
                   cursor.getHPos(), GraphicDevice::TEXT_VERTICAL_CENTER, cursor.getAngle());
      gd->setCurrentFont(fontBck);
    }
  }
  catch (std::bad_cast& e)
  {
    // Id-based event (less-efficient):
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
                   TextTools::toString(dynamic_cast<const Number<double>*>(b)->getValue()),
                   cursor.getHPos(), GraphicDevice::TEXT_VERTICAL_CENTER, cursor.getAngle());
      gd->setCurrentFont(fontBck);
    }
  }
}

void LabelInnerNodesTreeDrawingListener::afterDrawNode(const DrawNodeEvent& event)
{
  try
  {
    // Pointer-based event (efficient):
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
  catch (std::bad_cast& e)
  {
    // Id-based event (less-efficient):
    const TreeDrawing* td = event.getTreeDrawing();
    if (!td->isNodeCollapsed(event.getNodeId()))
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
    // Pointer-based event (efficient):
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
  catch (std::bad_cast& e)
  {
    // Id-based event (less-efficient):
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
