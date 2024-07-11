// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_GRAPHICS_CLADOGRAMPLOT_H
#define BPP_PHYL_GRAPHICS_CLADOGRAMPLOT_H

#include <Bpp/Exceptions.h>

#include "AbstractDendrogramPlot.h"

namespace bpp
{
class CladogramDrawBranchEvent :
  public DrawIBranchEvent
{
private:
  double orientation_;
  double length_;

public:
  CladogramDrawBranchEvent(const TreeDrawing* source, GraphicDevice* gd, const INode* node, double length_, const Cursor& cursor, short orientation);

public:
  Cursor getBranchCursor(double position) const;
};


/**
 * @brief Cladogram representation of trees.
 *
 * This representation is for trees without branch lengths.
 */
class CladogramPlot :
  public AbstractDendrogramPlot
{
private:
  double totalDepth_;
  double numberOfLeaves_;

public:
  CladogramPlot() :
    AbstractDendrogramPlot(), totalDepth_(0), numberOfLeaves_(0)
  {}

  virtual ~CladogramPlot() {}

  CladogramPlot* clone() const { return new CladogramPlot(*this); }

public:
  std::string getName() const { return "Cladogram"; }

  void setTree(const Tree* tree = 0);

  double getWidth() const { return totalDepth_; }
  double getHeight() const { return numberOfLeaves_; }

  void treeHasChanged()
  {
    if (hasTree())
    {
      totalDepth_ = static_cast<double>(TreeTemplateTools::getDepth(*getTree_()->getRootNode()));
      numberOfLeaves_ = static_cast<double>(getTree_()->getNumberOfLeaves());
    }
  }

private:
  void drawDendrogram_(GraphicDevice& gDevice) const;
  void recursivePlot_(GraphicDevice& gDevice, INode& node, double x, double& y, double hDirection, double vDirection, unsigned int& tipCounter) const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_GRAPHICS_CLADOGRAMPLOT_H
