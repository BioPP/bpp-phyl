// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_GRAPHICS_PHYLOGRAMPLOT_H
#define BPP_PHYL_GRAPHICS_PHYLOGRAMPLOT_H


#include "AbstractDendrogramPlot.h"

namespace bpp
{
class PhylogramDrawBranchEvent :
  public DrawIBranchEvent
{
private:
  double orientation_;

public:
  PhylogramDrawBranchEvent(const TreeDrawing* source, GraphicDevice* gd, const INode* node, const Cursor& cursor, short orientation);

public:
  Cursor getBranchCursor(double position) const;
};


/**
 * @brief Phylogram representation of trees.
 *
 * This representation is for trees with branch lengths.
 */
class PhylogramPlot :
  public AbstractDendrogramPlot
{
private:
  double totalDepth_;
  double numberOfLeaves_;

public:
  PhylogramPlot() :
    AbstractDendrogramPlot(), totalDepth_(0), numberOfLeaves_(0)
  {}

  virtual ~PhylogramPlot() {}

  PhylogramPlot* clone() const override { return new PhylogramPlot(*this); }

public:
  std::string getName() const override { return "Phylogram"; }

  void setTree(const Tree& tree) override;

  double getWidth() const override { return totalDepth_; }
  double getHeight() const override { return numberOfLeaves_; }

  void treeHasChanged() override
  {
    if (hasTree())
    {
      tree_->setVoidBranchLengths(0.);
      totalDepth_ = TreeTemplateTools::getHeight(*tree_->getRootNode());
      numberOfLeaves_ = static_cast<double>(tree_->getNumberOfLeaves());
    }
  }

private:
  void drawDendrogram_(GraphicDevice& gDevice) const override;

  void recursivePlot_(GraphicDevice& gDevice, INode& node, double x, double& y, double hDirection, double vDirection, unsigned int& tipCounter) const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_GRAPHICS_PHYLOGRAMPLOT_H
