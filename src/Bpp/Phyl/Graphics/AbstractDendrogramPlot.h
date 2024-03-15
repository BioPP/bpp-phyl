// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_GRAPHICS_ABSTRACTDENDROGRAMPLOT_H
#define BPP_PHYL_GRAPHICS_ABSTRACTDENDROGRAMPLOT_H


#include "AbstractTreeDrawing.h"

namespace bpp
{
/**
 * @brief Basic implementation of dendrogram plots.
 *
 * Dendrograms are oriented plots, with all the leaves on one side of the plot, and the root node at the opposite side.
 * This implementation offers to option for ploting form left to right or right to left. This will affect the direction
 * of plot annotations. The drawing can always be transformed using the regular translation/rotation operation on the
 * GraphicDevice.
 */
class AbstractDendrogramPlot :
  public AbstractTreeDrawing
{
private:
  short horOrientation_;
  short verOrientation_;

public:
  AbstractDendrogramPlot() :
    AbstractTreeDrawing(), horOrientation_(ORIENTATION_LEFT_TO_RIGHT), verOrientation_(ORIENTATION_TOP_TO_BOTTOM)
  {}

public:
  void setHorizontalOrientation(short orientation) { horOrientation_ = orientation; }
  void setVerticalOrientation(short orientation) { verOrientation_ = orientation; }

  short getHorizontalOrientation() const { return horOrientation_; }
  short getVerticalOrientation() const { return verOrientation_; }

  void plot(GraphicDevice& gDevice) const;

protected:
  virtual void drawDendrogram_(GraphicDevice& gDevice) const = 0;

public:
  static short ORIENTATION_LEFT_TO_RIGHT;
  static short ORIENTATION_RIGHT_TO_LEFT;
  static short ORIENTATION_TOP_TO_BOTTOM;
  static short ORIENTATION_BOTTOM_TO_TOP;
};
} // end of namespace bpp;
#endif // BPP_PHYL_GRAPHICS_ABSTRACTDENDROGRAMPLOT_H
