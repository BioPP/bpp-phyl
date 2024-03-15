// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractDendrogramPlot.h"

using namespace bpp;
using namespace std;

short AbstractDendrogramPlot::ORIENTATION_LEFT_TO_RIGHT = 1;
short AbstractDendrogramPlot::ORIENTATION_RIGHT_TO_LEFT = 2;
short AbstractDendrogramPlot::ORIENTATION_TOP_TO_BOTTOM = 3;
short AbstractDendrogramPlot::ORIENTATION_BOTTOM_TO_TOP = 4;

void AbstractDendrogramPlot::plot(GraphicDevice& gDevice) const
{
  gDevice.setCurrentPointSize(getDisplaySettings().pointSize);
  drawDendrogram_(gDevice);
}
