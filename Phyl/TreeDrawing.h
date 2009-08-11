//
// File: TreeDrawing.h
// Created by: Julien Dutheil
// Created on: Sun Oct 8 11:57 2006
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

#ifndef _TREEDRAWING_H_
#define _TREEDRAWING_H_

// From Utils:
#include <Utils/GraphicDevice.h>
#include <Utils/Point2D.h>

// From PhylLib:
#include "Tree.h"

namespace bpp
{

/**
 * @brief A set of options to tune the display of a TreeDrawing object.
 */
class TreeDrawingSettings
{
  public:
    bool drawLeafNames;
    //More options will be added in the future...
};

/**
 * @brief Basal interface for tree drawing classes.
 *
 * Basicly, a TreeDrawing object draw a graph of the tree and compute the coordinates
 * of each node on the graph. These coordinates may be retrieved by dedicated functions.
 * The drawing is performed on a GraphicDevice object.
 *
 * The TreeDrwing class is in charge of the tree reprensentation, and offer tools to retireve
 * the coordinates of nodes. Using these funcitons to plot annotation may turn to be unefficient
 * however, particularly for large trees, as they involve a search on the whole tree. For easier
 * tuning of the drawing extensions, the interface defines the drawProperty,
 * getSupportedDrawableProperties and isDrawable methods. These methods can be used to add features
 * to the plot. Adding new features can then be performed by subclassing an existing algorithm
 * and adding support for more properties.
 *
 * The TreeDrawing interface do not implies that the implementation works on a copy of the tree.
 * It takes a constant pointer toward the tree to plot. Depending on the implementation however,
 * the inheriting class may chose to store a copy of the tree for convenience. Refer to the
 * documentation of the specific implementation you are using for details.
 *
 */
class TreeDrawing
{
  public:
    TreeDrawing() {}
    virtual ~TreeDrawing() {}
  
  public:
    /**
     * @return A string describing this drawing algorithm.
     */
    virtual string getName() const = 0;

    /**
     * @return 'True' if a tree is attached to this instance.
     */
    virtual bool hasTree() const = 0;

    /**
     * @return A pointer toward the tree associated to this drawing.
     */
    virtual const Tree* getTree() const = 0;
    
    /**
     * @param A pointer toward the tree to associate with this drawing.
     */
    virtual void setTree(const Tree* tree) = 0;

    /**
     * @brief Set the 'horizontal' expansion unit.
     *
     * The effect of this expansion factor depends on the implementation of the interface.
     * @param xu The horizontal unit length.
     */
    virtual void setXUnit(double xu) = 0; 

    /**
     * @brief Set the 'vertical' expansion unit.
     *
     * The effect of this expansion factor depends on the implementation of the interface.
     * @param yu The vertical unit length.
     */
    virtual void setYUnit(double yu) = 0;

    /**
     * @return The horizontal unit length.
     */
    virtual double getXUnit() const = 0;
    
    /**
     * @return The vertical unit length.
     */
    virtual double getYUnit() const = 0;

    /**
     * @return The total width of the drawing, in X units.
     */
    virtual double getWidth() const = 0; 

    /**
     * @return The total height of the drawing, in Y units.
     */
    virtual double getHeight() const = 0; 

    /**
     * @brief Plot the tree onto the specified device.
     *
     * @param gDevice An object implementing the GraphicDevice interface.
     */
    virtual void plot(GraphicDevice& gDevice) const throw (Exception) = 0;

    /**
     * @brief Get the position of a node.
     *
     * @param nodeId The identifier of the node.
     * @return The localization of the node using the coordinate system of the last GraphicDevice used.
     * @throw NodeNotFoundException If the node does not correspond to a node in the tree.
     */
    virtual Point2D<double> getNodePosition(int nodeId) const throw (NodeNotFoundException) = 0;

    /**
     * @brief Get the node corresponding to a position on the device.
     *
     * @param position A position in the coordinates system of the last GraphicDevice used.
     * @return The corresponding node identifier if available.
     * @throw NodeNotFoundException If the node does not correspond to a node in the tree.
     */
    virtual int getNodeAt(const Point2D<double>& position) const throw (NodeNotFoundException) = 0;

    /**
     * @brief Plot a property on the tree.
     *
     * @param property The name of the property to plot.
     * @return True is the property could be drawn properly.
     */
    virtual bool drawProperty(GraphicDevice& gDevice, const string& property) const = 0;

    /**
     * @return The list of supported drawable properties by the plotting algorithm.
     */
    virtual const vector<string>& getSupportedDrawableProperties() const = 0;

    /**
     * @param property The name of the property to test.
     * @return True if the given property can be plotted on the tree.
     */
    virtual bool isDrawable(const string& property) const = 0;

    virtual void setDisplaySettings(TreeDrawingSettings& tds) = 0;
    virtual TreeDrawingSettings& getDisplaySettings() = 0;
    virtual const TreeDrawingSettings& getDisplaySettings() const = 0;
     
};

} //end of namespace bpp.

#endif //_TREEDRAWING_H_

