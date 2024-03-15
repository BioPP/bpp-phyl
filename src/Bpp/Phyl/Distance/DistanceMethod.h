// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_DISTANCE_DISTANCEMETHOD_H
#define BPP_PHYL_DISTANCE_DISTANCEMETHOD_H



// From bpp-core:
#include <Bpp/Clonable.h>
#include <Bpp/Exceptions.h>

// From the STL:
#include <string>

namespace bpp
{
class DistanceMatrix;
class Node;
class Tree;

/**
 * @brief General interface for distance-based phylogenetic reconstruction methods.
 */
class DistanceMethodInterface :
  public virtual Clonable
{
public:
  DistanceMethodInterface() {}
  virtual ~DistanceMethodInterface() {}

  virtual DistanceMethodInterface* clone() const override = 0;

public:
  /**
   * @brief Set the distance matrix to use.
   *
   * @param matrix The matrix to use.
   * @throw Exception In case an incorrect matrix is provided (eg smaller than 3).
   */
  virtual void setDistanceMatrix(const DistanceMatrix& matrix) = 0;

  /**
   * @brief Perform the clustering.
   */
  virtual void computeTree() = 0;

  /**
   * @return True if a tree has been computed.
   */
  virtual bool hasTree() const = 0;

  /**
   * @return A reference toward the computed tree. Throws an exception if no tree was computed.
   */
  virtual const Tree& tree() const = 0;

  /**
   * @return The name of the distance method.
   */
  virtual std::string getName() const = 0;

  /**
   * @param yn Enable/Disable verbose mode.
   */
  virtual void setVerbose(bool yn) = 0;

  /**
   * @return True if verbose mode is enabled.
   */
  virtual bool isVerbose() const = 0;
};

/**
 * @brief Interface for agglomerative distance methods.
 *
 * This interface does not contain any specific method and
 * is there only for "ontology" purposes. Specific methods
 * might be added later.
 */
class AgglomerativeDistanceMethodInterface :
  public virtual DistanceMethodInterface
{
public:
  AgglomerativeDistanceMethodInterface() {}
  virtual ~AgglomerativeDistanceMethodInterface() {}
  
  virtual AgglomerativeDistanceMethodInterface* clone() const override = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_DISTANCE_DISTANCEMETHOD_H
