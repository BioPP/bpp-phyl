// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_DISTANCE_BIONJ_H
#define BPP_PHYL_DISTANCE_BIONJ_H


#include "NeighborJoining.h"

namespace bpp
{
/**
 * @brief The BioNJ distance method.
 *
 * Reference:
 * Gascuel O.
 * BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data.
 * Mol Biol Evol. 1997 Jul;14(7):685-95.
 */
class BioNJ :
  public NeighborJoining
{
private:
  DistanceMatrix variance_;
  double lambda_;

public:
  /**
   * @brief Create a new BioNJ object instance and compute a tree from a distance matrix.
   *
   * @param rooted Tell if the output tree should be rooted.
   * @param positiveLengths Tell if negative lengths should be avoided.
   * @param verbose Allow to display extra information, like progress bars.
   */
  BioNJ(bool rooted = false, bool positiveLengths = false, bool verbose = true) :
    NeighborJoining(rooted, positiveLengths, verbose),
    variance_(0),
    lambda_(0) {}

  /**
   * @brief Create a new BioNJ object instance and compute a tree from a distance matrix.
   *
   * @param matrix Input distance matrix.
   * @param rooted Tell if the output tree should be rooted.
   * @param positiveLengths Tell if negative lengths should be avoided.
   * @param verbose Allow to display extra information, like progress bars.
   */
  BioNJ(const DistanceMatrix& matrix, bool rooted = false, bool positiveLengths = false, bool verbose = true) :
    NeighborJoining(rooted, positiveLengths, verbose),
    // Use the default constructor, because the other one call computeTree.
    variance_(matrix),
    lambda_(0)
  {
    setDistanceMatrix(matrix);
    outputPositiveLengths(positiveLengths);
    computeTree();
  }

  BioNJ* clone() const { return new BioNJ(*this); }

  virtual ~BioNJ() {}

public:
  std::string getName() const { return "BioNJ"; }

  void setDistanceMatrix(const DistanceMatrix& matrix)
  {
    NeighborJoining::setDistanceMatrix(matrix);
    variance_ = matrix;
  }
  void computeTree();
  double computeDistancesFromPair(const std::vector<size_t>& pair, const std::vector<double>& branchLengths, size_t pos);
};
} // end of namespace bpp.
#endif // BPP_PHYL_DISTANCE_BIONJ_H
