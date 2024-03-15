// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_DISTANCE_NEIGHBORJOINING_H
#define BPP_PHYL_DISTANCE_NEIGHBORJOINING_H


#include "AbstractAgglomerativeDistanceMethod.h"

namespace bpp
{
/**
 * @brief The neighbor joining distance method.
 *
 * Reference:
 * N Saitou and M Nei (1987), _Molecular Biology and Evolution_ 4(4) 406-25.
 */
class NeighborJoining :
  public AbstractAgglomerativeDistanceMethod
{
protected:
  std::vector<double> sumDist_;
  bool positiveLengths_;

public:
  /**
   * @brief Create a new NeighborJoining object instance, without performing any computation.
   *
   * @param rooted Tell if the output tree should be rooted.
   * @param positiveLengths Tell if negative lengths should be avoided.
   * @param verbose Allow to display extra information, like progress bars.
   */
  NeighborJoining(bool rooted = false, bool positiveLengths = false, bool verbose = true) :
    AbstractAgglomerativeDistanceMethod(verbose, rooted),
    sumDist_(),
    positiveLengths_(false)
  {}

  /**
   * @brief Create a new NeighborJoining object instance and compute a tree from a distance matrix.
   *
   * @param matrix Input distance matrix.
   * @param rooted Tell if the output tree should be rooted.
   * @param positiveLengths Tell if negative lengths should be avoided.
   * @param verbose Allow to display extra information, like progress bars.
   */
  NeighborJoining(const DistanceMatrix& matrix, bool rooted = false, bool positiveLengths = false, bool verbose = true) :
    AbstractAgglomerativeDistanceMethod(matrix, verbose, rooted),
    sumDist_(),
    positiveLengths_(positiveLengths)
  {
    sumDist_.resize(matrix.size());
    computeTree();
  }

  virtual ~NeighborJoining() {}

  NeighborJoining* clone() const { return new NeighborJoining(*this); }

public:
  std::string getName() const { return "NJ"; }

  virtual void setDistanceMatrix(const DistanceMatrix& matrix)
  {
    AbstractAgglomerativeDistanceMethod::setDistanceMatrix(matrix);
    sumDist_.resize(matrix.size());
  }

  virtual void outputPositiveLengths(bool yn) { positiveLengths_ = yn; }

protected:
  std::vector<size_t> getBestPair();
  std::vector<double> computeBranchLengthsForPair(const std::vector<size_t>& pair);
  double computeDistancesFromPair(const std::vector<size_t>& pair, const std::vector<double>& branchLengths, size_t pos);
  void finalStep(int idRoot);
};
} // end of namespace bpp.
#endif // BPP_PHYL_DISTANCE_NEIGHBORJOINING_H
