// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_DISTANCE_PGMA_H
#define BPP_PHYL_DISTANCE_PGMA_H


#include "../Tree/Tree.h"
#include "../Tree/TreeTemplate.h"
#include "AbstractAgglomerativeDistanceMethod.h"

namespace bpp
{
/**
 * @brief Inner data structure for WPGMA and UPGMA distance methods.
 */
struct PGMAInfos
{
  size_t numberOfLeaves;
  double time;
};

/**
 * @brief Compute WPGMA and UPGMA trees from a distance matrix.
 *
 * WPGMA = Weighted pair group method using arithmetic averaging,
 * is equivalent to the average linkage hierarchical clustering method.
 * The distance between two taxa is the average distance between all individuals in each taxa.
 * The unweighted version (named UPGMA), uses a weighted average, with the number of individuals in a group as a weight.
 */
class PGMA :
  public AbstractAgglomerativeDistanceMethod
{
protected:
  bool weighted_;

public:
  PGMA(bool weighted = true) :
    AbstractAgglomerativeDistanceMethod(true, true),
    weighted_(weighted) {}

  /**
   * @brief Create a (U/W)PGMA object instance.
   *
   * @param matrix Input distance matrix.
   * @param weighted Tell if we should perform Weighted or Unweighted pair group method.
   * @param verbose Allow to display extra information, like progress bars.
   */
  PGMA(const DistanceMatrix& matrix, bool weighted = true, bool verbose = true) :
    AbstractAgglomerativeDistanceMethod(matrix, verbose, true),
    weighted_(weighted)
  {
    computeTree();
  }
  virtual ~PGMA() {}

  PGMA* clone() const { return new PGMA(*this); }

public:
  std::string getName() const { return std::string(weighted_ ? "W" : "U") + "PGMA"; }

  void setDistanceMatrix(const DistanceMatrix& matrix)
  {
    AbstractAgglomerativeDistanceMethod::setDistanceMatrix(matrix);
  }

  void setWeighted(bool weighted) { weighted_ = weighted; }
  bool isWeighted() const { return weighted_; }

protected:
  std::vector<size_t> getBestPair();
  std::vector<double> computeBranchLengthsForPair(const std::vector<size_t>& pair);
  double computeDistancesFromPair(const std::vector<size_t>& pair, const std::vector<double>& branchLengths, size_t pos);
  void finalStep(int idRoot);
  virtual Node* getLeafNode(int id, const std::string& name);
  virtual Node* getParentNode(int id, Node* son1, Node* son2);
};
} // end of namespace bpp.
#endif // BPP_PHYL_DISTANCE_PGMA_H
