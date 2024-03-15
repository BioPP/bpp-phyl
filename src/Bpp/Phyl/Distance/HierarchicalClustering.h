// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_DISTANCE_HIERARCHICALCLUSTERING_H
#define BPP_PHYL_DISTANCE_HIERARCHICALCLUSTERING_H


#include "AbstractAgglomerativeDistanceMethod.h"

namespace bpp
{
class ClusterInfos
{
public:
  size_t numberOfLeaves;
  double length;

public:
  ClusterInfos() : numberOfLeaves(0),
    length(0) {}
};

/**
 * @brief Hierarchical clustering.
 *
 * This class implements the complete, single, average (= UPGMA), median, ward and centroid linkage methods.
 */
class HierarchicalClustering :
  public AbstractAgglomerativeDistanceMethod
{
public:
  static const std::string COMPLETE;
  static const std::string SINGLE;
  static const std::string AVERAGE;
  static const std::string MEDIAN;
  static const std::string WARD;
  static const std::string CENTROID;

protected:
  std::string method_;

public:
  /**
   * @brief Builds a new clustering object.
   *
   * @param method The linkage method to use. should be one of COMPLETE, SINGLE, AVERAGE, MEDIAN, WARD, CENTROID.
   * @param verbose Tell if some progress information should be displayed.
   */
  HierarchicalClustering(const std::string& method, bool verbose = false) :
    AbstractAgglomerativeDistanceMethod(verbose),
    method_(method) {}
  HierarchicalClustering(const std::string& method, const DistanceMatrix& matrix, bool verbose = false) :
    AbstractAgglomerativeDistanceMethod(matrix, verbose, true),
    method_(method)
  {
    computeTree();
  }

  virtual ~HierarchicalClustering() {}

  HierarchicalClustering* clone() const { return new HierarchicalClustering(*this); }

public:
  std::string getName() const { return "Hierarchical clustering: " + method_; }

protected:
  std::vector<size_t> getBestPair();
  std::vector<double> computeBranchLengthsForPair(const std::vector<size_t>& pair);
  double computeDistancesFromPair(const std::vector<size_t>& pair, const std::vector<double>& branchLengths, size_t pos);
  void finalStep(int idRoot);
  virtual Node* getLeafNode(int id, const std::string& name);
  virtual Node* getParentNode(int id, Node* son1, Node* son2);
};
} // end of namespace bpp.
#endif // BPP_PHYL_DISTANCE_HIERARCHICALCLUSTERING_H
