// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_TREE_BIPARTITIONLIST_H
#define BPP_PHYL_TREE_BIPARTITIONLIST_H

#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Utils/MapTools.h>

#include "Tree.h"

// From the STL:
#include <map>
#include <algorithm>

namespace bpp
{
class Node;
template<class N> class TreeTemplate;

/**
 * @brief This class deals with the bipartitions defined by trees
 *
 * Any branch of a tree defines a bipartition, i.e. two non-overlapping subsets of leaves
 * whose union is the whole set of leaves. A tree topology can therefore be represented by
 * the list of bipartitions its branches define - the so-called matrix representation.
 * Coding trees this way is useful for comparing topologies, calculating topological distances,
 * producing consensus trees or super-trees, calculating bootstrap support.
 *
 * A BipartitionList includes a set of element names (typically leaf names) and a vector of arrays of bits
 * (int*, not bitsets, are used to allow for dynamic memory allocation). Each array of bit codes for one bipartition.
 * Each bit in an array of bit corresponds to one element, so the order of element names matter.
 * Bits set to zero versus bits set to one define the two partitions of elements.
 * A BipartitionList is called sorted if its elements (leaf names) are in alphabetic order (recommended).
 *
 * BipartitionList objects are typically created from a tree, in which case elements are leaf names.
 * Note that BipartitionList is an unrooted object: a rooted or unrooted versions of the same tree will
 * yield the same BipartitionList in which the root location is ignored. Bipartitions can be accessed as
 * arrays of bits (e.g. getBitBipartition), or as map<string, bool>, in which keys are leaf names and
 * true/false values define the two partitions (e.g. getBipartition, addBipartition).
 *
 * @see Tree
 * @see BipartitionTools
 * @see TreeTools
 */
class BipartitionList :
  public virtual Clonable
{
private:
  std::vector<int*> bitBipartitionList_;
  std::vector<std::string> elements_;
  bool sorted_;

public:
  /**
   * @brief The main constructor
   *
   * @param tr The tree to be coded as bipartitions
   * @param sorted Tells whether leave names should be alphabetically sorted (recommended)
   * @param index An output optional vector to keep trace of the nodes id underlying each bipartition.
   */
  BipartitionList(const Tree& tr, bool sorted = true, std::vector<int>* index = 0);

  /**
   * @brief An alternative constructor in which elements and bipartitions are passed directly
   *
   * @param elements Leaf names
   * @param bipl The list of bit-encoded bipartitions
   */
  BipartitionList(const std::vector<std::string>& elements, const std::vector<int*>& bipl);

  /**
   * @brief Copy-constructor
   */
  BipartitionList(const BipartitionList& bipl);

  /**
   * @brief Assignment operator
   */
  BipartitionList& operator=(const BipartitionList& bipl);

  virtual ~BipartitionList();

  BipartitionList* clone() const { return new BipartitionList(*this); }

public:
  size_t getNumberOfElements() const { return elements_.size(); }

  const std::vector<std::string>& getElementNames() const { return elements_; }

  size_t getNumberOfBipartitions() const { return bitBipartitionList_.size(); }

  const std::vector<int*>& getBitBipartitionList() const { return bitBipartitionList_; }

  std::map<std::string, bool> getBipartition(size_t i) const;

  int* getBitBipartition(size_t i);

  bool haveSameElementsThan(std::map<std::string, bool>& bipart) const;

  void addBipartition(std::map<std::string, bool>& bipart, bool checkElements = 1);

  void deleteBipartition(size_t i);

  bool isSorted() const { return sorted_; }

  void sortElements();

  bool containsBipartition(std::map<std::string, bool>& bipart, bool checkElements = 1) const;

  bool areIdentical(size_t k1, size_t k2) const;

  void removeRedundantBipartitions();

  /**
   * @brief Tells whether 2 bipartitions from the list are compatible
   *
   * Let A=A1|A2 and B=B1|B2 be two bipartitions (such that A1 U A2 = B1 U B2 = set of leaves)
   * A and B are said compatible if (A1 contains B1 and B2 contains A2) or (A1 contains B2
   * and B1 contains A2) or (B1 contains A1 and A2 contains B2) or (B2 contains A1 and A2
   * contains B1). Only compatible bipartitions can belong to the same tree.
   */
  bool areCompatible(size_t k1, size_t k2) const;

  /**
   * @brief Tells whether all bipartitions in the list are compatible with each other
   */
  bool areAllCompatible() const;

  /**
   * @brief Tells whether all bipartitions in the list are compatible with a given bipartition
   *
   * @param bipart A map representing a bipartition.
   * @param checkElements Check the correspondence of element sets or not.
   */
  bool areAllCompatibleWith(std::map<std::string, bool>& bipart, bool checkElements = true) const;

  /**
   * @brief Removes bipartitions corresponding to external branches (1 vs n-1)
   */
  void removeTrivialBipartitions();

  /**
   * @brief Adds bipartitions corresponding to external branches if missing
   */
  void addTrivialBipartitions(bool checkExisting);


  /**
   * @brief Replaces ones by zeros and zeros by ones in the ith bipartition
   */
  void flip(size_t i);

  /**
   * @brief Returns the size of the smallest of the two partitions (e.g. 1 for external branches)
   */
  size_t getPartitionSize(size_t i) const;

  /**
   * @brief Sort bipartitions by partition size
   */
  void sortByPartitionSize();

  /**
   * @brief Translate into a tree
   */
  std::unique_ptr<TreeTemplate<Node>> toTree() const;

  /**
   * @brief Create a matrix representation of the bifurcations.
   *
   * Each row corresponds to an element, each column to a bipartition.
   *
   * NB: using RowMatrix<bool> leads to unexplained compilation error...
   *
   * @return A boolean matrix.
   */
  RowMatrix<int> toMatrix() const;

private:
  std::vector<std::string> buildBitBipartitions(const Node* nd, std::vector<int*>& bitbip, const std::vector<std::string>& elements, size_t* cpt, std::vector<int>* index) const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_TREE_BIPARTITIONLIST_H
