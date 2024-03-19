// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_TREE_BIPARTITIONTOOLS_H
#define BPP_PHYL_TREE_BIPARTITIONTOOLS_H


#include "BipartitionList.h"

// From bpp-seq:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

namespace bpp
{
/**
 * @brief This class provides tools related to the BipartitionList class
 *
 * BipartitionTools includes functions dealing with dynamic arrays of bits
 * and tools for dealing with bipartitions from distinct lists.
 *
 * @see Tree
 * @see BipartitionList
 * @see TreeTools
 */
class BipartitionTools
{
public:
  /**
   * @brief
   *
   * Unit length (in bits) of arrays of bits. Must be a multiple of CHAR_BIT*sizeof(int).
   * Default value is 64.
   */
  static int LWORD;

public:
  BipartitionTools() {}
  virtual ~BipartitionTools() {}

public:
  /**
   * @brief Sets bit number num of bit array list to one
   *
   * Note that no control of memory allocation is made
   */
  static void bit1(int* list, int num);

  /**
   * @brief Sets bit number num of bit array plist to zero
   *
   * Note that no control of memory allocation is made
   */
  static void bit0(int* list, int num);

  /**
   * @brief bit-wise logical AND between two arrays of bits
   *
   * (1 AND 1 = 1; 1 AND 0 = 0 AND 1 = 0 AND 0 = 0)
   *
   * Note that no control of memory allocation is made
   *
   * param list1 first array of bit
   * param list2 second array of bit
   * param listet resulting array of bit
   * param len number of int over which the operation is performed
   */
  static void bitAnd(int* listet, int* list1, int* list2, size_t len);

  /**
   * @brief bit-wise logical OR between two arrays of bits
   *
   * (1 OR 1 = 1 OR 0 = 0 OR 1 = 1; 0 OR 0 = 0)
   *
   * Note that no control of memory allocation is made
   *
   * param list1 first array of bit
   * param list2 second array of bit
   * param listou resulting array of bit
   * param len number of int over which the operation is performed
   */
  static void bitOr(int* listou, int* list1, int* list2, size_t len);

  /**
   * @brief bit-wise logical NOT
   *
   * (NOT 1 = 0; NOT 0 = 1)
   *
   * Note that no control of memory allocation is made
   *
   * param list input array of bit
   * param listnot resulting array of bit
   * param len number of int over which the operation is performed
   */
  static void bitNot(int* listnon, int* list, size_t len);

  /**
   * @brief Tells whether bit number num in bit array list is one
   */
  static bool testBit(int* list, int num);

  /**
   * @brief Makes one BipartitionList out of several
   *
   * The input BipartitionList objects must share the same set of elements. This will be checked or not
   * depending on checkElements
   */
  static std::unique_ptr<BipartitionList> mergeBipartitionLists(const std::vector<std::unique_ptr<BipartitionList>>& vecBipartL, bool checkElements = true);

  /**
   * @brief Construct a BipartitionList containing two bipartitions taken from distinct input lists
   */
  static std::unique_ptr<BipartitionList> buildBipartitionPair(const BipartitionList& bipartL1, size_t i1, const BipartitionList& bipartL2, size_t i2, bool checkElements = true);

  /**
   * @brief Tells whether two bipartitions from distinct lists are identical
   */
  static bool areIdentical(const BipartitionList& bipart1, size_t i1, const BipartitionList& bipart2, size_t i2, bool checkElements = true);

  /**
   * @brief Tells whether two bipartitions from distinct lists are compatible
   */
  static bool areCompatible(const BipartitionList& bipart1, size_t i1, const BipartitionList& bipart2, size_t i2, bool checkElements = true);

  /**
   * @brief Create a sequence data set corresponding to the Matrix Representation of the input BipartitionList objects
   *
   * The input BipartitionList objects can have distinct sets of elements - missing data will be represented as 'N'.
   * The output alignment (DNA sequences including only A, C and N)) is ready for maximum parsimony analysis
   * according to the MRP supertree method.
   */
  static std::unique_ptr<VectorSiteContainer> MRPEncode(
      const std::vector<std::unique_ptr<BipartitionList>>& vecBipartL);

  /**
   * @brief Create a sequence data set corresponding to the Matrix Representation of the input BipartitionList objects and accomodates multilabel trees
   *
   * The input BipartitionList objects can have distinct sets of elements - missing data will be represented as 'N'.
   * The output alignment (DNA sequences including only A, C and N)) is ready for maximum parsimony analysis
   * according to the MRP supertree method.
   */
  static std::unique_ptr<VectorSiteContainer> MRPEncodeMultilabel(
      const std::vector<std::unique_ptr<BipartitionList>>& vecBipartL);
};
} // end of namespace bpp.
#endif // BPP_PHYL_TREE_BIPARTITIONTOOLS_H
