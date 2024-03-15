// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_IOTREEFACTORY_H
#define BPP_PHYL_IO_IOTREEFACTORY_H


#include "IoTree.h"

namespace bpp
{
/**
 * @brief Utilitary class for creating tree readers and writers.
 *
 * @see IOSequenceFactory
 * @see IODistanceMatrixFactory
 */
class IOTreeFactory
{
public:
  static const std::string NEWICK_FORMAT;
  static const std::string NEXUS_FORMAT;
  static const std::string NHX_FORMAT;

public:
  /**
   * @brief Creates a new factory object.
   *
   * Example:
   * @code
   * ITree * treeReader = IOTreeFactory().createReader(IOTreeFactory::NEWICK);
   * Tree * tree = treeReader->read("file.dnd");
   * delete treeReader;
   * @endcode
   */
  IOTreeFactory() {}
  virtual ~IOTreeFactory() {}

  /**
   * @brief Get a new dynamically created ITree object.
   *
   * @param format The input file format.
   * @return A pointer toward a new ITree object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual ITree* createReader(const std::string& format);

  /**
   * @brief Get a new dynamically created OTree object.
   *
   * @param format The output file format.
   * @return A pointer toward a new OTree object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual OTree* createWriter(const std::string& format);
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_IOTREEFACTORY_H
