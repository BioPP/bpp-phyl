// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_IO_IODISTANCEMATRIXFACTORY_H
#define BPP_PHYL_LEGACY_IO_IODISTANCEMATRIXFACTORY_H

#include <Bpp/Exceptions.h>

#include "IoDistanceMatrix.h"

// From the STL:
#include <string>

namespace bpp
{
/**
 * @brief Utilitary class for creating distance matrix readers and writers.
 *
 * @see IOSequenceFactory
 * @see IOTreeFactory
 */
class IODistanceMatrixFactory
{
public:
  static const std::string PHYLIP_FORMAT;

public:
  /**
   * @brief Creates a new factory object.
   *
   * Example:
   * @code
   * IDistanceMatrix * matReader = IODistanceMatrixFactory().createReader(IODistanceMatrixFactory::PHYLIP);
   * DistanceMatrix * matrix = matReader->read("file.ph");
   * delete matReader;
   * @endcode
   */
  IODistanceMatrixFactory() {}
  virtual ~IODistanceMatrixFactory() {}

  /**
   * @brief Get a new dynamically created IDistanceMatrix object.
   *
   * @param format The input file format, and whether names should be
   *      only less than 10 characters, or not (false=10 characters max).
   * @param extended format (default false).
   * @return A pointer toward a new IDistanceMatrix object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual IDistanceMatrix* createReader(const std::string& format, bool extended = false);

  /**
   * @brief Get a new dynamically created ODistanceMatrix object.
   *
   * @param format The output file format, and whether names should be
   *        only less than 10 characters, or not (false=10 characters max).
   * @param extended format (default false).
   * @return A pointer toward a new ODistanceMatrix object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual ODistanceMatrix* createWriter(const std::string& format, bool extended = false);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_IO_IODISTANCEMATRIXFACTORY_H
