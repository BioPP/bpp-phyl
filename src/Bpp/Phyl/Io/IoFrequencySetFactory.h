// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_IOFREQUENCYSETFACTORY_H
#define BPP_PHYL_IO_IOFREQUENCYSETFACTORY_H

#include <Bpp/Exceptions.h>

#include "../Model/SubstitutionModel.h"
#include "IoFrequencySet.h"

// From the STL:
#include <string>

namespace bpp
{
/**
 * @brief Utilitary class for creating frequencies set readers and writers.
 *
 * @see IOSequenceFactory
 * @see IOTreeFactory
 */
class IOFrequencySetFactory
{
public:
  static const std::string BPPO_FORMAT;

public:
  /**
   * @brief Creates a new factory object.
   *
   * Example:
   * @code
   * IFrequencySet * freqReader = IOFrequencySetFactory().createReader(IOFrequencySetFactory::BPP_FORMAT);
   * FrequencySet * freqset = freqReader->read(...);
   * delete freqReader;
   * @endcode
   */
  IOFrequencySetFactory() {}
  virtual ~IOFrequencySetFactory() {}

  /**
   * @brief Get a new dynamically created IFrequencySet object.
   *
   * @param format The input file format.
   * @return A pointer toward a new IFrequencySet object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual IFrequencySet* createReader(const std::string& format);

  /**
   * @brief Get a new dynamically created OFrequencySet object.
   *
   * @param format The output file format.
   * @return A pointer toward a new OFrequencySet object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual OFrequencySet* createWriter(const std::string& format);
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_IOFREQUENCYSETFACTORY_H
