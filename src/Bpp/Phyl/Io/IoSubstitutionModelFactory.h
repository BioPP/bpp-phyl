// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_IOSUBSTITUTIONMODELFACTORY_H
#define BPP_PHYL_IO_IOSUBSTITUTIONMODELFACTORY_H

#include <Bpp/Exceptions.h>

#include "../Model/SubstitutionModel.h"
#include "IoSubstitutionModel.h"

// From the STL:
#include <string>

namespace bpp
{
/**
 * @brief Utilitary class for creating substitution model readers and
 * writers.
 *
 * @see IOSequenceFactory
 * @see IOTreeFactory
 */
class IOSubstitutionModelFactory
{
public:
  static const std::string BPPO_FORMAT;

public:
  /**
   * @brief Creates a new factory object.
   *
   * Example:
   * @code
   * ISubstitutionModel * modReader = IOSubstitutionModelFactory().createReader(IOSubstitutionModelFactory::BPP_FORMAT);
   * SubstitutionModel * model = modReader->read(...);
   * delete modReader;
   * @endcode
   */
  IOSubstitutionModelFactory() {}
  virtual ~IOSubstitutionModelFactory() {}

  /**
   * @brief Get a new dynamically created ISubstitutionModel object.
   *
   * @param format The input file format.
   * @return A pointer toward a new ISubstitutionModel object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual ISubstitutionModel* createReader(const std::string& format);

  /**
   * @brief Get a new dynamically created OSubstitutionModel object.
   *
   * @param format The output file format.
   * @return A pointer toward a new OSubstitutionModel object.
   * @throw Exception If the format name do not match any available format.
   */
  virtual OSubstitutionModel* createWriter(const std::string& format);
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_IOSUBSTITUTIONMODELFACTORY_H
