// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_BPPOMULTITREEWRITERFORMAT_H
#define BPP_PHYL_IO_BPPOMULTITREEWRITERFORMAT_H


#include "IoTreeFactory.h"

namespace bpp
{
/**
 * @brief Tree I/O in BppO format.
 *
 * Creates a new IMultiTree object according to
 * distribution description syntax (see the Bio++ Program Suite
 * manual for a detailed description of this syntax).
 *
 */
class BppOMultiTreeWriterFormat :
  public virtual IOFormat
{
private:
  std::map<std::string, std::string> unparsedArguments_;
  int warningLevel_;

public:
  BppOMultiTreeWriterFormat(int warningLevel) :
    unparsedArguments_(), warningLevel_(warningLevel) {}

  virtual ~BppOMultiTreeWriterFormat() {}

public:
  const std::string getFormatName() const { return "BppO"; }

  const std::string getFormatDescription() const { return "Bpp Options format."; }

  const std::string getDataType() const { return "Trees writer"; }

  /**
   * @brief Read a OMultiTree object from a string.
   *
   * @param description A string describing the reader in the keyval syntax.
   * @return A new OMultiTree object according to options specified.
   * @throw Exception if an error occured.
   */
  OMultiTree* readOMultiTree(const std::string& description);

  /**
   * @return The arguments and their unparsed values from the last call of the read function, if there are any.
   */
  virtual const std::map<std::string, std::string>& getUnparsedArguments() const { return unparsedArguments_; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_BPPOMULTITREEWRITERFORMAT_H
