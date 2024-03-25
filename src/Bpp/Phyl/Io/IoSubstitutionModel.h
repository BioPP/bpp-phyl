// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_IOSUBSTITUTIONMODEL_H
#define BPP_PHYL_IO_IOSUBSTITUTIONMODEL_H


#include "../Model/SubstitutionModel.h"

// From bpp-core:
#include <Bpp/Exceptions.h>
#include <Bpp/Io/IoFormat.h>
#include <Bpp/Io/OutputStream.h>

// From bpp-seq:
#include <Bpp/Seq/Container/AlignmentData.h>

namespace bpp
{
class TransitionModelInterface;
class SubstitutionModelInterface;

/**
 * @brief General interface for model I/O.
 */
class IoSubstitutionModel :
  public virtual IOFormat
{
public:
  IoSubstitutionModel() {}
  virtual ~IoSubstitutionModel() {}

public:
  virtual const std::string getDataType() const { return "Substitution Model"; }
};

/**
 * @brief General interface for model readers.
 */

class ISubstitutionModel :
  public virtual IoSubstitutionModel
{
public:
  ISubstitutionModel() {}
  virtual ~ISubstitutionModel() {}

public:
  /**
   * @brief Read a substitution model from a string.
   *
   * @param alphabet         The alphabet to use in the model.
   * @param modelDescription A string describing the model in the format.
   * @param data             A pointer toward a AlignedValuesContainer, which can be used to initial some parameters like frequencies.
   * @param parseArguments Attempt to parse function arguments. If not, only store them and use default values instead.
   * @return A new SubstitutionModel object according to options specified.
   * @throw Exception if an error occurred.
   */
  virtual std::unique_ptr<SubstitutionModelInterface> readSubstitutionModel(
      std::shared_ptr<const Alphabet> alphabet,
      const std::string& modelDescription,
      const AlignmentDataInterface& data,
      bool parseArguments = true) = 0;

  /**
   * @return The arguments and their unparsed values from the last call of the read function, if there are any.
   */
  virtual const std::map<std::string, std::string>& getUnparsedArguments() const = 0;
};

/**
 * @brief General interface for distance matrix writers.
 */
class OSubstitutionModel :
  public virtual IoSubstitutionModel
{
public:
  OSubstitutionModel() {}
  virtual ~OSubstitutionModel() {}

public:
  /**
   * @brief Write a substitution model to a stream.
   *
   * @param model A substitution model object;
   * @param out The output stream;
   * @param globalAliases parameters linked to global alias.
   * @param writtenNames is the vector of the written
   * parameters so far [in, out];
   * @throw Exception if an error occurred.
   */
  virtual void write(
      const BranchModelInterface& model,
      OutputStream& out,
      std::map<std::string, std::string>& globalAliases,
      std::vector<std::string>& writtenNames) const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_IOSUBSTITUTIONMODEL_H
