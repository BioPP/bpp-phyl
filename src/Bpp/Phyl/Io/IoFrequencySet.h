// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_IOFREQUENCYSET_H
#define BPP_PHYL_IO_IOFREQUENCYSET_H


#include "../Model/FrequencySet/FrequencySet.h"

// From bpp-core:
#include <Bpp/Exceptions.h>
#include <Bpp/Io/IoFormat.h>
#include <Bpp/Io/OutputStream.h>

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainer.h>

namespace bpp
{
/**
 * @brief General interface for model I/O.
 */
class IoFrequencySet :
  public virtual IOFormat
{
public:
  IoFrequencySet() {}
  virtual ~IoFrequencySet() {}

public:
  virtual const std::string getDataType() const { return "Frequencies Set"; }
};

/**
 * @brief General interface for distance matrix readers.
 */
class IFrequencySet :
  public virtual IoFrequencySet
{
public:
  IFrequencySet() {}
  virtual ~IFrequencySet() {}

public:
  /**
   * @brief Read a frequencies set from a string.
   *
   * @param alphabet         The alphabet to use in the model.
   * @param freqDescription  A string describing the frequencies set.
   * @param data             A SiteContainer with the data to use to initialize frequency parameters. Can be set to 0.
   * @param parseArguments   Attempt to parse function arguments. If not, only store them and use default values instead.
   * @return A new FrequencySet object according to options specified.
   * @throw Exception if an error occurred.
   */
  virtual std::unique_ptr<FrequencySetInterface> readFrequencySet(
      std::shared_ptr<const Alphabet> alphabet,
      const std::string& freqDescription,
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
class OFrequencySet :
  public virtual IoFrequencySet
{
public:
  OFrequencySet() {}
  virtual ~OFrequencySet() {}

public:
  /**
   * @brief Write a substitution model to a stream.
   *
   * @param pfreqset A frequency set object;
   * @param out The output stream;
   * @param globalAliases parameters linked to global alias. The
   * output will be "name=alias_name";
   * @param writtenNames is the vector of the written
   * parameters so far [in, out];
   */
  virtual void writeFrequencySet(
      const FrequencySetInterface& pfreqset,
      OutputStream& out,
      std::map<std::string, std::string>& globalAliases,
      std::vector<std::string>& writtenNames) const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_IOFREQUENCYSET_H
