// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_BPPOFREQUENCYSETFORMAT_H
#define BPP_PHYL_IO_BPPOFREQUENCYSETFORMAT_H


#include "IoFrequencySetFactory.h"

// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

namespace bpp
{
/**
 * @brief Frequencies set I/O in BppO format.
 *
 * Allow to create a new frequencies set object according to model description syntax
 * (see the Bio++ Progam Suite manual for a detailed description of this syntax).
 */
class BppOFrequencySetFormat :
  public virtual IFrequencySet,
  public virtual OFrequencySet
{
public:
  static unsigned char DNA;
  static unsigned char RNA;
  static unsigned char NUCLEOTIDE;
  static unsigned char PROTEIN;
  static unsigned char CODON;
  static unsigned char WORD;
  static unsigned char ALL;

private:
  unsigned char alphabetCode_;
  bool verbose_;
  std::map<std::string, std::string> unparsedArguments_;
  std::shared_ptr<const GeneticCode> geneticCode_;
  int warningLevel_;

public:
  BppOFrequencySetFormat(unsigned char alphabetCode, bool verbose, int warn) :
    alphabetCode_(alphabetCode),
    verbose_(verbose),
    unparsedArguments_(),
    geneticCode_(0),
    warningLevel_(warn)
  {}

  BppOFrequencySetFormat(const BppOFrequencySetFormat& format) :
    alphabetCode_(format.alphabetCode_),
    verbose_(format.verbose_),
    unparsedArguments_(format.unparsedArguments_),
    geneticCode_(format.geneticCode_),
    warningLevel_(format.warningLevel_)
  {}

  BppOFrequencySetFormat& operator=(const BppOFrequencySetFormat& format)
  {
    alphabetCode_      = format.alphabetCode_;
    verbose_           = format.verbose_;
    unparsedArguments_ = format.unparsedArguments_;
    geneticCode_       = format.geneticCode_;
    warningLevel_      = format.warningLevel_;
    return *this;
  }

  virtual ~BppOFrequencySetFormat() {}

public:
  const std::string getFormatName() const override { return "BppO"; }

  const std::string getFormatDescription() const override { return "Bpp Options format."; }

  /**
   * @brief Set the genetic code to use in case a codon frequencies set should be built.
   *
   * @param gCode The genetic code to use.
   */
  void setGeneticCode(std::shared_ptr<const GeneticCode> gCode)
  {
    geneticCode_ = gCode;
  }

  std::unique_ptr<FrequencySetInterface> readFrequencySet(
      std::shared_ptr<const Alphabet> alphabet,
      const std::string& freqDescription,
      const AlignmentDataInterface& data,
      bool parseArguments = true) override;

  const std::map<std::string, std::string>& getUnparsedArguments() const override
  {
    return unparsedArguments_;
  }

  void writeFrequencySet(
      const FrequencySetInterface& freqset,
      OutputStream& out,
      std::map<std::string, std::string>& globalAliases,
      std::vector<std::string>& writtenNames) const override;

  void setVerbose(bool verbose) { verbose_ = verbose; }

private:
  void initialize_(FrequencySetInterface& freqSet, const AlignmentDataInterface& data);
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_BPPOFREQUENCYSETFORMAT_H
