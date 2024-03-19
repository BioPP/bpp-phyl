// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_BPPOSUBSTITUTIONMODELFORMAT_H
#define BPP_PHYL_IO_BPPOSUBSTITUTIONMODELFORMAT_H


#include "IoSubstitutionModelFactory.h"
#include "../Model/MixedTransitionModel.h"

// From bpp-seq
#include <Bpp/Seq/Container/AlignmentData.h>
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

namespace bpp
{
/**
 * @brief Substitution model I/O in BppO format.
 *
 * Creates a new substitution model object according to model description syntax
 * (see the Bio++ Progam Suite manual for a detailed description of this syntax).
 */
class BppOSubstitutionModelFormat :
  public ISubstitutionModel,
  public OSubstitutionModel
{
public:
  static unsigned char DNA;
  static unsigned char RNA;
  static unsigned char NUCLEOTIDE;
  static unsigned char PROTEIN;
  static unsigned char CODON;
  static unsigned char WORD;
  static unsigned char BINARY;
  static unsigned char INTEGER;
  static unsigned char ALL;

protected:
  unsigned char alphabetCode_;
  bool allowCovarions_;
  bool allowMixed_;
  bool allowGaps_;
  bool verbose_;
  std::map<std::string, std::string> unparsedArguments_;
  std::shared_ptr<const GeneticCode> geneticCode_;
  int warningLevel_;

public:
  /**
   * @brief Create a new BppOSubstitutionModelFormat object.
   *
   * @param alphabetCode     Bit saying which alphabets are allowed in the model specification.
   * @param allowCovarions   Tell is a covarion model can be returned.
   * @param allowMixed       Tell is a mixture model can be returned.
   * @param allowGaps        Tell is a gap model can be returned.
   * @param verbose          Tell if the construction is verbose.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   */
  BppOSubstitutionModelFormat(unsigned char alphabetCode, bool allowCovarions, bool allowMixed, bool allowGaps, bool verbose, int warn) :
    alphabetCode_(alphabetCode),
    allowCovarions_(allowCovarions),
    allowMixed_(allowMixed),
    allowGaps_(allowGaps),
    verbose_(verbose),
    unparsedArguments_(),
    geneticCode_(0),
    warningLevel_(warn)
  {}

  BppOSubstitutionModelFormat(const BppOSubstitutionModelFormat& format) :
    alphabetCode_(format.alphabetCode_),
    allowCovarions_(format.allowCovarions_),
    allowMixed_(format.allowMixed_),
    allowGaps_(format.allowGaps_),
    verbose_(format.verbose_),
    unparsedArguments_(format.unparsedArguments_),
    geneticCode_(format.geneticCode_),
    warningLevel_(format.warningLevel_)
  {}

  BppOSubstitutionModelFormat& operator=(const BppOSubstitutionModelFormat& format)
  {
    alphabetCode_      = format.alphabetCode_;
    allowCovarions_    = format.allowCovarions_;
    allowMixed_        = format.allowMixed_;
    allowGaps_         = format.allowGaps_;
    verbose_           = format.verbose_;
    unparsedArguments_ = format.unparsedArguments_;
    geneticCode_       = format.geneticCode_;
    warningLevel_      = format.warningLevel_;
    return *this;
  }

  virtual ~BppOSubstitutionModelFormat() {}

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

  std::unique_ptr<SubstitutionModelInterface> readSubstitutionModel(
    std::shared_ptr<const Alphabet> alphabet,
    const std::string& modelDescription,
    const AlignmentDataInterface& data,
    bool parseArguments = true) override;

  const std::map<std::string, std::string>& getUnparsedArguments() const override
  {
    return unparsedArguments_;
  }

  /**
   * @brief Write a substitution model to a stream.
   *
   * @param model A substitution model object;
   * @param out The output stream;
   * @param globalAliases parameters linked to global alias. The
   * output will be "name=alias_name";
   * @param writtenNames is the vector of the written
   * parameters so far [in, out];
   * @throw Exception If an error occured.
   */
  void write(const BranchModelInterface& model,
      OutputStream& out,
      std::map<std::string, std::string>& globalAliases,
      std::vector<std::string>& writtenNames) const override;

  void setVerbose(bool verbose) { verbose_ = verbose;}

private:
  std::unique_ptr<SubstitutionModelInterface> readWord_(
    std::shared_ptr<const Alphabet> alphabet,
    const std::string& modelDescription,
    const AlignmentDataInterface& data);

  void writeMixed_(const MixedTransitionModelInterface& model,
      OutputStream& out,
      std::map<std::string, std::string>& globalAliases,
      std::vector<std::string>& writtenNames) const;

protected:
  /**
   * @brief Finish parsing of parameters, taking care of aliases.
   */
  void updateParameters_(
    BranchModelInterface& model,
    std::map<std::string, std::string>& args);

  /**
   * @brief Set parameter initial values of a given model according to options.
   *
   * Parameters actually depends on the model passed as argument.
   * See getSubstitutionModel for more information.
   *
   * This function is mainly for internal usage, you're probably looking for the getSubstitutionModel or getSubstitutionModelSet function.
   *
   * @param model                   The model to set.
   * @param data   A pointer toward the AlignmentDataInterface for which the substitution model is designed.
   *               The alphabet associated to the data must be of the same type as the one specified for the model.
   *               May be equal to NULL, but in this case use_observed_freq option will be unavailable.
   * @throw Exception if an error occured.
   */
  void initialize_(
    BranchModelInterface& model,
    const AlignmentDataInterface& data);
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_BPPOSUBSTITUTIONMODELFORMAT_H
