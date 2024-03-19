// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_BPPOTRANSITIONMODELFORMAT_H
#define BPP_PHYL_IO_BPPOTRANSITIONMODELFORMAT_H


#include "../Model/MixedTransitionModel.h"
#include "BppOSubstitutionModelFormat.h"

namespace bpp
{
/**
 * @brief Transition model I/O in BppO format.
 *
 * Creates a new transition model object according to model description syntax
 * (see the Bio++ Progam Suite manual for a detailed description of this syntax).
 */
class BppOTransitionModelFormat :
  public BppOSubstitutionModelFormat
{
private:
  std::unique_ptr<MixedTransitionModelInterface> readMixed_(
    std::shared_ptr<const Alphabet> alphabet,
    const std::string& modelDescription,
    const AlignmentDataInterface& data);

public:
  /**
   * @brief Create a new BppOTransitionModelFormat object.
   *
   * @param alphabetCode     Bit saying which alphabets are allowed in the model specification.
   * @param allowCovarions   Tell is a covarion model can be returned.
   * @param allowMixed       Tell is a mixture model can be returned.
   * @param allowGaps        Tell is a gap model can be returned.
   * @param verbose          Tell if the construction is verbose.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   */
  BppOTransitionModelFormat(unsigned char alphabetCode, bool allowCovarions, bool allowMixed, bool allowGaps, bool verbose, int warn) :
    BppOSubstitutionModelFormat(alphabetCode, allowCovarions, allowMixed, allowGaps, verbose, warn)
  {}

  BppOTransitionModelFormat(const BppOTransitionModelFormat& format) :
    BppOSubstitutionModelFormat(format)
  {}

  BppOTransitionModelFormat& operator=(const BppOTransitionModelFormat& format)
  {
    BppOSubstitutionModelFormat::operator=(format);
    return *this;
  }

  virtual ~BppOTransitionModelFormat() {}

public:
  std::unique_ptr<TransitionModelInterface> readTransitionModel(
    std::shared_ptr<const Alphabet> alphabet,
    const std::string& modelDescription,
    const AlignmentDataInterface& data,
    bool parseArguments = true);
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_BPPOTRANSITIONMODELFORMAT_H
