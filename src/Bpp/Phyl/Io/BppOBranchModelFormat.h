// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_BPPOBRANCHMODELFORMAT_H
#define BPP_PHYL_IO_BPPOBRANCHMODELFORMAT_H


#include "BppOTransitionModelFormat.h"

namespace bpp
{
/**
 * @brief Branch model I/O in BppO format.
 *
 * Creates a new branch model object according to model description syntax
 * (see the Bio++ Progam Suite manual for a detailed description of this syntax).
 *
 */

class BppOBranchModelFormat :
  public BppOTransitionModelFormat
{
public:
  /**
   * @brief Create a new BppOBranchModelFormat object.
   *
   * @param alphabetCode     Bit saying which alphabets are allowed in the model specification.
   * @param allowCovarions   Tell is a covarion model can be returned.
   * @param allowMixed       Tell is a mixture model can be returned.
   * @param allowGaps        Tell is a gap model can be returned.
   * @param verbose          Tell if the construction is verbose.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   */
  BppOBranchModelFormat(unsigned char alphabetCode, bool allowCovarions, bool allowMixed, bool allowGaps, bool verbose, int warn) :
    BppOTransitionModelFormat(alphabetCode, allowCovarions, allowMixed, allowGaps, verbose, warn)
  {}

  BppOBranchModelFormat(const BppOTransitionModelFormat& format) :
    BppOTransitionModelFormat(format)
  {}

  BppOBranchModelFormat& operator=(const BppOBranchModelFormat& format)
  {
    BppOTransitionModelFormat::operator=(format);
    return *this;
  }

public:
  std::unique_ptr<BranchModelInterface> readBranchModel(
    std::shared_ptr<const Alphabet> alphabet,
    const std::string& modelDescription,
    const AlignmentDataInterface& data,
    bool parseArguments = true);
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_BPPOBRANCHMODELFORMAT_H
