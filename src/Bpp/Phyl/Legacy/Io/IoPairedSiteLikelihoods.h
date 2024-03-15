// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_IO_IOPAIREDSITELIKELIHOODS_H
#define BPP_PHYL_LEGACY_IO_IOPAIREDSITELIKELIHOODS_H



// From the STL:
#include <vector>
#include <string>
#include <sstream>

// From Bio++
#include <Bpp/Io/IoFormat.h>
#include "../Likelihood/PairedSiteLikelihoods.h"

namespace bpp
{
/**
 * @brief Base class for I/O on paired-site likelihoods.
 *
 * @author Nicolas Rochette
 */
class IOPairedSiteLikelihoods : public virtual IOFormat
{};


/**
 * @brief This class provides I/O for the Tree-Puzzle/RAxML (phylip-like) paired-site likelihoods format.
 *
 * @author Nicolas Rochette
 */
class IOTreepuzzlePairedSiteLikelihoods : public virtual IOPairedSiteLikelihoods
{
public:
  /**
   * @brief Read paired-site likelihoods from a Treepuzzle/RAxML-formatted stream.
   *
   * @throw Exception If the format is not recognized.
   */
  static PairedSiteLikelihoods readPairedSiteLikelihoods(std::istream& is);

  /**
   * @brief Read paired-site likelihoods from a Treepuzzle/RAxML-formatted file.
   *
   * @throw Exception If the format is not recognized.
   */
  static PairedSiteLikelihoods readPairedSiteLikelihoods(const std::string& path);

  /**
   * @brief Write paired-site likelihoods to a stream.
   *
   * @param psl The PairedSiteLikelihoods object to write.
   * @param os The output stream.
   * @param delim The delimiter between model names and likelihoods. The defaut is a tab but two spaces might be used.
   */
  static void writePairedSiteLikelihoods(const PairedSiteLikelihoods& psl, std::ostream& os, const std::string& delim = "\t");

  /**
   * @brief Write paired-site likelihoods to a file.
   *
   * @param psl The PairedSiteLikelihoods object to write.
   * @param path The path of the output file.
   * @param delim The delimiter between model names and likelihoods. (The defaut is a tab but two spaces might be used.)
   */
  static void writePairedSiteLikelihoods(const PairedSiteLikelihoods& psl, const std::string& path, const std::string& delim = "\t");
};


/**
 * @brief This class provides input for the Phyml paired-site likelihoods format.
 *
 * @author Nicolas Rochette
 */
class IOPhymlPairedSiteLikelihoods : public virtual IOPairedSiteLikelihoods
{
public:
  /**
   * @brief Read paired-site likelihoods from a Phyml-formatted stream.
   * @throw Exception If the format is not recognized.
   */
  static std::vector<double> readPairedSiteLikelihoods(std::istream& is);

  /**
   * @brief Read Phyml paired-site likelihoods from a Phyml-formatted file.
   * @throw Exception If the format is not recognized.
   */
  static std::vector<double> readPairedSiteLikelihoods(const std::string& path);
};
} // namespace bpp
#endif // BPP_PHYL_LEGACY_IO_IOPAIREDSITELIKELIHOODS_H
