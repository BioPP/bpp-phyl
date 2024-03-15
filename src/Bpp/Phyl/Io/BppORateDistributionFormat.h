// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_BPPORATEDISTRIBUTIONFORMAT_H
#define BPP_PHYL_IO_BPPORATEDISTRIBUTIONFORMAT_H

#include <Bpp/Io/BppODiscreteDistributionFormat.h>


namespace bpp
{
/**
 * @brief Rate Distribution I/O in BppO format.
 *
 * Creates a new discrete distribution object according to
 * distribution description syntax (see the Bio++ Progam Suite
 * manual for a detailed description of this syntax).
 *
 * Rate distributions are normalized and have a mean of 1, so that branch lengths are measured in mean number of substitutions per site.
 *
 * @see BppODiscreteDistribtution for a more generic parser.
 *
 */
class BppORateDistributionFormat :
  public BppODiscreteDistributionFormat
{
private:
  bool allowConstant_;

public:
  /**
   * @brief Build a new BppORateDistributionFormat object.
   *
   * @param allowConstant Is contant distribution allowed.
   */
  BppORateDistributionFormat(bool allowConstant) :
    BppODiscreteDistributionFormat(),
    allowConstant_(allowConstant)
  {}

  virtual ~BppORateDistributionFormat() {}

public:
  std::unique_ptr<DiscreteDistributionInterface> readDiscreteDistribution(const std::string& distDescription, bool parseArguments);

  void writeDiscreteDistribution(const DiscreteDistributionInterface& dist,
                                 OutputStream& out,
                                 std::map<std::string, std::string>& globalAliases,
                                 std::vector<std::string>& writtenNames) const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_BPPORATEDISTRIBUTIONFORMAT_H
