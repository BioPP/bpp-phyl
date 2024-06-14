// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PARTITIONSEQUENCEEVOLUTION_H
#define BPP_PHYL_LIKELIHOOD_PARTITIONSEQUENCEEVOLUTION_H


#include "MultiProcessSequenceEvolution.h"

namespace bpp
{
/**
 * @brief Sequence evolution framework based on a mixture of
 * substitution processes
 *
 * @see MultiProcessSequenceEvolution
 */

class PartitionSequenceEvolution :
  public MultiProcessSequenceEvolution
{
private:
  /**
   * @brief vector of the substitution process numbers along the sequence.
   */
  std::vector<size_t> vProc_;

  size_t vSize_;

  /**
   * @brief On the reverse, for each process number, the vector
   * of the sites where it is used.
   *
   * Convenient for process specific site patterns.
   */

  std::map<size_t, std::vector<size_t>> mProcPos_;

public:
  /*
   * @brief constructor
   *
   * @param the used SubstitutionProcessCollection
   * @param A vector of the number of the processes along the sequence.
   */
  PartitionSequenceEvolution(
      std::shared_ptr<SubstitutionProcessCollection> processColl,
      std::vector<size_t>& posProc);

  PartitionSequenceEvolution(const PartitionSequenceEvolution& mlc) :
    MultiProcessSequenceEvolution(mlc),
    vProc_(mlc.vProc_),
    vSize_(mlc.vSize_),
    mProcPos_(mlc.mProcPos_) {}

  PartitionSequenceEvolution& operator=(const PartitionSequenceEvolution& mlc)
  {
    MultiProcessSequenceEvolution::operator=(mlc);
    vProc_ = mlc.vProc_;
    vSize_ = mlc.vSize_;

    mProcPos_ = mlc.mProcPos_;

    return *this;
  }

  virtual ~PartitionSequenceEvolution() {}

  PartitionSequenceEvolution* clone() const { return new PartitionSequenceEvolution(*this); }

public:
  size_t getNumberOfSites() const
  {
    return vSize_;
  }

  const std::vector<size_t>& getProcessNumbersPerSite() const
  {
    return vProc_;
  }

  /*
   * @brief Get Substitution Process Number for a given site.
   * @param i the index of the site
   *
   */
  size_t getSubstitutionProcessNumber(size_t i) const
  {
    if (i >= vSize_)
      throw IndexOutOfBoundsException("PartitionSequenceEvolution::getSubstitutionProcess", i, 0, vSize_);
    return vProc_[i];
  }

  std::map<size_t, std::vector<size_t>>& mapOfProcessSites()
  {
    return mProcPos_;
  }

  const std::map<size_t, std::vector<size_t>>& mapOfProcessSites() const
  {
    return mProcPos_;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PARTITIONSEQUENCEEVOLUTION_H
