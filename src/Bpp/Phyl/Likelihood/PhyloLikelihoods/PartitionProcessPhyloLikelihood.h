// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PARTITIONPROCESSPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PARTITIONPROCESSPHYLOLIKELIHOOD_H


#include "../PartitionSequenceEvolution.h"
#include "AlignedPhyloLikelihoodProduct.h"
#include "SequencePhyloLikelihood.h"
#include "SingleProcessPhyloLikelihood.h"

// From bpp-seq:
#include <Bpp/Seq/Container/AlignmentData.h>

namespace bpp
{
/**
 * @brief Likelihood framework based on a partition of a sequence in
 * simple likelihoods.
 *
 * @see AbstractSequencePhyloLikelihood
 */

struct ProcPos
{
  size_t nProc;
  size_t pos;
};


class PartitionProcessPhyloLikelihood :
  public AbstractSequencePhyloLikelihood,
  public AbstractPhyloLikelihoodSet
{
private:
  /**
   * @brief to avoid the dynamic casts
   */
  std::shared_ptr<PartitionSequenceEvolution> mSeqEvol_;

  /**
   * vector of couples <number of process, site> specific to
   * this partition process.
   */
  std::vector<ProcPos> vProcPos_;

  /**
   * map nbe of phylolikelihood : created AlignmentDataInterface<std::string>
   */
  std::map<size_t, std::shared_ptr<AlignmentDataInterface>> mData_;

  /**
   * AlignedLikelihoodCalculation to store DF nodes
   */
  mutable std::shared_ptr<AlignedLikelihoodCalculation> likCal_;

public:
  PartitionProcessPhyloLikelihood(
      Context& context,
      std::shared_ptr<PartitionSequenceEvolution> processSeqEvol,
      size_t nSeqEvol = 0);

  PartitionProcessPhyloLikelihood(
      Context& context,
      std::shared_ptr<const AlignmentDataInterface> data,
      std::shared_ptr<PartitionSequenceEvolution> processSeqEvol,
      size_t nSeqEvol = 0,
      size_t nData = 0);

  PartitionProcessPhyloLikelihood(
      std::shared_ptr<const AlignmentDataInterface> data,
      std::shared_ptr<PartitionSequenceEvolution> processSeqEvol,
      std::shared_ptr<CollectionNodes> collNodes,
      size_t nSeqEvol = 0,
      size_t nData = 0);

protected:
  PartitionProcessPhyloLikelihood(const PartitionProcessPhyloLikelihood& lik) :
    AbstractPhyloLikelihood(lik),
    AbstractAlignedPhyloLikelihood(lik),
    AbstractSingleDataPhyloLikelihood(lik),
    AbstractParametrizable(""),
    AbstractSequencePhyloLikelihood(lik),
    AbstractPhyloLikelihoodSet(lik),
    mSeqEvol_(lik.mSeqEvol_),
    vProcPos_(lik.vProcPos_),
    mData_(lik.mData_),
    likCal_(lik.likCal_)
  {}

  PartitionProcessPhyloLikelihood& operator=(const PartitionProcessPhyloLikelihood& lik)
  {
    AbstractSequencePhyloLikelihood::operator=(lik);
    AbstractPhyloLikelihoodSet::operator=(lik);
    mSeqEvol_ = lik.mSeqEvol_;
    vProcPos_ = lik.vProcPos_;
    mData_    = lik.mData_;
    likCal_   = lik.likCal_;
    return *this;
  }

  PartitionProcessPhyloLikelihood* clone() const override { return new PartitionProcessPhyloLikelihood(*this); }

public:
  virtual ~PartitionProcessPhyloLikelihood() {}

  /**
   * @brief Set the dataset for which the likelihood must be evaluated.
   *
   * @param data The data set to use.
   * @param nData the number of the data (optional, default = 0)
   */

  void setData(std::shared_ptr<const AlignmentDataInterface> data, size_t nData = 0) override;

  /*
   * @brief Get PhyloLikelihood Number for a given site.
   * @param siteIndex the index of the site
   *
   */
  std::shared_ptr<const SingleProcessPhyloLikelihood> getPhyloLikelihoodForASite(size_t siteIndex) const
  {
    return dynamic_pointer_cast<const SingleProcessPhyloLikelihood>(getPhyloLikelihood(vProcPos_[siteIndex].nProc));
  }

  /*
   * @brief Get LikelihoodCaluclationSingleProcess for a given site.
   * @param siteIndex the index of the site
   *
   */
  std::shared_ptr<LikelihoodCalculationSingleProcess> getLikelihoodCalculationForASite(size_t siteIndex) const
  {
    return dynamic_pointer_cast<const SingleProcessPhyloLikelihood>(getPhyloLikelihood(vProcPos_[siteIndex].nProc))->getLikelihoodCalculationSingleProcess();
  }

  /**
   * @name The Likelihood interface.
   *
   * @{
   */
  std::shared_ptr<const AlignmentDataInterface> getData() const override
  {
    return getPhyloContainer()->getData(getPhyloContainer()->getNumbersOfPhyloLikelihoods()[0]);
  }

  const std::vector<ProcPos>& getProcessSiteRelations() const
  {
    return vProcPos_;
  }

  LikelihoodCalculation& likelihoodCalculation() const override
  {
    return *likCal_;
  }

  std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation() const override
  {
    return likCal_;
  }

  AlignedLikelihoodCalculation& alignedLikelihoodCalculation() const override
  {
    return *likCal_;
  }

  std::shared_ptr<AlignedLikelihoodCalculation> getAlignedLikelihoodCalculation() const override
  {
    return likCal_;
  }

  size_t getNumberOfSites() const override
  {
    return mSeqEvol_->getNumberOfSites();
  }

  std::shared_ptr<PartitionSequenceEvolution> getPartitionSequenceEvolution() const
  {
    return mSeqEvol_;
  }

private:
  void makeLikCal_();
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PARTITIONPROCESSPHYLOLIKELIHOOD_H
