//
// File: PartitionProcessPhyloLikelihood.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: samedi 16 mai 2015, ÃÂ  13h 34
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PARTITIONPROCESSPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PARTITIONPROCESSPHYLOLIKELIHOOD_H


#include "../PartitionSequenceEvolution.h"
#include "ProductOfAlignedPhyloLikelihood.h"
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
  public AbstractSetOfPhyloLikelihood
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
  std::map<size_t, std::shared_ptr<AlignmentDataInterface> > mData_;

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
    AbstractSetOfPhyloLikelihood(lik),
    mSeqEvol_(lik.mSeqEvol_),
    vProcPos_(lik.vProcPos_),
    mData_(lik.mData_),
    likCal_(lik.likCal_)
  {}

  PartitionProcessPhyloLikelihood& operator=(const PartitionProcessPhyloLikelihood& lik)
  {
    AbstractSequencePhyloLikelihood::operator=(lik);
    AbstractSetOfPhyloLikelihood::operator=(lik);
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
   * @param nData the number of the data (optionnal, default = 0)
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

private:
  void makeLikCal_();
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PARTITIONPROCESSPHYLOLIKELIHOOD_H
