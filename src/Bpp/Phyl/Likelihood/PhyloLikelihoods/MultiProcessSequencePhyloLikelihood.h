// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_MULTIPROCESSSEQUENCEPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_MULTIPROCESSSEQUENCEPHYLOLIKELIHOOD_H

#include <Bpp/Numeric/AbstractParametrizable.h>

#include "../DataFlow/LikelihoodCalculationSingleProcess.h"
#include "../MultiProcessSequenceEvolution.h"
#include "SequencePhyloLikelihood.h"

// From bpp-seq:
#include <Bpp/Seq/Container/AlignmentData.h>

using namespace std;

namespace bpp
{
/**
 * @brief Partial implementation of the Likelihood interface for
 * multiple processes.
 *
 * This class uses several LikelihoodTreeCalculation instances to
 * compute a the global likelihood of a unique data set, as well as a
 * collection of SubstitutionProcess.
 *
 * It implements the Function interface and manages the parameters of
 * all substitution processes.
 */
class MultiProcessSequencePhyloLikelihood :
  public AbstractParametrizableSequencePhyloLikelihood
{
private:
  /**
   * @brief to avoid the dynamic casts
   */
  std::shared_ptr<MultiProcessSequenceEvolution> mSeqEvol_;

protected:
  /**
   * vector of pointers towards LikelihoodCalculationSingleProcess, used
   * for the global likelihood.
   */
  mutable std::vector<std::shared_ptr<LikelihoodCalculationSingleProcess>> vLikCal_;

public:
  MultiProcessSequencePhyloLikelihood(
      std::shared_ptr<const AlignmentDataInterface> data,
      std::shared_ptr<MultiProcessSequenceEvolution> processSeqEvol,
      std::shared_ptr<CollectionNodes> collNodes,
      size_t nSeqEvol = 0,
      size_t nData = 0);

protected:
  MultiProcessSequencePhyloLikelihood(const MultiProcessSequencePhyloLikelihood& mpspl) :
    AbstractParametrizableSequencePhyloLikelihood(mpspl),
    mSeqEvol_(mpspl.mSeqEvol_),
    vLikCal_(mpspl.vLikCal_)
  {}

  MultiProcessSequencePhyloLikelihood& operator=(const MultiProcessSequencePhyloLikelihood& mpspl)
  {
    AbstractParametrizableSequencePhyloLikelihood::operator=(mpspl);
    mSeqEvol_ = mpspl.mSeqEvol_;
    vLikCal_ = mpspl.vLikCal_;
    return *this;
  }

public:
  virtual ~MultiProcessSequencePhyloLikelihood() {}

public:
  /**
   * @name The Likelihood interface.
   *
   * @{
   */
  std::shared_ptr<const AlignmentDataInterface> getData() const override
  {
    return vLikCal_[0]->getData();
  }

  std::shared_ptr<const Alphabet> getAlphabet() const override
  {
    return vLikCal_[0]->stateMap().getAlphabet();
  }

  /**
   * @}
   */

public:
  std::shared_ptr<LikelihoodCalculationSingleProcess> getLikelihoodCalculationForAProcess(size_t p)
  {
    return vLikCal_[p];
  }

  /**
   * @brief Get the likelihood for a site for a process.
   *
   * @param site The site index to analyse.
   * @param p the process index in the given order.
   * @return The likelihood for site <i>site</i>.
   */
  DataLik getLikelihoodForASiteForAProcess(size_t site, size_t p) const
  {
    return vLikCal_[p]->getLikelihoodForASite(site);
  }

  VVdouble getLikelihoodPerSitePerProcess() const;

  /*
   *@brief return the posterior probabilities of subprocess on each site.
   *
   *@return 2D-vector sites x states
   */

  virtual VVdouble getPosteriorProbabilitiesPerSitePerProcess() const = 0;

  bool isInitialized() const override
  {
    for (auto& lik : vLikCal_)
    {
      if (!lik->isInitialized())
        return false;
    }

    return true;
  }

  /**
   * @brief Set the dataset for which the likelihood must be evaluated.
   *
   * @param sites The data set to use.
   * @param nData the number of the data (optionnal, default 0).
   */

  void setData(std::shared_ptr<const AlignmentDataInterface> sites, size_t nData = 0) override;

  /**
   * @brief Return the number of process used for computation.
   */
  size_t getNumberOfSubstitutionProcess() const { return vLikCal_.size(); }

  /** @} */
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_MULTIPROCESSSEQUENCEPHYLOLIKELIHOOD_H
