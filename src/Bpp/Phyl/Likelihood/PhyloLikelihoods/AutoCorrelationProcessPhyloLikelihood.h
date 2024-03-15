// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_AUTOCORRELATIONPROCESSPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_AUTOCORRELATIONPROCESSPHYLOLIKELIHOOD_H


#include "../AutoCorrelationSequenceEvolution.h"
#include "HmmPhyloEmissionProbabilities.h"
#include "MultiProcessSequencePhyloLikelihood.h"

// From bpp-seq:
#include <Bpp/Seq/Container/AlignmentData.h>

#include "HmmLikelihood_DF.h"
#include <Bpp/Numeric/Hmm/AutoCorrelationTransitionMatrix.h>

namespace bpp
{
/**
 * @brief Likelihood framework based on an auto-correlation of simple likelihoods.
 *
 * The resulting likelihood is the likelihood of the given
 * AutoCorrelation with the site emission probabilities proportional
 * to the computed likelihoods of the process.
 */

class AutoCorrelationProcessPhyloLikelihood :
  public MultiProcessSequencePhyloLikelihood
{
private:
  std::shared_ptr<HmmPhyloEmissionProbabilities> Hpep_;

  mutable std::shared_ptr<HmmLikelihood_DF> hmm_;

public:
  AutoCorrelationProcessPhyloLikelihood(
    std::shared_ptr<const AlignmentDataInterface> data,
    std::shared_ptr<AutoCorrelationSequenceEvolution> processSeqEvol,
    std::shared_ptr<CollectionNodes> collNodes,
    size_t nSeqEvol = 0,
    size_t nData = 0);

protected:
  AutoCorrelationProcessPhyloLikelihood(const AutoCorrelationProcessPhyloLikelihood& mlc) :
    AbstractPhyloLikelihood(mlc),
    AbstractAlignedPhyloLikelihood(mlc),
    AbstractSingleDataPhyloLikelihood(mlc),
    AbstractSequencePhyloLikelihood(mlc),
    AbstractParametrizable(""),
    MultiProcessSequencePhyloLikelihood(mlc),
    Hpep_(mlc.Hpep_),
    hmm_(mlc.hmm_)
  {}

  AutoCorrelationProcessPhyloLikelihood& operator=(const AutoCorrelationProcessPhyloLikelihood& mlc)
  { 
    MultiProcessSequencePhyloLikelihood::operator=(mlc);
    Hpep_ = mlc.Hpep_;
    hmm_ = mlc.hmm_;
    return *this;
  }

public:
  virtual ~AutoCorrelationProcessPhyloLikelihood() {}

  AutoCorrelationProcessPhyloLikelihood* clone() const override
  { 
    return new AutoCorrelationProcessPhyloLikelihood(*this);
  }

public:
  void setNamespace(const std::string& nameSpace) override;

  void fireParameterChanged(const ParameterList& parameters) override;

  /**
   * @name The likelihood functions.
   *
   * @{
   */
  LikelihoodCalculation& likelihoodCalculation () const override
  {
    return *hmm_;
  }
  
  std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation () const override
  {
    return hmm_;
  }

  AlignedLikelihoodCalculation& alignedLikelihoodCalculation () const override
  {
    return *hmm_;
  }

  std::shared_ptr<AlignedLikelihoodCalculation> getAlignedLikelihoodCalculation () const override
  {
    return hmm_;
  }

  /**
   * @brief return the posterior probabilities of subprocess on each site.
   *
   * @return MatrixXd sites x states
   */
  VVdouble getPosteriorProbabilitiesPerSitePerProcess() const override
  {
    VVdouble pp;
    auto mat = hmm_->getHiddenStatesPosteriorProbabilities().transpose();
    copyEigenToBpp(mat, pp);
    return pp;
  }

  const Eigen::MatrixXd& getHmmTransitionMatrix() const
  {
    return hmm_->getHmmTransitionMatrix();
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_AUTOCORRELATIONPROCESSPHYLOLIKELIHOOD_H
