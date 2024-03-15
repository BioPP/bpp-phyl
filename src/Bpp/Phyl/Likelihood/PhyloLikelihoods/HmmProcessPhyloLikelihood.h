// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMPROCESSPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMPROCESSPHYLOLIKELIHOOD_H


#include "../HmmSequenceEvolution.h"
#include "HmmPhyloEmissionProbabilities.h"
#include "MultiProcessSequencePhyloLikelihood.h"

// From SeqLib:
#include <Bpp/Seq/Container/AlignmentData.h>


#include "HmmLikelihood_DF.h"
#include <Bpp/Numeric/Hmm/HmmTransitionMatrix.h>


namespace bpp
{
/**
 * @brief Likelihood framework based on a hmm of simple likelihoods
 *
 * The resulting likelihood is the likelihood of the given Hmm with
 * the site emission probabilities proportional to the computed
 * likelihoods of the process.
 */
class HmmProcessPhyloLikelihood :
  public MultiProcessSequencePhyloLikelihood
{
private:
  std::shared_ptr<HmmPhyloEmissionProbabilities> Hpep_;

  /**
   * @brief LikelihoodCalculation in context of HMM.
   */
  mutable std::shared_ptr<HmmLikelihood_DF> hmm_;

public:
  HmmProcessPhyloLikelihood(
    std::shared_ptr<const AlignmentDataInterface> data,
    std::shared_ptr<HmmSequenceEvolution> processSeqEvol,
    std::shared_ptr<CollectionNodes> collNodes,
    size_t nSeqEvol = 0,
    size_t nData = 0);

protected:
  HmmProcessPhyloLikelihood(const HmmProcessPhyloLikelihood& mlc) :
    AbstractPhyloLikelihood(mlc),
    AbstractAlignedPhyloLikelihood(mlc),
    AbstractSingleDataPhyloLikelihood(mlc),
    AbstractSequencePhyloLikelihood(mlc),
    AbstractParametrizable(mlc),
    MultiProcessSequencePhyloLikelihood(mlc),
    Hpep_(std::shared_ptr<HmmPhyloEmissionProbabilities>(mlc.Hpep_->clone())),
    hmm_(std::shared_ptr<HmmLikelihood_DF>(mlc.hmm_->clone())) {}

  HmmProcessPhyloLikelihood* clone() const override { return new HmmProcessPhyloLikelihood(*this); }

public:

  virtual ~HmmProcessPhyloLikelihood() {}

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

  /*
   *@brief return the posterior probabilities of subprocess on each site.
   *
   *@return MatrixXd sites x states
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
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMPROCESSPHYLOLIKELIHOOD_H
