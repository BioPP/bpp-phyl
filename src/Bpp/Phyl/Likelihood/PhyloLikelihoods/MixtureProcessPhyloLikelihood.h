// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_MIXTUREPROCESSPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_MIXTUREPROCESSPHYLOLIKELIHOOD_H


#include "../DataFlow/Simplex_DF.h"
#include "../MixtureSequenceEvolution.h"
#include "MultiProcessSequencePhyloLikelihood.h"

// From SeqLib:
#include <Bpp/Seq/Container/AlignmentData.h>


namespace bpp
{
/**
 * @brief Likelihood framework based on a mixture of simple likelihoods
 *
 * The resulting likelihood is the mean value of
 * the SinglePhyloLikelihoods, ponderated with parametrized probabilities
 * (through a Simplex).
 *
 * @see MultiProcessSequencePhyloLikelihood
 */

class MixtureProcessPhyloLikelihood :
  public MultiProcessSequencePhyloLikelihood
{
private:
  /**
   * @brief to avoid the dynamic casts
   */
  std::shared_ptr<MixtureSequenceEvolution> mSeqEvol_;

  /**
   * DF simplex for computation
   */
  std::shared_ptr<ConfiguredSimplex> simplex_;

  /**
   * Aligned LikelihoodCalculation to store DF nodes
   */
  mutable std::shared_ptr<AlignedLikelihoodCalculation> likCal_;

public:
  MixtureProcessPhyloLikelihood(
      std::shared_ptr<const AlignmentDataInterface> data,
      std::shared_ptr<MixtureSequenceEvolution> processSeqEvol,
      std::shared_ptr<CollectionNodes> collNodes,
      size_t nSeqEvol = 0,
      size_t nData = 0);

protected:
  MixtureProcessPhyloLikelihood(const MixtureProcessPhyloLikelihood& mlc) = default;

  MixtureProcessPhyloLikelihood& operator=(const MixtureProcessPhyloLikelihood& mlc) = default;

  MixtureProcessPhyloLikelihood* clone() const override
  {
    return new MixtureProcessPhyloLikelihood(*this);
  }

public:
  /**
   * @brief return the probability of a subprocess
   *
   * @param i the index of the subprocess
   */
  double getSubProcessProb(size_t i) const
  {
    return simplex_->targetValue()->prob(i);
  }

  /**
   * @name The likelihood functions.
   *
   * @{
   */
  LikelihoodCalculation& likelihoodCalculation () const override
  {
    return *likCal_;
  }

  std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation () const override
  {
    return likCal_;
  }

  AlignedLikelihoodCalculation& alignedLikelihoodCalculation () const override
  {
    return *likCal_;
  }

  std::shared_ptr<AlignedLikelihoodCalculation> getAlignedLikelihoodCalculation () const override
  {
    return likCal_;
  }

  /**
   * @brief return the posterior probabilities of subprocess on each site.
   *
   * @return 2D-vector sites x states
   */
  VVdouble getPosteriorProbabilitiesPerSitePerProcess() const override;

  /**
   * @}
   */
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_MIXTUREPROCESSPHYLOLIKELIHOOD_H
