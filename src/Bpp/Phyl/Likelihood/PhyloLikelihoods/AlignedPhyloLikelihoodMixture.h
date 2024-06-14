// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_MIXTUREOFALIGNEDPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_MIXTUREOFALIGNEDPHYLOLIKELIHOOD_H


#include "AlignedPhyloLikelihoodSet.h"

// From SeqLib:
#include <Bpp/Seq/Container/AlignmentData.h>

#include "../DataFlow/Simplex_DF.h"

namespace bpp
{
/**
 * @brief Likelihood framework based on a mixture of aligned likelihoods
 *
 * The resulting likelihood is the mean value of
 * the AlignedPhyloLikelihoods, ponderated with parametrized probabilities
 * (through a Simplex).
 */
class AlignedPhyloLikelihoodMixture :
  public AbstractAlignedPhyloLikelihoodSet
{
private:
  /**
   * DF simplex
   */
  std::shared_ptr<ConfiguredSimplex> simplex_;

  /**
   * Aligned LikelihoodCalculation to store DF nodes
   */
  mutable std::shared_ptr<AlignedLikelihoodCalculation> likCal_;

public:
  AlignedPhyloLikelihoodMixture(
      Context& context,
      std::shared_ptr<PhyloLikelihoodContainer> pC,
      const std::vector<size_t>& nPhylo,
      bool inCollection = true);

  virtual ~AlignedPhyloLikelihoodMixture() {}

protected:
  AlignedPhyloLikelihoodMixture(const AlignedPhyloLikelihoodMixture& mlc) :
    AbstractPhyloLikelihood(mlc),
    AbstractParametrizable(""),
    AbstractPhyloLikelihoodSet(mlc),
    AbstractAlignedPhyloLikelihood(mlc),
    AbstractAlignedPhyloLikelihoodSet(mlc),
    simplex_(mlc.simplex_),
    likCal_(mlc.likCal_)
  {}

  AlignedPhyloLikelihoodMixture& operator=(const AlignedPhyloLikelihoodMixture& mlc)
  {
    AbstractAlignedPhyloLikelihoodSet::operator=(mlc);
    simplex_ = mlc.simplex_;
    likCal_ = mlc.likCal_;
    return *this;
  }

  AlignedPhyloLikelihoodMixture* clone() const override
  {
    return new AlignedPhyloLikelihoodMixture(*this);
  }

protected:
  void fireParameterChanged(const ParameterList& parameters) override;

public:
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
   * @brief Get the probabilities of the simplex
   */
  Vdouble getPhyloProbabilities() const;

  /**
   * @brief Get the probability of a phylolikelihood
   */
  double getPhyloProb(size_t index) const;

  /**
   * @brief Set the probabilities of the simplex
   */
  void setPhyloProb(Simplex const& simplex);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_MIXTUREOFALIGNEDPHYLOLIKELIHOOD_H
