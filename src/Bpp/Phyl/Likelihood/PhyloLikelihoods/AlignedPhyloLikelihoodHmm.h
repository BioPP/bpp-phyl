// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMOFALIGNEDPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMOFALIGNEDPHYLOLIKELIHOOD_H


#include "HmmLikelihood_DF.h"
#include "AlignedPhyloLikelihoodSet.h"

// From Numeric
#include <Bpp/Numeric/Hmm/FullHmmTransitionMatrix.h>

#include "HmmPhyloAlphabet.h"
#include "HmmPhyloEmissionProbabilities.h"


namespace bpp
{
/**
 * @brief Likelihood framework based on a hmm of simple likelihoods
 *
 * The resulting likelihood is the likelihood of the given Hmm with
 * the site emission probabilities proportional to the computed
 * likelihoods of the process.
 *
 *
 */
class AlignedPhyloLikelihoodHmm :
  public AbstractAlignedPhyloLikelihoodSet
{
private:
  std::shared_ptr<HmmPhyloAlphabet> hma_;

  std::shared_ptr<FullHmmTransitionMatrix> htm_;

  std::shared_ptr<HmmPhyloEmissionProbabilities> hpep_;

  mutable std::shared_ptr<HmmLikelihood_DF> hmm_;

public:
  AlignedPhyloLikelihoodHmm(
      Context& context,
      std::shared_ptr<PhyloLikelihoodContainer> pC,
      const std::vector<size_t>& nPhylo,
      bool inCollection = true);

  virtual ~AlignedPhyloLikelihoodHmm() {}

protected:
  AlignedPhyloLikelihoodHmm(const AlignedPhyloLikelihoodHmm& mlc) :
    AbstractPhyloLikelihood(mlc),
    AbstractParametrizable(""),
    AbstractPhyloLikelihoodSet(mlc),
    AbstractAlignedPhyloLikelihood(mlc),
    AbstractAlignedPhyloLikelihoodSet(mlc),
    hma_(mlc.hma_),
    htm_(mlc.htm_),
    hpep_(mlc.hpep_),
    hmm_(mlc.hmm_)
  {}

  AlignedPhyloLikelihoodHmm& operator=(const AlignedPhyloLikelihoodHmm& mlc)
  {
    AbstractAlignedPhyloLikelihoodSet::operator=(mlc);
    hma_ = mlc.hma_;
    htm_ = mlc.htm_;
    hpep_ = mlc.hpep_;
    hmm_ = mlc.hmm_;
    return *this;
  }

  AlignedPhyloLikelihoodHmm* clone() const override
  {
    return new AlignedPhyloLikelihoodHmm(*this);
  }

public:
  void setNamespace(const std::string& nameSpace) override;

  void fireParameterChanged(const ParameterList& parameters) override;

  const HmmPhyloAlphabet hmmStateAlphabet() const
  {
    return *hma_;
  }

  std::shared_ptr<const HmmPhyloAlphabet> getHmmStateAlphabet() const
  {
    return hma_;
  }

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
  VVdouble getPosteriorProbabilitiesPerSitePerAligned() const
  {
    VVdouble pp;
    auto mat = hmm_->getHiddenStatesPosteriorProbabilities().transpose();
    copyEigenToBpp(mat, pp);
    return pp;
  }

  Vdouble getPosteriorProbabilitiesForASitePerAligned(size_t site) const
  {
    Vdouble pp;
    auto mat = hmm_->getHiddenStatesPosteriorProbabilitiesForASite(site);
    copyEigenToBpp(mat, pp);
    return pp;
  }

  const Eigen::MatrixXd& getHmmTransitionMatrix() const
  {
    return hmm_->getHmmTransitionMatrix();
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMOFALIGNEDPHYLOLIKELIHOOD_H
