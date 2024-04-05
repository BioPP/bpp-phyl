// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_AUTOCORRELATIONOFALIGNEDPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_AUTOCORRELATIONOFALIGNEDPHYLOLIKELIHOOD_H


#include "HmmPhyloEmissionProbabilities.h"
#include "AlignedPhyloLikelihoodSet.h"

// From Numeric

#include "HmmLikelihood_DF.h"
#include <Bpp/Numeric/Hmm/AutoCorrelationTransitionMatrix.h>


namespace bpp
{
/**
 * @brief Likelihood framework based on a hmm of simple likelihoods
 *
 * The resulting likelihood is the likelihood of the given Hmm with
 * the site emission probabilities proportional to the computed
 * likelihoods of the process.
 */
class AlignedPhyloLikelihoodAutoCorrelation :
  public AbstractAlignedPhyloLikelihoodSet
{
private:
  std::shared_ptr<HmmPhyloAlphabet> hma_;

  std::shared_ptr<AutoCorrelationTransitionMatrix> htm_;

  std::shared_ptr<HmmPhyloEmissionProbabilities> hpep_;

  mutable std::shared_ptr<HmmLikelihood_DF> hmm_;

public:
  AlignedPhyloLikelihoodAutoCorrelation(
      Context& context,
      std::shared_ptr<PhyloLikelihoodContainer> pC,
      const std::vector<size_t>& nPhylo,
      bool inCollection = true);

protected:
  AlignedPhyloLikelihoodAutoCorrelation(const AlignedPhyloLikelihoodAutoCorrelation& mlc) :
    AbstractPhyloLikelihood(mlc),
    AbstractParametrizable(mlc),
    AbstractPhyloLikelihoodSet(mlc),
    AbstractAlignedPhyloLikelihood(mlc),
    AbstractAlignedPhyloLikelihoodSet(mlc),
    hma_(mlc.hma_),
    htm_(mlc.htm_),
    hpep_(mlc.hpep_),
    hmm_(mlc.hmm_)
  {}

  AlignedPhyloLikelihoodAutoCorrelation& operator=(const AlignedPhyloLikelihoodAutoCorrelation& mlc)
  {
    AbstractAlignedPhyloLikelihoodSet::operator=(mlc);
    hma_ = mlc.hma_;
    htm_ = mlc.htm_;
    hpep_ = mlc.hpep_;
    hmm_ = mlc.hmm_;
    return *this;
  }

  AlignedPhyloLikelihoodAutoCorrelation* clone() const
  {
    return new AlignedPhyloLikelihoodAutoCorrelation(*this);
  }

public:
  virtual ~AlignedPhyloLikelihoodAutoCorrelation() {}

public:
  void setNamespace(const std::string& nameSpace);

  void fireParameterChanged(const ParameterList& parameters);

  /**
   * @name The likelihood functions.
   *
   * @{
   */
  LikelihoodCalculation& likelihoodCalculation () const
  {
    return *hmm_;
  }

  std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation () const
  {
    return hmm_;
  }

  AlignedLikelihoodCalculation& alignedLikelihoodCalculation () const
  {
    return *hmm_;
  }

  std::shared_ptr<AlignedLikelihoodCalculation> getAlignedLikelihoodCalculation () const
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
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_AUTOCORRELATIONOFALIGNEDPHYLOLIKELIHOOD_H
