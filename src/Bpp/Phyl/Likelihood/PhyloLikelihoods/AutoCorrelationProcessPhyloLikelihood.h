//
// File: AutoCorrelationProcessPhyloLikelihood.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: lundi 23 septembre 2013, ÃÂ  22h 56
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
