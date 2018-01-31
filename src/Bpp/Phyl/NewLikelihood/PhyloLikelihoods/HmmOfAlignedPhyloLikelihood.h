//
// File: HmmOfAlignedPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: mercredi 28 octobre 2015, à 23h 06
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef _HMM_OF_ALIGNED_PHYLOLIKELIHOOD_H_
#define _HMM_OF_ALIGNED_PHYLOLIKELIHOOD_H_


#include "SetOfAlignedPhyloLikelihood.h"
#include "HmmPhyloEmissionProbabilities.h"

// From Numeric
#include <Bpp/Numeric/Hmm/HmmLikelihood.h>
#include <Bpp/Numeric/Hmm/LogsumHmmLikelihood.h>
#include <Bpp/Numeric/Hmm/FullHmmTransitionMatrix.h>


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
class HmmOfAlignedPhyloLikelihood :
  public SetOfAlignedPhyloLikelihood
{
private:
  std::unique_ptr<HmmPhyloAlphabet> hma_;

  std::unique_ptr<FullHmmTransitionMatrix> htm_;

  std::unique_ptr<HmmPhyloEmissionProbabilities> hpep_;

  mutable std::unique_ptr<LogsumHmmLikelihood> hmm_;

public:
  HmmOfAlignedPhyloLikelihood(PhyloLikelihoodContainer* pC, const std::vector<size_t>& nPhylo);

  HmmOfAlignedPhyloLikelihood(const HmmOfAlignedPhyloLikelihood& mlc) :
    AbstractPhyloLikelihood(mlc),
    AbstractAlignedPhyloLikelihood(mlc),
    SetOfAlignedPhyloLikelihood(mlc),
    hma_(std::unique_ptr<HmmPhyloAlphabet>(mlc.hma_->clone())),
    htm_(std::unique_ptr<FullHmmTransitionMatrix>(mlc.htm_->clone())),
    hpep_(std::unique_ptr<HmmPhyloEmissionProbabilities>(mlc.hpep_->clone())),
    hmm_(std::unique_ptr<LogsumHmmLikelihood>(mlc.hmm_->clone()))
  {}

  HmmOfAlignedPhyloLikelihood& operator=(const HmmOfAlignedPhyloLikelihood& mlc)
  {
    SetOfAlignedPhyloLikelihood::operator=(mlc);

    hma_ = std::unique_ptr<HmmPhyloAlphabet>(mlc.hma_->clone());
    htm_ = std::unique_ptr<FullHmmTransitionMatrix>(mlc.htm_->clone());
    hpep_ = std::unique_ptr<HmmPhyloEmissionProbabilities>(mlc.hpep_->clone());
    hmm_ = std::unique_ptr<LogsumHmmLikelihood>(mlc.hmm_->clone());

    return *this;
  }

  virtual ~HmmOfAlignedPhyloLikelihood() {}

  HmmOfAlignedPhyloLikelihood* clone() const { return new HmmOfAlignedPhyloLikelihood(*this); }

public:
  void setNamespace(const std::string& nameSpace);

  void fireParameterChanged(const ParameterList& parameters);

  /**
   * @name The likelihood functions.
   *
   * @{
   */
  double getLogLikelihood() const
  {
    return hmm_->getLogLikelihood();
  }

  double getDLogLikelihood(const std::string& variable) const
  {
    return hmm_->getDLogLikelihood();
  }

  double getD2LogLikelihood(const std::string& variable) const
  {
    return hmm_->getD2LogLikelihood();
  }

  /**
   * @brief Get the likelihood for a site.
   *
   * @param site The site index to analyse.
   * @return The likelihood for site <i>site</i>.
   */
  double getLikelihoodForASite(size_t site) const
  {
    return hmm_->getLikelihoodForASite(site);
  }

  double getLogLikelihoodForASite(size_t site) const
  {
    return log(hmm_->getLikelihoodForASite(site));
  }

  double getDLogLikelihoodForASite(const std::string& variable, size_t site) const
  {
    return hmm_->getDLogLikelihoodForASite(site);
  }

  double getD2LogLikelihoodForASite(const std::string& variable, size_t site) const
  {
    return hmm_->getD2LogLikelihoodForASite(site);
  }

  Vdouble getLikelihoodPerSite() const
  {
    return hmm_->getLikelihoodForEachSite();
  }

  VVdouble getPosteriorProbabilitiesPerSitePerAligned() const
  {
    VVdouble pp;
    hmm_->getHiddenStatesPosteriorProbabilities(pp, false);
    return pp;
  }

  Vdouble getPosteriorProbabilitiesForASitePerAligned(size_t site) const
  {
    return hmm_->getHiddenStatesPosteriorProbabilitiesForASite(site);
  }

  const HmmTransitionMatrix& getHmmTransitionMatrix() const
  {
    return hmm_->getHmmTransitionMatrix();
  }

protected:
  void computeDLogLikelihood_(const std::string& variable) const
  {
    hmm_->getFirstOrderDerivative(variable);

    dValues_[variable]= std::nan("");
  }


  void computeD2LogLikelihood_(const std::string& variable) const
  {
    hmm_->getSecondOrderDerivative(variable);

    d2Values_[variable]= std::nan("");
  }

  ParameterList getNonDerivableParameters() const;
  
  /*
   * @}
   */
};
} // end of namespace bpp.

#endif  // _HMM_OF_ALIGNED_PHYLOLIKELIHOOD_H_
