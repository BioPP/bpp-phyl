//
// File: AutoCorrelationPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: lundi 23 septembre 2013, à 22h 56
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

#ifndef _AUTOCORRELATIONPHYLOLIKELIHOOD_H_
#define _AUTOCORRELATIONPHYLOLIKELIHOOD_H_


#include "MultiProcessPhyloLikelihood.h"

// From SeqLib:
#include <Bpp/Seq/Container/SiteContainer.h>

// From Numeric
#include <Bpp/Numeric/Hmm/HmmLikelihood.h>
#include <Bpp/Numeric/Hmm/LogsumHmmLikelihood.h>
#include <Bpp/Numeric/Hmm/AutoCorrelationTransitionMatrix.h>



namespace bpp
{
namespace newlik
{
/**
 * @brief Likelihood framework based on an auto-correlation of simple likelihoods.
 *
 * The resulting likelihood is the likelihood of the given
 * AutoCorrelation with the site emission probabilities proportional
 * to the computed likelihoods of the process.
 *
 *
 */

  
class AutoCorrelationPhyloLikelihood :
    public MultiProcessPhyloLikelihood
{
private:
  std::auto_ptr<LogsumHmmLikelihood> Hmm_;
  
public:
  AutoCorrelationPhyloLikelihood(
    SubstitutionProcessCollection* processColl,
    char recursivity,
    bool verbose = true,
    bool patterns = true);

  AutoCorrelationPhyloLikelihood(
    const SiteContainer& data,
    SubstitutionProcessCollection* processColl,
    char recursivity,
    size_t nData = 0,
    bool verbose = true,
    bool patterns = true);

  AutoCorrelationPhyloLikelihood(const AutoCorrelationPhyloLikelihood& mlc) :
    MultiProcessPhyloLikelihood(mlc),
    Hmm_(std::auto_ptr<LogsumHmmLikelihood>(mlc.Hmm_->clone())) {}

  AutoCorrelationPhyloLikelihood& operator=(const AutoCorrelationPhyloLikelihood& mlc)
  {
    MultiProcessPhyloLikelihood::operator=(mlc);
    Hmm_ = std::auto_ptr<LogsumHmmLikelihood>(mlc.Hmm_->clone());
    return *this;
  }

  virtual ~AutoCorrelationPhyloLikelihood() {}

  AutoCorrelationPhyloLikelihood* clone() const { return new AutoCorrelationPhyloLikelihood(*this); }

public:
  void fireParameterChanged(const ParameterList& parameters);

  ParameterList getDerivableParameters() const { return getBranchLengthsParameters(); }

  ParameterList getNonDerivableParameters() const;

  /**
   * @name The likelihood functions.
   *
   * @{
   */

  double getLogLikelihood() const
  {
    return Hmm_->getLogLikelihood();  
  }

  double getDLogLikelihood() const
  {
    return Hmm_->getDLogLikelihood();
  }

  double getD2LogLikelihood() const
  {
    return Hmm_->getD2LogLikelihood();  
  }

  /**
   * @brief Get the likelihood for a site.
   *
   * @param site The site index to analyse.
   * @return The likelihood for site <i>site</i>.
   */

  double getLikelihoodForASite(size_t site) const
  {
    return Hmm_->getLikelihoodForASite(site);
  }

  Vdouble getLikelihoodForEachSite() const
  {
    return Hmm_->getLikelihoodForEachSite();
  }

  VVdouble getPosteriorProbabilitiesForEachSiteForEachProcess() const
  {
    VVdouble pp;
    Hmm_->getHiddenStatesPosteriorProbabilities(pp, false);
    return pp;
  }

  const HmmTransitionMatrix& getHmmTransitionMatrix() const
  {
    return Hmm_->getHmmTransitionMatrix();
  }
  
protected:
  void computeDLogLikelihood_(const std::string& variable) const
  {
    Hmm_->getFirstOrderDerivative(variable);
  }


  void computeD2LogLikelihood_(const std::string& variable) const
  {
    Hmm_->getSecondOrderDerivative(variable);
  }


  /*
   * @}
   */



};
} // end of namespace newlik.
} // end of namespace bpp.

#endif  // _AUTOCORRELATIONPHYLOLIKELIHOOD_H_

