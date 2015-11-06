//
// File: AutoCorrelationOfAlignedPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: vendredi 6 novembre 2015, à 00h 12
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

#ifndef _AUTO_CORRELATION_OF_ALIGNED_PHYLOLIKELIHOOD_H_
#define _AUTO_CORRELATION_OF_ALIGNED_PHYLOLIKELIHOOD_H_


#include "SetOfAlignedPhyloLikelihood.h"
#include "HmmPhyloEmissionProbabilities.h"

// From Numeric
#include <Bpp/Numeric/Hmm/HmmLikelihood.h>
#include <Bpp/Numeric/Hmm/LogsumHmmLikelihood.h>
#include <Bpp/Numeric/Hmm/AutoCorrelationTransitionMatrix.h>


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

  
  class AutoCorrelationOfAlignedPhyloLikelihood :
    public SetOfAlignedPhyloLikelihood
  {
  private:
    std::auto_ptr<HmmPhyloAlphabet> hma_;
      
    std::auto_ptr<AutoCorrelationTransitionMatrix> htm_;
    
    std::auto_ptr<HmmPhyloEmissionProbabilities> hpep_;

    mutable std::auto_ptr<LogsumHmmLikelihood> hmm_;
  
  public:
    AutoCorrelationOfAlignedPhyloLikelihood(PhyloLikelihoodContainer* pC, const std::vector<size_t>& nPhylo);

    AutoCorrelationOfAlignedPhyloLikelihood(const AutoCorrelationOfAlignedPhyloLikelihood& mlc) :
      AbstractPhyloLikelihood(mlc),
      AbstractAlignedPhyloLikelihood(mlc),
      SetOfAlignedPhyloLikelihood(mlc),
      hma_(std::auto_ptr<HmmPhyloAlphabet>(mlc.hma_->clone())),
      htm_(std::auto_ptr<AutoCorrelationTransitionMatrix>(mlc.htm_->clone())),
      hpep_(std::auto_ptr<HmmPhyloEmissionProbabilities>(mlc.hpep_->clone())),
      hmm_(std::auto_ptr<LogsumHmmLikelihood>(mlc.hmm_->clone()))
    {}

    AutoCorrelationOfAlignedPhyloLikelihood& operator=(const AutoCorrelationOfAlignedPhyloLikelihood& mlc)
    {
      SetOfAlignedPhyloLikelihood::operator=(mlc);
        
      hma_=std::auto_ptr<HmmPhyloAlphabet>(mlc.hma_->clone());
      htm_=std::auto_ptr<AutoCorrelationTransitionMatrix>(mlc.htm_->clone());
      hpep_=std::auto_ptr<HmmPhyloEmissionProbabilities>(mlc.hpep_->clone());
      hmm_=std::auto_ptr<LogsumHmmLikelihood>(mlc.hmm_->clone());

      return *this;
    }

    virtual ~AutoCorrelationOfAlignedPhyloLikelihood() {}
      
    AutoCorrelationOfAlignedPhyloLikelihood* clone() const { return new AutoCorrelationOfAlignedPhyloLikelihood(*this); }

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

    double getDLogLikelihood() const
    {
      return hmm_->getDLogLikelihood();
    }

    double getD2LogLikelihood() const
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

    double getDLogLikelihoodForASite(size_t site) const
    {
      return hmm_->getDLogLikelihoodForASite(site);
    }

    double getD2LogLikelihoodForASite(size_t site) const
    {
      return hmm_->getD2LogLikelihoodForASite(site);
    }

    Vdouble getLikelihoodForEachSite() const
    {
      return hmm_->getLikelihoodForEachSite();
    }

    VVdouble getPosteriorProbabilitiesForEachSiteForEachAligned() const
    {
      VVdouble pp;
      hmm_->getHiddenStatesPosteriorProbabilities(pp, false);
      return pp;
    }

    const HmmTransitionMatrix& getHmmTransitionMatrix() const
    {
      return hmm_->getHmmTransitionMatrix();
    }
  
  protected:
    void computeDLogLikelihood_(const std::string& variable) const
    {
      hmm_->getFirstOrderDerivative(variable);
    }


    void computeD2LogLikelihood_(const std::string& variable) const
    {
      hmm_->getSecondOrderDerivative(variable);
    }


    /*
     * @}
     */



  };
} // end of namespace bpp.

#endif  // _AUTO_CORRELATION_OF_ALIGNED_PHYLOLIKELIHOOD_H_

