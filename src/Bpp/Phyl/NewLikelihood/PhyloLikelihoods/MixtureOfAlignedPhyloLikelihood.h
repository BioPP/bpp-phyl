//
// file: MixtureOfAlignedPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: mercredi 7 octobre 2015, à 14h 03
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

#ifndef _MIXTURE_ALIGNED_PHYLOLIKELIHOOD_H_
#define _MIXTURE_ALIGNED_PHYLOLIKELIHOOD_H_


#include "SetOfAlignedPhyloLikelihood.h"

// From SeqLib:
#include <Bpp/Seq/Container/SiteContainer.h>

#include <Bpp/Numeric/Prob/Simplex.h>

namespace bpp
{
/**
 * @brief Likelihood framework based on a mixture of aligned likelihoods
 *
 * The resulting likelihood is the mean value of
 * the AlignedPhyloLikelihoods, ponderated with parametrized probabilities
 * (through a Simplex).
 *
 */

  class MixtureOfAlignedPhyloLikelihood :
    public SetOfAlignedPhyloLikelihood
  {
  private:
    /**
     * @brief Simplex of the probabilities of the AlignedPhylos, in
     * the same orders as in the
     * SetOfAbstractPhyloLikelihood::nPhylo.
     *
     */

    Simplex simplex_;

  public:
    MixtureOfAlignedPhyloLikelihood(PhyloLikelihoodContainer* pC, const std::vector<size_t>& nPhylo);
      
    MixtureOfAlignedPhyloLikelihood(const MixtureOfAlignedPhyloLikelihood& mlc);

    MixtureOfAlignedPhyloLikelihood& operator=(const MixtureOfAlignedPhyloLikelihood& mlc);

    ~MixtureOfAlignedPhyloLikelihood() {};
      
    MixtureOfAlignedPhyloLikelihood* clone() const
    {
      return new MixtureOfAlignedPhyloLikelihood(*this);
    }

    void setPhyloProb(const Simplex& si);

  protected:
      
    void computeDLogLikelihood_(const std::string& variable) const
    {
      SetOfAbstractPhyloLikelihood::computeDLogLikelihood_(variable);
    }

    void computeD2LogLikelihood_(const std::string& variable) const
    {
      SetOfAbstractPhyloLikelihood::computeD2LogLikelihood_(variable);
    }

    void fireParameterChanged(const ParameterList& parameters);

  public:

    /**
     * @brief return the probability of the phylolikelihood of a given
     * index.
     *
     */

    const std::vector<double>& getPhyloProbabilities() const
    {
      return simplex_.getFrequencies();
    }
    
        
    /**
     * @brief return the probability of the phylolikelihood of a given
     * index.
     *
     */
    
    double getPhyloProb(size_t index)
    {
      return simplex_.prob(index);
    }
    
    /**
     *
     * @name Inherited from PhyloLikelihood
     *
     * @{
     */
      
    /**
     * @name The site likelihood functions.
     *
     */

    double getLogLikelihood() const;

    double getDLogLikelihood() const;

    double getD2LogLikelihood() const;
    
    /**
     * @brief Get the logarithm of the likelihood for any site.
     *
     * @return The logarithm of the likelihood of the data at this site.
     */

    double getLikelihoodForASite(size_t site) const;
    
    double getLogLikelihoodForASite(size_t site) const;

    /**
     * @brief Get the derivates of the LogLikelihood at a Site
     *
     */

    double getDLogLikelihoodForASite(size_t site) const;
      
    double getD2LogLikelihoodForASite(size_t site) const;
      
    /*
     * @}
     */
  };
} // end of namespace bpp.

#endif  // _MIXTURE_ALIGNED_PHYLOLIKELIHOOD_H_

