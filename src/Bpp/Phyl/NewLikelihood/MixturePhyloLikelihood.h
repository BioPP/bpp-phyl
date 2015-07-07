//
// File: MixturePhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: jeudi 11 juillet 2013, à 14h 05
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

#ifndef _MIXTUREPHYLOLIKELIHOOD_H_
#define _MIXTUREPHYLOLIKELIHOOD_H_


#include "MultiProcessPhyloLikelihood.h"
#include "MixtureSequenceEvolution.h"

// From SeqLib:
#include <Bpp/Seq/Container/SiteContainer.h>


namespace bpp
{
/**
 * @brief Likelihood framework based on a mixture of simple likelihoods
 *
 * The resulting likelihood is the mean value of
 * the SinglePhyloLikelihoods, ponderated with parametrized probabilities
 * (through a Simplex).
 *
 * @see MultiProcessPhyloLikelihood
 */

    class MixturePhyloLikelihood :
      public MultiProcessPhyloLikelihood
    {
    private:
      /**
       * @brief to avoid the dynamic casts
       *
       */

      MixtureSequenceEvolution& mSeqEvol_;
      
    public:
      MixturePhyloLikelihood(
        const SiteContainer& data,
        MixtureSequenceEvolution& processSeqEvol,
        size_t nSeqEvol = 0,
        size_t nData = 0,
        bool verbose = true,
        bool patterns = true);

      MixturePhyloLikelihood(const MixturePhyloLikelihood& mlc) :
        MultiProcessPhyloLikelihood(mlc),
        mSeqEvol_(mlc.mSeqEvol_)
      {}

      MixturePhyloLikelihood& operator=(const MixturePhyloLikelihood& mlc)
      {
        MultiProcessPhyloLikelihood::operator=(mlc);
        mSeqEvol_=mlc.mSeqEvol_;
        
        return *this;
      }

      virtual ~MixturePhyloLikelihood() {}

      MixturePhyloLikelihood* clone() const { return new MixturePhyloLikelihood(*this); }

    public:
      /**
       * @brief return the probability of a  subprocess
       *
       * @param i the index of the subprocess
       */
  
      double getSubProcessProb(size_t i) const
      {
        return mSeqEvol_.getSubProcessProb(i);
      }

     /**
       * @name The likelihood functions.
       *
       * @{
       */

      double getLogLikelihood() const;

      double getDLogLikelihood() const;
  
      double getD2LogLikelihood() const;

      double getDLogLikelihoodForASite(size_t site) const;
  
      double getD2LogLikelihoodForASite(size_t site) const;
  
      double getLikelihoodForASite(size_t site) const;
  
      VVdouble getPosteriorProbabilitiesForEachSiteForEachProcess() const;

    protected:
  
      void computeDLogLikelihood_(const std::string& variable) const;

      void computeD2LogLikelihood_(const std::string& variable) const;

      /*
       * @}
       */
    };
} // end of namespace bpp.

#endif  // _MIXTURELIKELIHOODCOLLECTION_H_

