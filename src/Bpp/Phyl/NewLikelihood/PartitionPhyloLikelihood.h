//
// File: PartitionPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: samedi 16 mai 2015, à 13h 34
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

#ifndef _PARTITION_PHYLOLIKELIHOOD_H_
#define _PARTITION_PHYLOLIKELIHOOD_H_


#include "SequencePhyloLikelihood.h"
#include "PartitionSequenceEvolution.h"
#include "SumOfDataPhyloLikelihood.h"

// From SeqLib:
#include <Bpp/Seq/Container/SiteContainer.h>

namespace bpp
{
/**
 * @brief Likelihood framework based on a partition of a sequence in
 * simple likelihoods.
 *
 * @see AbstractSequencePhyloLikelihood
 */

    struct ProcPos
    {
      size_t nProc;
      size_t pos;
    };
    
      
    class PartitionPhyloLikelihood :
      public SequencePhyloLikelihood,
      public SumOfDataPhyloLikelihood
    {
    private:
      /**
       * @brief to avoid the dynamic casts
       *
       */

      PartitionSequenceEvolution& mSeqEvol_;

      /**
       * vector of couples number of process, sites specific to
       * this process.
       *
       */

      std::vector<ProcPos> vProcPos_;
      
    public:
      PartitionPhyloLikelihood(
        PartitionSequenceEvolution& processSeqEvol,
        size_t nSeqEvol = 0,
        bool verbose = true,
        bool patterns = true);

      PartitionPhyloLikelihood(
        const SiteContainer& data,
        PartitionSequenceEvolution& processSeqEvol,
        size_t nSeqEvol = 0,
        size_t nData = 0,
        bool verbose = true,
        bool patterns = true);

      PartitionPhyloLikelihood(const PartitionPhyloLikelihood& lik) :
        SequencePhyloLikelihood(lik),
        SumOfDataPhyloLikelihood(lik),
        mSeqEvol_(lik.mSeqEvol_),
        vProcPos_(lik.vProcPos_)
      {
      }

      PartitionPhyloLikelihood& operator=(const PartitionPhyloLikelihood& lik)
      {
        SequencePhyloLikelihood::operator=(lik);
        SumOfDataPhyloLikelihood::operator=(lik);
        
        mSeqEvol_=lik.mSeqEvol_;
        vProcPos_=lik.vProcPos_;
 
        return *this;
      }

      virtual ~PartitionPhyloLikelihood()
      {
      }

      PartitionPhyloLikelihood* clone() const { return new PartitionPhyloLikelihood(*this); }

      /**
       * @brief Set the dataset for which the likelihood must be evaluated.
       *
       * @param data The data set to use.
       * @param nData the number of the data (optionnal, default = 0)
       */
  
      void setData(const SiteContainer& data, size_t nData = 0);

      /**
       * @name The Likelihood interface.
       *
       * @{
       */
      
      const Alphabet* getAlphabet() const {
        return getSingleDataPhylolikelihood(vProcPos_[0].nProc)->getAlphabet();
      }

      const SiteContainer* getData() const
      {
        return SumOfDataPhyloLikelihood::getData((size_t)0);
      }
      
      char getRecursivity() const 
      {
        return getSingleDataPhylolikelihood(vProcPos_[0].nProc)->getRecursivity();
      }
      
      /**
       *
       * @}
       */
      
      /**
       * @name The likelihood functions.
       *
       * @{
       */

      double getLikelihoodForASite(size_t site) const;

      double getDLogLikelihoodForASite(size_t site) const;
  
      double getD2LogLikelihoodForASite(size_t site) const;
  
      /**
       *
       * @}
       */
      
    };
} // end of namespace bpp.

#endif  // _PARTITION_PHYLOLIKELIHOOD_H_

