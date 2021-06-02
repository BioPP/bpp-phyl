//
// File: MultiProcessSequencePhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: jeudi 11 juillet 2013, à 21h 51
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

#ifndef _MULTI_PROCESSSEQUENCEPHYLOLIKELIHOOD_H_
#define _MULTI_PROCESSSEQUENCEPHYLOLIKELIHOOD_H_

#include "SequencePhyloLikelihood.h"
#include "../MultiProcessSequenceEvolution.h"

#include "../DataFlow/LikelihoodCalculationSingleProcess.h"

#include <Bpp/Numeric/AbstractParametrizable.h>

// From SeqLib:
#include <Bpp/Seq/Container/AlignedValuesContainer.h>

using namespace std;

namespace bpp
{
/**
 * @brief Partial implementation of the Likelihood interface for
 * multiple processes.
 *
 * This class uses several LikelihoodTreeCalculation instances to
 * compute a the global likelihood of a unique data set, as well as a
 * collection of SubstitutionProcess.
 *
 * It implements the Function interface and manages the parameters of
 * all substitution processes.
 */

    class MultiProcessSequencePhyloLikelihood :
      public AbstractSequencePhyloLikelihood
    {
    private:
      /**
       * @brief to avoid the dynamic casts
       *
       */

      MultiProcessSequenceEvolution& mSeqEvol_;

    protected:
      /**
       * vector of pointers towards LikelihoodCalculationSingleProcess, used
       * for the global likelihood.
       */
  
      mutable std::vector<std::shared_ptr<LikelihoodCalculationSingleProcess>> vLikCal_;

    public:
      MultiProcessSequencePhyloLikelihood(
        const AlignedValuesContainer& data,
        MultiProcessSequenceEvolution& processSeqEvol,
        CollectionNodes& collNodes,
        size_t nSeqEvol = 0, 
        size_t nData = 0,
        bool verbose = true,
        bool patterns = true);

      MultiProcessSequencePhyloLikelihood(const MultiProcessSequencePhyloLikelihood& lik) :
        AbstractPhyloLikelihood(lik),
        AbstractAlignedPhyloLikelihood(lik),
        AbstractSequencePhyloLikelihood(lik),
        mSeqEvol_(lik.mSeqEvol_),
        vLikCal_(lik.vLikCal_)
      {
      }

    public:
      /**
       * @name The Likelihood interface.
       *
       * @{
       */

      const AlignedValuesContainer* getData() const {
        return vLikCal_[0]->getData();
      }

      const Alphabet* getAlphabet() const {
        return vLikCal_[0]->getStateMap().getAlphabet();
      }

      /*
       * @}
       */
      
    public:

      /**
       * @brief Get the likelihood for a site for a process.
       *
       * @param site The site index to analyse.
       * @param p the process index in the given order.
       * @return The likelihood for site <i>site</i>.
       */

      std::shared_ptr<LikelihoodCalculationSingleProcess> getLikelihoodCalculationForAProcess(size_t p)
      {
        return vLikCal_[p];
      }
      
      double getLikelihoodForASiteForAProcess(size_t site, size_t p) const
      {
        return vLikCal_[p]->getLikelihoodForASite(site);
      }

      VVdouble getLikelihoodPerSitePerProcess() const;

      /*
       *@brief return the posterior probabilities of subprocess on each site.
       *
       *@return 2D-vector sites x states
       */

      virtual VVdouble getPosteriorProbabilitiesPerSitePerProcess() const = 0;

      bool isInitialized() const 
      {
        for (auto& lik : vLikCal_)
          if (!lik->isInitialized())
            return false;
      
        return true;
      }
  
      /**
       * @brief Set the dataset for which the likelihood must be evaluated.
       *
       * @param sites The data set to use.
       * @param nData the number of the data (optionnal, default 0).
       */
  
      void setData(const AlignedValuesContainer& sites, size_t nData = 0);

      /**
       * @brief Return the number of process used for computation.
       */
  
      size_t getNumberOfSubstitutionProcess() const { return vLikCal_.size(); }

      /** @} */
    };
} // end of namespace bpp.

#endif  // _MULTIPROCESSSEQUENCEPHYLOLIKELIHOOD_H_

