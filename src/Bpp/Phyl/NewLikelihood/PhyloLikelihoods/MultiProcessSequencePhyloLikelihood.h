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

#ifndef _MULTIPROCESSSEQUENCEPHYLOLIKELIHOOD_H_
#define _MULTIPROCESSSEQUENCEPHYLOLIKELIHOOD_H_

#include "SequencePhyloLikelihood.h"
#include "../MultiProcessSequenceEvolution.h"

#include "../LikelihoodTreeCalculation.h"

#include <Bpp/Numeric/AbstractParametrizable.h>

// From SeqLib:
#include <Bpp/Seq/Container/SiteContainer.h>

using namespace std;

namespace bpp
{
/**
 * @brief Partial implementation of the Likelihood interface for multiple processes.
 *
 * This class uses several LikelihoodTreeCalculation instances to
 * compute a the global likelihood of the data set, as well as a
 * collection of SubstitutionProcess.
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
       * vector of pointers towards TreelikelihoodCalculations, used
       * for the global likelihood.
       */
  
      std::vector<LikelihoodTreeCalculation*> vpTreelik_;

    public:
      MultiProcessSequencePhyloLikelihood(
        const SiteContainer& data,
        MultiProcessSequenceEvolution& processSeqEvol,
        size_t nSeqEvol = 0, 
        size_t nData = 0,
        bool verbose = true,
        bool patterns = true);

      MultiProcessSequencePhyloLikelihood(const MultiProcessSequencePhyloLikelihood& lik) :
        AbstractPhyloLikelihood(lik),
        AbstractAlignedPhyloLikelihood(lik),
        AbstractSequencePhyloLikelihood(lik),
        mSeqEvol_(lik.mSeqEvol_),
        vpTreelik_()
      {
        for (size_t i = 0; i < lik.vpTreelik_.size(); i++)
        {
          vpTreelik_.push_back(lik.vpTreelik_[i]->clone());
        }
      }

      MultiProcessSequencePhyloLikelihood& operator=(const MultiProcessSequencePhyloLikelihood& lik)
      {
        AbstractSequencePhyloLikelihood::operator=(lik);
        mSeqEvol_=lik.mSeqEvol_;

        for (size_t i = 0; i < vpTreelik_.size(); i++)
        {
          if (vpTreelik_[i])
            delete vpTreelik_[i];
        }

        vpTreelik_.empty();

        for (size_t i = 0; i < lik.vpTreelik_.size(); i++)
        {
          vpTreelik_.push_back(lik.vpTreelik_[i]->clone());
        }

        return *this;
      }

      /**
       * @brief Abstract class destructor
       *
       */
      virtual ~MultiProcessSequencePhyloLikelihood()
      {
        for (size_t i = 0; i < vpTreelik_.size(); i++)
        {
          if (vpTreelik_[i])
            delete vpTreelik_[i];
        }
        vpTreelik_.empty();
      }

    public:
      /**
       * @name The Likelihood interface.
       *
       * @{
       */

      const SiteContainer* getData() const {
        return vpTreelik_[0]->getData();
      }

      const Alphabet* getAlphabet() const {
        return vpTreelik_[0]->getAlphabet();
      }

      /*
       * @}
       */
      
    public:
      /**
       * @brief To be defined in inheriting classes.
       *
       */

      void updateLikelihood() const
      {
        if (computeLikelihoods_)
        {
          for (size_t i = 0; i < vpTreelik_.size(); i++)
            vpTreelik_[i]->updateLikelihood();
        }
      }
      
      void computeLikelihood() const
      {
        if (computeLikelihoods_)
        {
          for (size_t i = 0; i < vpTreelik_.size(); i++)
            vpTreelik_[i]->computeTreeLikelihood();

          computeLikelihoods_=false;
        }
        
      }

      /**
       * @brief sets using log in all likelihood arrays.
       *
       */
      
      void setUseLog(bool useLog)
      {
        for (size_t i = 0; i < vpTreelik_.size(); i++)
          vpTreelik_[i]->setAllUseLog(useLog);
      }

      /**
       * @brief Get the likelihood for a site.
       *
       * @param site The site index to analyse.
       * @return The likelihood for site <i>site</i>.
       */

      virtual double getLikelihoodForASite(size_t site) const = 0;

      virtual double getLogLikelihoodForASite(size_t site) const = 0;

      /**
       * @brief Get the likelihood for a site for a process.
       *
       * !!! No check on the update of the computing
       *
       * @param site The site index to analyse.
       * @param p the process index in the given order.
       * @return The likelihood for site <i>site</i>.
       */

      double getLikelihoodForASiteForAProcess(size_t site, size_t p) const
      {
        return vpTreelik_[p]->getLikelihoodForASite(site);
      }

      /**
       * @brief Compute the first derivative of the likelihood for a process.
       *
       * @param variable the name of the variable.
       * @param p the process index in the given order.
       * @return The likelihood for site <i>site</i>.
       */

      void computeDLogLikelihoodForAProcess(std::string& variable, size_t p) const;

      /**
       * @brief Get the first derivative of the likelihood for a site for
       * a process.
       *
       * @param site The site index to analyse.
       * @param p the process index in the given order.
       * @return The likelihood for site <i>site</i>.
       */

      virtual double getDLogLikelihoodForASiteForAProcess(size_t site, size_t p) const
      {
        return vpTreelik_[p]->getDLogLikelihoodForASite(site);
      }

      /**
       * @brief Compute the second derivative of the likelihood for a process.
       *
       * @param variable the name of the variable.
       * @param p the process index in the given order.
       * @return The likelihood for site <i>site</i>.
       */
  
      virtual void computeD2LogLikelihoodForAProcess(std::string& variable, size_t p) const;

      /**
       * @brief Get the second derivative of the likelihood for a site for
       * a process.
       *
       * @param site The site index to analyse.
       * @param p the process index in the given order.
       * @return The likelihood for site <i>site</i>.
       */

      double getD2LogLikelihoodForASiteForAProcess(size_t site, size_t p) const
      {
        return vpTreelik_[p]->getD2LogLikelihoodForASite(site);
      }

      VVdouble getLikelihoodForEachSiteForEachProcess() const;

      virtual VVdouble getPosteriorProbabilitiesForEachSiteForEachProcess() const = 0;

  
      /**
       * @brief Set the dataset for which the likelihood must be evaluated.
       *
       * @param sites The data set to use.
       * @param nData the number of the data (optionnal, default 0).
       */
  
      void setData(const SiteContainer& sites, size_t nData = 0);

      /**
       * @brief Return the number of process used for computation.
       */
  
      size_t getNumberOfSubstitutionProcess() const { return vpTreelik_.size(); }
      
      /** @} */
    };
} // end of namespace bpp.

#endif  // _MULTIPROCESSSEQUENCEPHYLOLIKELIHOOD_H_

