//
// File: SingleDataPhyloLikelihood.h
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

#ifndef _SINGLEDATAPHYLOLIKELIHOOD_H_
#define _SINGLEDATAPHYLOLIKELIHOOD_H_

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/SiteContainer.h>

// from bpp-core

#include <Bpp/Numeric/AbstractParametrizable.h>

#include "PhyloLikelihood.h"

namespace bpp
{
    /**
     * @brief The SingleDataPhyloLikelihood interface, for phylogenetic likelihood.
     *
     * This interface defines the common methods needed to compute a likelihood
     * from a sequence alignement, usually involving one or more phylogenetic trees.
     */
    class SingleDataPhyloLikelihood :
      virtual public PhyloLikelihood
    {
    public:
      SingleDataPhyloLikelihood() {}
      virtual ~SingleDataPhyloLikelihood() {}

      virtual SingleDataPhyloLikelihood* clone() const = 0;

    public:

      /**
       *
       * @name The data functions
       *
       * @{
       */
      
      /**
       * @brief Set the dataset for which the likelihood must be evaluated.
       *
       * @param sites The data set to use.
       * @param nData the number of the data
       */
      
      virtual void setData(const SiteContainer& sites, size_t nData = 0) = 0;
    
      /**
       * @brief Get the dataset for which the likelihood must be evaluated.
       *
       * @return A pointer toward the site container where the sequences are stored.
       */
      virtual const SiteContainer* getData() const = 0;

      /**
       * @brief Get the number of dataset concerned.
       *
       */
      virtual size_t getNData() const = 0;

      /**
       * @brief Get the number of sites in the dataset.
       *
       * @return the number of sites in the dataset.
       */
      virtual size_t getNumberOfSites() const = 0;

      /**
       * @brief Get the number the states.
       *
       */

      virtual size_t getNumberOfStates() const = 0;

      /**
       * @brief Get the alphabet associated to the dataset.
       *
       * @return the alphabet associated to the dataset.
       */    
      virtual const Alphabet* getAlphabet() const = 0;
 
      /**
       * @return the recurvisity used for the computations.
       *
       */    
      virtual char getRecursivity() const = 0;
    
      /**
       * @}
       */
      
      /**
       * @name The likelihood functions.
       *
       * @{
       */

      /**
       * @brief Get the likelihood for a site.
       *
       * @param site The site index to analyse.
       * @return The likelihood for site <i>site</i>.
       */
      virtual double getLikelihoodForASite(size_t site) const = 0;

      /**
       * @brief Get the log likelihood for a site, and its derivatives.
       *
       * @param site The site index to analyse.
       * @return The (D)log likelihood for site <i>site</i>.
       */
      
      virtual double getLogLikelihoodForASite(size_t site) const
      {
        return log(getLikelihoodForASite(site));
      }

      virtual double getDLogLikelihoodForASite(size_t site) const = 0;

      virtual double getD2LogLikelihoodForASite(size_t site) const = 0;

/**
       * @brief Get the likelihood for each site.
       *
       * @return A vector with all likelihoods for each site.
       */
      virtual Vdouble getLikelihoodForEachSite() const = 0;

      /** @} */
    };

    
    class AbstractSingleDataPhyloLikelihood :
      public virtual SingleDataPhyloLikelihood
    {
    protected:

      size_t nbSites_;
      size_t nbStates_;

      /**
       * @brief Number of the concerned data.
       *
       **/
      
      size_t nData_;

    public:
      AbstractSingleDataPhyloLikelihood(size_t nbSites, size_t nbStates, size_t nData = 0) :
        nbSites_(nbSites),
        nbStates_(nbStates),
        nData_(nData)
      {}

      
      AbstractSingleDataPhyloLikelihood(const AbstractSingleDataPhyloLikelihood& asd) :
        nbSites_(asd.nbSites_),
        nbStates_(asd.nbStates_),
        nData_(asd.nData_)
      {
      }
      
      virtual ~AbstractSingleDataPhyloLikelihood() {}

      AbstractSingleDataPhyloLikelihood* clone() const = 0;
      
      AbstractSingleDataPhyloLikelihood& operator=(const AbstractSingleDataPhyloLikelihood& asd)
      {
        nbSites_=asd.nbSites_;
        nbStates_=asd.nbStates_;
        
        nData_=asd.nData_;
        
        return *this;
      }

      void setData(const SiteContainer& sites, size_t nData = 0)
      {
        nbSites_ = sites.getNumberOfSites();
        nbStates_ = sites.getAlphabet()->getSize();
        nData_=nData;
        initialize();
      }

      size_t getNData() const
      {
        return nData_;
      }

      void setNData(size_t nData)
      {
        nData_=nData;
      }

      size_t getNumberOfSites() const { return nbSites_; }

      size_t getNumberOfStates() const { return nbStates_; }

      Vdouble getLikelihoodForEachSite() const
      {
        Vdouble l(getNumberOfSites());
        for (unsigned int i = 0; i < l.size(); ++i)
        {
          l[i] = getLikelihoodForASite(i);
        }
        return l;
      }
      
    };
      
    
    class AbstractParametrizableSingleDataPhyloLikelihood :
      public AbstractSingleDataPhyloLikelihood,
      public AbstractParametrizable
    {
    public:
      AbstractParametrizableSingleDataPhyloLikelihood(size_t nData = 0) :
        AbstractSingleDataPhyloLikelihood(0, 0, nData),
        AbstractParametrizable("")
      {}

      AbstractParametrizableSingleDataPhyloLikelihood(const AbstractParametrizableSingleDataPhyloLikelihood& asd) :
        AbstractSingleDataPhyloLikelihood(asd),
        AbstractParametrizable(asd)
      {
      }
      
      virtual ~AbstractParametrizableSingleDataPhyloLikelihood() {}

      AbstractParametrizableSingleDataPhyloLikelihood* clone() const = 0;
      
      AbstractParametrizableSingleDataPhyloLikelihood& operator=(const AbstractParametrizableSingleDataPhyloLikelihood& asd)
      {
        AbstractSingleDataPhyloLikelihood::operator=(*this);
        AbstractParametrizable::operator=(*this);
        
        return *this;
      }

      
    };
      
} //end of namespace bpp.

#endif  //_SINGLEDATAPHYLOLIKELIHOOD_H_

