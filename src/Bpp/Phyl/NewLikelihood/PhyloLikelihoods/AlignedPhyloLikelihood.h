//
// File: AlignedPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: jeudi 17 septembre 2015, à 13h 31
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

#ifndef _ALIGNEDPHYLOLIKELIHOOD_H_
#define _ALIGNEDPHYLOLIKELIHOOD_H_

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>

// from bpp-core

#include <Bpp/Numeric/AbstractParametrizable.h>

#include "PhyloLikelihood.h"

namespace bpp
{
  /**
   * @brief The AlignedPhyloLikelihood interface, for phylogenetic likelihood.
   *
   * This interface defines the common methods needed to compute a
   * likelihood from an alignement of data (involving several process
   * & alphabets, if needed), where all sites are independent (eq it
   * does not match for HMM phylolikelihoods).
   *
   */
  
  class AlignedPhyloLikelihood :
    virtual public PhyloLikelihood
  {
  public:
    AlignedPhyloLikelihood() {}
    virtual ~AlignedPhyloLikelihood() {}
    
    virtual AlignedPhyloLikelihood* clone() const = 0;
    
  public:
    
    /**
     *
     * @name The data functions
     *
     * @{
     */
    
    /**
     * @brief Get the number of sites in the dataset.
     *
     * @return the number of sites in the dataset.
     */

    virtual size_t getNumberOfSites() const = 0;

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

    virtual double getLogLikelihoodForASite(size_t site) const = 0;
    
    virtual double getDLogLikelihoodForASite(size_t site) const = 0;
    
    virtual double getD2LogLikelihoodForASite(size_t site) const = 0;
    
    /**
     * @brief Get the likelihood for each site.
     *
     * @return A vector with all likelihoods for each site.
     */
    virtual Vdouble getLikelihoodForEachSite() const = 0;

/** @} */

    friend class SetOfAlignedPhyloLikelihood;
    
  };

      
      
  class AbstractAlignedPhyloLikelihood :
    public virtual AlignedPhyloLikelihood,
    public virtual AbstractPhyloLikelihood
  {
  protected:

    size_t nbSites_;

  public:
    AbstractAlignedPhyloLikelihood(size_t nbSites) :
    AbstractPhyloLikelihood(),
    nbSites_(nbSites)
    {}
      
    AbstractAlignedPhyloLikelihood(const AbstractAlignedPhyloLikelihood& asd) :
    AbstractPhyloLikelihood(asd),
    nbSites_(asd.nbSites_)
    {
    }
      
    virtual ~AbstractAlignedPhyloLikelihood() {}

    AbstractAlignedPhyloLikelihood* clone() const = 0;
      
    AbstractAlignedPhyloLikelihood& operator=(const AbstractAlignedPhyloLikelihood& asd)
    {
      AbstractPhyloLikelihood::operator=(*this);
      
      nbSites_=asd.nbSites_;
        
      return *this;
    }

    size_t getNumberOfSites() const { return nbSites_; }

    void setNumberOfSites(size_t nbSites)
    {
      nbSites_ = nbSites;
    }

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
      
} //end of namespace bpp.

#endif  //_ALIGNEDPHYLOLIKELIHOOD_H_

