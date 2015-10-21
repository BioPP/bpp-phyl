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

#include "AlignedPhyloLikelihood.h"

namespace bpp
{
    /**
     * @brief The SingleDataPhyloLikelihood interface, for phylogenetic likelihood.
     *
     * This interface defines the common methods needed to compute a
     * likelihood from aligned sequences.
     *
     */
  
  class SingleDataPhyloLikelihood :
    virtual public AlignedPhyloLikelihood
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
     * @}
     */
    
  };

  
  class AbstractSingleDataPhyloLikelihood :
    public SingleDataPhyloLikelihood,
    virtual public AbstractAlignedPhyloLikelihood
  {
  protected:
    size_t nbStates_;

    /**
     * @brief Number of the concerned data.
     *
     **/
    
    size_t nData_;
    
  public:
    AbstractSingleDataPhyloLikelihood(size_t nbSites, size_t nbStates, size_t nData = 0) :
      AbstractAlignedPhyloLikelihood(nbSites),
      nbStates_(nbStates),
      nData_(nData)
    {}
    
    
    AbstractSingleDataPhyloLikelihood(const AbstractSingleDataPhyloLikelihood& asd) :
      AbstractAlignedPhyloLikelihood(asd),
      nbStates_(asd.nbStates_),
      nData_(asd.nData_)
    {
    }
    
    virtual ~AbstractSingleDataPhyloLikelihood() {}
    
    AbstractSingleDataPhyloLikelihood* clone() const = 0;
    
    AbstractSingleDataPhyloLikelihood& operator=(const AbstractSingleDataPhyloLikelihood& asd)
    {
      AbstractAlignedPhyloLikelihood::operator=(asd);
      nbStates_=asd.nbStates_;
      
      nData_=asd.nData_;
      
      return *this;
    }
    
    virtual void setData(const SiteContainer& sites, size_t nData = 0)
    {
      setNumberOfSites(sites.getNumberOfSites());
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
    
    size_t getNumberOfStates() const { return nbStates_; }
    
  };
      
} //end of namespace bpp.

#endif  //_SINGLEDATAPHYLOLIKELIHOOD_H_

