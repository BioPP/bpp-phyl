//
// File: PhyloLikelihoodContainer.h
// Created by: Laurent Guéguen
// Created on: mercredi 7 octobre 2015, à 22h 34
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

#ifndef _PHYLOLIKELIHOOD_CONTAINER_H_
#define _PHYLOLIKELIHOOD_CONTAINER_H_

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainer.h>

#include "PhyloLikelihood.h"
#include "SingleDataPhyloLikelihood.h"

namespace bpp
{
  /**
   * @brief The PhyloLikelihoodContainer, assigns numbers to
   * Phylolikelihoods.
   *
   * It owns the PhyloLikelihoods
   *
   */
  
  class PhyloLikelihoodContainer:
    virtual public Clonable
  {
  protected:
    std::map<size_t, PhyloLikelihood*>  mPhylo_;

  public:
    PhyloLikelihoodContainer():
      mPhylo_()
    {}
    
    PhyloLikelihoodContainer(const PhyloLikelihoodContainer& lik):
    mPhylo_()
    {
      for (std::map<size_t, PhyloLikelihood*>::const_iterator it=lik.mPhylo_.begin(); it!=lik.mPhylo_.end(); it++)
        mPhylo_[it->first]=it->second->clone();
    }

    PhyloLikelihoodContainer& operator=(const PhyloLikelihoodContainer& lik)
    {
      for (std::map<size_t, PhyloLikelihood*>::const_iterator it=mPhylo_.begin(); it!=mPhylo_.end(); it++)
        delete it->second;
      
      mPhylo_.clear();
      
      for (std::map<size_t, PhyloLikelihood*>::const_iterator it=lik.mPhylo_.begin(); it!=lik.mPhylo_.end(); it++)
        mPhylo_[it->first]=it->second->clone();
      
      return *this;
    }

    PhyloLikelihoodContainer* clone() const
    {
      return new PhyloLikelihoodContainer(*this);
    }
    
    /**
     * @brief Abstract class destructor
     *
     */
    
    virtual ~PhyloLikelihoodContainer()
    {
      for (std::map<size_t, PhyloLikelihood*>::iterator it=mPhylo_.begin(); it!=mPhylo_.end(); it++)
        delete it->second;
      mPhylo_.clear();
    }

  public:

    virtual void addPhyloLikelihood(size_t pos, PhyloLikelihood* Ap)
    {
      if (mPhylo_.find(pos)!=mPhylo_.end())
        throw Exception("PhyloLikelihoodContainer::addPhylolikelihood: map number already used : " + TextTools::toString(pos));
      mPhylo_[pos]=Ap;
    }

    bool hasPhyloLikelihood(size_t pos) const
    {
      return (mPhylo_.find(pos)!=mPhylo_.end());
    }

    virtual const PhyloLikelihood* operator[](size_t pos) const
    {
      std::map<size_t, PhyloLikelihood*>::const_iterator it=mPhylo_.find(pos);
      return (it!=mPhylo_.end()?it->second:0);
    }

    virtual PhyloLikelihood* operator[](size_t pos)
    {
      std::map<size_t, PhyloLikelihood*>::iterator it=mPhylo_.find(pos);
      return (it!=mPhylo_.end()?it->second:0);
    }
    
    size_t getSize() const
    {
      return mPhylo_.size();
    }

    std::vector<size_t> getNumbersOfPhyloLikelihoods() const
    {
      std::vector<size_t> vnum;
      
      for (std::map<size_t, PhyloLikelihood*>::const_iterator it=mPhylo_.begin(); it!=mPhylo_.end(); it++)
        vnum.push_back(it->first);

      return vnum;
    }
    
    /**
     * @brief Set the dataset for which the likelihood must be
     * evaluated, iff the pointed PhyloLikelihood is a SingleDataPhyloLikelihood
     *
     * @param nPhyl The number of the Likelihood.
     * @param sites The data set to use.
     */

    void setData(const SiteContainer& sites, size_t nPhyl)
    {
      std::map<size_t, PhyloLikelihood*>::iterator it=mPhylo_.find(nPhyl);
      if (it!=mPhylo_.end())
      {
        SingleDataPhyloLikelihood* sdp=dynamic_cast<SingleDataPhyloLikelihood*>(it->second);
        if (sdp)
          sdp->setData(sites);
      }
    }
    
    
    /**
     * @brief Get the dataset for which the likelihood must be evaluated.
     *
     * @return A pointer toward the site container where the sequences are stored.
     */

    const SiteContainer* getData(size_t nPhyl) const
    {
      std::map<size_t, PhyloLikelihood*>::const_iterator it=mPhylo_.find(nPhyl);
      if (mPhylo_.find(nPhyl)!=mPhylo_.end())
      {
        const SingleDataPhyloLikelihood* sdp=dynamic_cast<const SingleDataPhyloLikelihood*>(it->second);
        if (sdp)
          sdp->getData();
      }
      return 0;
    }
      

  };
  
  
} //end of namespace bpp.

#endif  //_PHYLOLIKELIHOOD_CONTAINER_H_

