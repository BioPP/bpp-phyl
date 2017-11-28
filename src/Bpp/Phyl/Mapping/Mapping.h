//
// File: Mapping.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 09:51 2005
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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

#ifndef MAPPING_H_
#define MAPPING_H_

#include <Bpp/Clonable.h>

#include "../Tree/PhyloTree.h"
#include "PhyloBranchMapping.h"

#include <Bpp/Graph/AssociationTreeGraphImplObserver.h>

//From the STL:
#include <vector>
#include <memory>

namespace bpp
{

  /**
   * @brief General interface for storing mapping data.
   *
   */
  class Mapping:
    public virtual Clonable
  {

  public:
    Mapping() {}
    virtual ~Mapping() {}

    Mapping* clone() const = 0;

  public:
    
    /**
     * @return The number of sites mapped.
     */
    virtual size_t getNumberOfSites() const = 0;

    /**
     * @return The PhyloBranchMapping to a given index
     *
     */
    
    virtual const PhyloBranch& getBranch(unsigned int branchId) const = 0;

    virtual PhyloBranch& getBranch(unsigned int branchId) = 0;

    // virtual unsigned int getBranchIndex(const std::shared_ptr<PhyloBranch> branch) const = 0;
    
    /**
     * @return The number of branches mapped.
     */
    
    virtual size_t getNumberOfBranches() const = 0;
    
    /**
     * @param index The site index.
     * @return The site position corresponding to the index.
     */
    
    virtual int getSitePosition(size_t index) const = 0;
    
    /**
     * @brief Set the position of a given site.
     *
     * @warning No index checking is performed, use with care!
     * @param index The site index.
     * @param position The position of the site.
     */
    
    virtual void setSitePosition(size_t index, int position) = 0;
  };




  /**
   * @brief Partial implementation of the mapping interface.
   * 
   *  Here, there is no information about site compression.
   *
   */
  
  class AbstractMapping:
    virtual public Mapping
  {
  protected:
    std::vector<int> sitesPositions_;
    size_t nbSites_;

  public:
    AbstractMapping(): sitesPositions_(), nbSites_(0)
    {
    }

    AbstractMapping(size_t nb): sitesPositions_(), nbSites_(nb)
    {      
    }

    AbstractMapping(const AbstractMapping& absm):
      sitesPositions_(absm.sitesPositions_),
      nbSites_(absm.nbSites_)
    {
    }

    AbstractMapping& operator=(const AbstractMapping& absm)
    {
      sitesPositions_ = absm.sitesPositions_;
      nbSites_        = absm.nbSites_;
      
      return *this;
    }

    AbstractMapping* clone() const = 0;
  
    virtual ~AbstractMapping() {}
  
  public:
  
    int getSitePosition(size_t index) const
    {
      return (sitesPositions_.size()==0?(int)(index+1):sitesPositions_[index]);
    }
    
    void setSitePosition(size_t index, int position)
    {
      if (sitesPositions_.size()==0)
        sitesPositions_.resize(nbSites_);
      
      sitesPositions_[index] = position;
    }
		
    size_t getNumberOfSites() const { return nbSites_; }

    virtual void setNumberOfSites(size_t numberOfSites)
    {
      nbSites_ = numberOfSites;
      sitesPositions_.clear();
    }

  };

} //end of namespace bpp.

#endif //_MAPPING_H_

