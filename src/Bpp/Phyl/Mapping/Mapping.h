// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_MAPPING_H
#define BPP_PHYL_MAPPING_MAPPING_H

#include <Bpp/Clonable.h>
#include <Bpp/Graph/AssociationTreeGraphImplObserver.h>

#include "../Tree/PhyloTree.h"
#include "PhyloBranchMapping.h"

// From the STL:
#include <vector>
#include <memory>

namespace bpp
{
/**
 * @brief General interface for storing mapping data.
 *
 */
class MappingInterface :
  public virtual Clonable
{
public:
  MappingInterface() {}
  virtual ~MappingInterface() {}

  MappingInterface* clone() const override = 0;

public:
  /**
   * @return The number of sites mapped.
   */
  virtual size_t getNumberOfSites() const = 0;

  /**
   * @return The PhyloBranchMapping to a given index
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

class AbstractMapping :
  public virtual MappingInterface
{
protected:
  std::vector<int> sitesPositions_;
  size_t nbSites_;

public:
  AbstractMapping() : sitesPositions_(), nbSites_(0)
  {}

  AbstractMapping(size_t nb) : sitesPositions_(), nbSites_(nb)
  {}

  AbstractMapping(const AbstractMapping& absm) = default;

  AbstractMapping& operator=(const AbstractMapping& absm) = default;

  virtual ~AbstractMapping() {}

public:

  int getSitePosition(size_t index) const
  {
    return sitesPositions_.size() == 0 ? (int)(index + 1) : sitesPositions_[index];
  }

  void setSitePosition(size_t index, int position)
  {
    if (sitesPositions_.size() == 0)
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
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_MAPPING_H
