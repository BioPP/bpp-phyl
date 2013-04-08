//
// File: ModelIterator.h
// Created by: Julien Dutheil
// Created on: Tue Jul 17 09:50 2012
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

#ifndef _MODELITERATOR_H_
#define _MODELITERATOR_H_

#include "../Model/SubstitutionModel.h"
#include "ParametrizableTree.h"

namespace bpp
{
/**
 * @brief An iterator over a set of branches, specified by their node ids.
 */
class BranchIterator
{
public:
  virtual ~BranchIterator() {}

public:
  /**
   * @return The id of the next node in the set.
   */
  virtual int next() throw (Exception) = 0;
  /**
   * @return True if there is at least another node in the set.
   */
  virtual bool hasNext() const = 0;
};

/**
 * @brief An iterator over a set of sites, speicfied by their position.
 *
 * In most cases, the position will reflect the index of an inner array used for likelihood storage.
 */
class SiteIterator
{
public:
  virtual ~SiteIterator() {}

public:
  /**
   * @return The position of the next site in the set.
   */
  virtual unsigned int next() throw (Exception) = 0;
  /**
   * @return True is there is at least another site in the set.
   */
  virtual bool hasNext() const = 0;
};

/**
 * @brief A pair of SubstitutionModel / SiteIterator.
 */
class ConstBranchModelDescription
{
public:
  virtual ~ConstBranchModelDescription() {}

public:
  virtual SiteIterator* getNewSiteIterator() const = 0;
};

/**
 * @brief Iterates through all models used for all sites on a given branch.
 */
class ConstBranchModelIterator
{
public:
  virtual ~ConstBranchModelIterator() {}

public:
  virtual ConstBranchModelDescription* next() throw (Exception) = 0;
  virtual bool hasNext() const = 0;
};

/**
 * @brief A pair of SubstitutionModel / BranchIterator.
 */
class ConstSiteModelDescription
{
public:
  virtual ~ConstSiteModelDescription() {}

public:
  virtual BranchIterator* getNewBranchIterator() const = 0;
};

/**
 * @brief Iterates through all models used for all branches on a given site.
 */
class ConstSiteModelIterator
{
public:
  virtual ~ConstSiteModelIterator() {}

public:
  virtual ConstSiteModelDescription* next() throw (Exception) = 0;
  virtual bool hasNext() const = 0;
};

/**
 * @brief A very simple branch iterator.
 *
 * The constructor takes a vector of nodes id to iterate over.
 */
class SimpleBranchIterator :
  public BranchIterator
{
private:
  std::vector<int> nodesId_;
  unsigned int index_;

public:
  SimpleBranchIterator(const std::vector<int>& nodesId) :
    nodesId_(nodesId),
    index_(0) {}

public:
  int next() throw (Exception)
  {
    if (!hasNext())
      throw Exception("AbstractTreeLikelihood::SimpleBranchIterator::next(). No more branch in the set.");
    return nodesId_[index_++];
  }

  bool hasNext() const { return index_ < nodesId_.size(); }
};

/**
 * @brief A very simple site iterator.
 *
 * This iterator loops over a continuous range of sites.
 * The constructor takes as input the number of sites to iterate over,
 * and optionally an offset argument, specifying the index of the first site.
 */
class SimpleSiteIterator :
  public SiteIterator
{
private:
  unsigned int maxIndex_;
  unsigned int index_;
  unsigned int offset_;

public:
  SimpleSiteIterator(unsigned int nbSites, unsigned int offset = 0) :
    maxIndex_(nbSites),
    index_(0),
    offset_(offset) {}

public:
  unsigned int next() throw (Exception)
  {
    if (!hasNext())
      throw Exception("SimpleSiteIterator::next(). No more site in the set.");
    return offset_ + index_++;
  }

  bool hasNext() const { return index_ < maxIndex_; }
};

/**
 * @name Branch iterator for models without site partition.
 *
 * @{
 */
class ConstNoPartitionBranchModelDescription :
  public ConstBranchModelDescription
{
private:
  unsigned int nbSites_;

public:
  ConstNoPartitionBranchModelDescription(const SubstitutionModel* model, unsigned int nbSites) :
    nbSites_(nbSites) {}

  ConstNoPartitionBranchModelDescription(const ConstNoPartitionBranchModelDescription& bmd) :
    nbSites_(bmd.nbSites_)
  {}

  ConstNoPartitionBranchModelDescription& operator=(const ConstNoPartitionBranchModelDescription& bmd)
  {
    nbSites_ = bmd.nbSites_;
    return *this;
  }

public:
  SiteIterator* getNewSiteIterator() const { return new SimpleSiteIterator(nbSites_); }
};

class ConstNoPartitionBranchModelIterator :
  public ConstBranchModelIterator
{
private:
  ConstNoPartitionBranchModelDescription branchModelDescription_;
  unsigned int index_;

public:
  ConstNoPartitionBranchModelIterator(const SubstitutionModel* model, unsigned int nbSites) :
    branchModelDescription_(model, nbSites),
    index_(0) {}

public:
  ConstNoPartitionBranchModelDescription* next() throw (Exception)
  {
    if (!hasNext())
      throw Exception("ConstNoPartitionBranchModelDescription::next(). No more branch in the set.");
    index_++;
    return &branchModelDescription_;
  }

  bool hasNext() const { return index_ == 0; }
};

class ConstNoPartitionSiteModelDescription :
  public ConstSiteModelDescription
{
private:
  std::vector<int> nodesId_;

public:
  ConstNoPartitionSiteModelDescription(const SubstitutionModel* model, const std::vector<int> nodesId) :
    nodesId_(nodesId) {}

  ConstNoPartitionSiteModelDescription(const ConstNoPartitionSiteModelDescription& smd) :
    nodesId_(smd.nodesId_)
  {}

  ConstNoPartitionSiteModelDescription& operator=(const ConstNoPartitionSiteModelDescription& smd)
  {
    nodesId_ = smd.nodesId_;
    return *this;
  }

public:
  BranchIterator* getNewBranchIterator() const { return new SimpleBranchIterator(nodesId_); }
};

/** @} */

class ConstHomogeneousSiteModelIterator :
  public ConstSiteModelIterator
{
  private:
    ConstNoPartitionSiteModelDescription siteModelDescription_;
    unsigned int index_;

  public:
    ConstHomogeneousSiteModelIterator(const ParametrizableTree& tree, const SubstitutionModel* model) :
     siteModelDescription_(model, tree.getBranchesId()), index_(0) {}

  public:
    ConstSiteModelDescription* next() throw (Exception)
    {
      if (!hasNext())
        throw Exception("ConstHomogeneousSiteModelIterator::next(). No more site in the set.");
      index_++;
      return &siteModelDescription_;
    }

    bool hasNext() const { return index_ == 0; }
};


} // end of namespace bpp.

#endif // _MODELITERATOR_H_
