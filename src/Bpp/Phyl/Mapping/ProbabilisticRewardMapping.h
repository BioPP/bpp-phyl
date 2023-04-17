//
// File: ProbabilisticRewardMapping.h
// Authors:
//   Laurent Guéguen
// Created: lundi 20 novembre 2017, ÃÂ  16h 55
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef BPP_PHYL_MAPPING_PROBABILISTICREWARDMAPPING_H
#define BPP_PHYL_MAPPING_PROBABILISTICREWARDMAPPING_H

#include <Bpp/Text/TextTools.h>

#include "../Likelihood/DataFlow/DataFlowCWise.h"
#include "../Tree/PhyloTreeExceptions.h"
#include "PhyloBranchReward.h"
#include "Reward.h"
#include "RewardMapping.h"

// From the STL:
#include <vector>

namespace bpp
{
/**
 * @brief Data storage class for probabilistic rewards mappings.
 *
 * A 'probabilistic' mapping contains an single value for each branch and each site.
 * This number is an average reward.
 * Probabilistic was coined there by opposition to the'stochastic'
 * mapping, where a path (sequence of rewards along the branch) is
 * available for each branch and site.
 */
class ProbabilisticRewardMapping :
  public AbstractRewardMapping,
  public AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchReward>
{
public:
  typedef AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchReward> mapTree;

private:
  /**
   * @brief Links between sites and patterns.
   *
   * The size of this vector is equal to the number of sites in the container,
   * each element corresponds to a site in the container and points to the
   * corresponding column in the count array of the root node.
   * If the container contains no repeated site, there will be a strict
   * equivalence between each site and the likelihood array of the root node.
   * However, if this is not the case, some pointers may point toward the same
   * element in the likelihood array.
   */
  PatternType rootPatternLinks_;

  bool usePatterns_;

  size_t numberOfDistinctSites_;

public:
  typedef mapTree::EdgeIterator EdgeIterator;
  typedef mapTree::NodeIterator NodeIterator;

public:
  /**
   * @brief Build a new ProbabilisticRewardMapping object.
   *
   * @param tree The tree object to use. It will be cloned for internal use.
   * @param numberOfSites The number of sites to map.
   */
  ProbabilisticRewardMapping(const PhyloTree& tree, size_t numberOfSites) :
    AbstractMapping(numberOfSites), AbstractRewardMapping(), mapTree(tree), rootPatternLinks_(), usePatterns_(false), numberOfDistinctSites_(numberOfSites)
  {
    std::unique_ptr<mapTree::EdgeIterator> nIT = allEdgesIterator();
    for ( ; !nIT->end(); nIT->next())
    {
      (**nIT)->setNumberOfSites(numberOfDistinctSites_);
    }
  }

  /**
   * @brief the same with rootPatternLinks
   */
  ProbabilisticRewardMapping(const PhyloTree& tree, const PatternType& rootpatterns, size_t nbDistinctSites) :
    AbstractMapping(size_t(rootpatterns.size())), AbstractRewardMapping(), mapTree(tree), rootPatternLinks_(rootpatterns), usePatterns_(true), numberOfDistinctSites_(nbDistinctSites)
  {
    std::unique_ptr<mapTree::EdgeIterator> nIT = allEdgesIterator();
    for ( ; !nIT->end(); nIT->next())
    {
      (**nIT)->setNumberOfSites(numberOfDistinctSites_);
    }
  }

  /**
   * @brief Build a new ProbabilisticRewardMapping object.
   *
   * @param tree The tree object to use. It will be cloned for internal use.
   */
  ProbabilisticRewardMapping(const PhyloTree& tree) :
    AbstractMapping(), AbstractRewardMapping(), mapTree(tree), rootPatternLinks_(), usePatterns_(false), numberOfDistinctSites_(0)
  {}


  ProbabilisticRewardMapping* clone() const { return new ProbabilisticRewardMapping(*this); }

  ProbabilisticRewardMapping(const ProbabilisticRewardMapping& prm) :
    AbstractMapping(prm), AbstractRewardMapping(prm), mapTree(prm), rootPatternLinks_(prm.rootPatternLinks_), usePatterns_(prm.usePatterns_), numberOfDistinctSites_(prm.numberOfDistinctSites_)
  {}

  ProbabilisticRewardMapping& operator=(const ProbabilisticRewardMapping& prm)
  {
    AbstractMapping::operator=(prm);
    mapTree::operator=(prm);

    rootPatternLinks_ = prm.rootPatternLinks_;
    usePatterns_ = prm.usePatterns_;
    numberOfDistinctSites_ = prm.numberOfDistinctSites_;
    return *this;
  }

  virtual ~ProbabilisticRewardMapping() {}

public:
  /*
   * @brief From Mapping interface
   *
   * @{
   */
  const PhyloBranch& getBranch(unsigned int branchIndex) const
  {
    return *getEdge(branchIndex);
  }

  PhyloBranch& getBranch(unsigned int branchIndex)
  {
    return *getEdge(branchIndex);
  }

  size_t getNumberOfBranches() const
  {
    return getNumberOfEdges();
  }

  size_t getNumberOfDistinctSites() const
  {
    return numberOfDistinctSites_;
  }

  /*
   * @}
   */

  /**
   * @brief Retrieve the rewards, with compressed site positions.
   *
   */
  virtual double getReward(int branchId, size_t site) const
  {
    return getEdge((uint)branchId)->getSiteReward(getSiteIndex(site));
  }

  virtual void setNumberOfSites(size_t numberOfSites);

  /**
   * @brief Direct access to rewards, with COMPRESSED
   * site positions (ie site indexes)
   *
   * @warning No index checking is performed, use with care!
   */
  virtual double operator()(uint branchId, size_t siteIndex) const
  {
    return (*getEdge(branchId))(siteIndex);
  }

  /**
   * @brief Direct access to rewards, with COMPRESSED
   * site positions (ie site indexes).
   *
   * @warning No index checking is performed, use with care!
   */
  virtual double& operator()(uint branchId, size_t siteIndex)
  {
    return (*getEdge(branchId))(siteIndex);
  }

  /**
   * @brief Does it use site patterns?
   */
  bool usePatterns() const
  {
    return usePatterns_;
  }

  /**
   * @brief returns the vector of site patterns
   */
  const PatternType& getPatterns() const
  {
    return rootPatternLinks_;
  }

  const size_t getSiteIndex(size_t site) const
  {
    return rootPatternLinks_[Eigen::Index(site)];
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_PROBABILISTICREWARDMAPPING_H
