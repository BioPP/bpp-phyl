// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
    AbstractRewardMapping(), mapTree(tree), rootPatternLinks_(), usePatterns_(false), numberOfDistinctSites_(numberOfSites)
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
    AbstractRewardMapping(), mapTree(tree), rootPatternLinks_(rootpatterns), usePatterns_(true), numberOfDistinctSites_(nbDistinctSites)
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
    AbstractRewardMapping(), mapTree(tree), rootPatternLinks_(), usePatterns_(false), numberOfDistinctSites_(0)
  {}


  ProbabilisticRewardMapping* clone() const override { return new ProbabilisticRewardMapping(*this); }

  ProbabilisticRewardMapping(const ProbabilisticRewardMapping& prm) :
    AbstractRewardMapping(prm), mapTree(prm), rootPatternLinks_(prm.rootPatternLinks_), usePatterns_(prm.usePatterns_), numberOfDistinctSites_(prm.numberOfDistinctSites_)
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
  const PhyloBranch& getBranch(unsigned int branchIndex) const override
  {
    return *getEdge(branchIndex);
  }

  PhyloBranch& getBranch(unsigned int branchIndex) override
  {
    return *getEdge(branchIndex);
  }

  size_t getNumberOfBranches() const override
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

  virtual void setNumberOfSites(size_t numberOfSites) override;

  /**
   * @brief Direct access to rewards, with COMPRESSED
   * site positions (ie site indexes)
   *
   * @warning No index checking is performed, use with care!
   */
  virtual double operator()(uint branchId, size_t siteIndex) const override
  {
    return (*getEdge(branchId))(siteIndex);
  }

  /**
   * @brief Direct access to rewards, with COMPRESSED
   * site positions (ie site indexes).
   *
   * @warning No index checking is performed, use with care!
   */
  virtual double& operator()(uint branchId, size_t siteIndex) override
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
