// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_PROBABILISTICSUBSTITUTIONMAPPING_H
#define BPP_PHYL_MAPPING_PROBABILISTICSUBSTITUTIONMAPPING_H

#include <Bpp/Text/TextTools.h>

#include "../Likelihood/DataFlow/DataFlowCWise.h"
#include "SubstitutionMapping.h"

// From the STL:
#include <vector>

namespace bpp
{
/**
 * @brief Data storage class for probabilistic substitution mappings.
 *
 * A 'probabilistic' mapping contains an single value for each branch
 * and each site. This number can be an average number of
 * substitutions, optionally waited, or a probability of observing a
 * certain number of substitutions. Probabilistic was coined there by
 * opposition to the'stochastic' mapping, where a path (number of
 * susbstitutions + there position along the branch) is available for
 * each branch and site. The probabilistic mapping can however be
 * extended to contain a matrix will all types of substitutions,
 * instead of their total number.
 */
class ProbabilisticSubstitutionMapping :
  public AbstractSubstitutionMapping,
  public AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchMapping>
{
public:
  typedef AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchMapping> mapTree;

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
   * @brief Build a new ProbabilisticSubstitutionMapping object.
   *
   * @param tree The PhyloTree object to use.
   * @param nbTypes the number of types
   * @param numberOfSites The number of sites to map.
   */
  ProbabilisticSubstitutionMapping(const PhyloTree& tree, size_t nbTypes, size_t numberOfSites) :
    AbstractMapping(numberOfSites), AbstractSubstitutionMapping(), mapTree(tree), rootPatternLinks_(), usePatterns_(false), numberOfDistinctSites_(numberOfSites)
  {
    setNumberOfSubstitutionTypes(nbTypes);

    std::unique_ptr<EdgeIterator> nIT = allEdgesIterator();
    for ( ; !nIT->end(); nIT->next())
    {
      (**nIT)->setNumberOfSitesAndTypes(numberOfDistinctSites_, nbTypes);
    }
  }

  /**
   * @brief the same with rootPatternLinks
   */
  ProbabilisticSubstitutionMapping(const PhyloTree& tree, size_t nbTypes, const PatternType& rootpatterns, size_t nbDistinctSites) :
    AbstractMapping(size_t(rootpatterns.size())), AbstractSubstitutionMapping(), mapTree(tree), rootPatternLinks_(rootpatterns), usePatterns_(true), numberOfDistinctSites_(nbDistinctSites)
  {
    setNumberOfSubstitutionTypes(nbTypes);
    std::unique_ptr<EdgeIterator> nIT = allEdgesIterator();
    for ( ; !nIT->end(); nIT->next())
    {
      (**nIT)->setNumberOfSitesAndTypes(numberOfDistinctSites_, nbTypes);
    }
  }

  /**
   * @brief Build a new ProbabilisticSubstitutionMapping object.
   *
   * @param tree The tree object to use. It will be cloned for internal use.
   */
  ProbabilisticSubstitutionMapping(const PhyloTree& tree) :
    AbstractMapping(), AbstractSubstitutionMapping(), mapTree(tree), rootPatternLinks_(), usePatterns_(false), numberOfDistinctSites_(0)
  {}

  ProbabilisticSubstitutionMapping* clone() const override
  {
    return new ProbabilisticSubstitutionMapping(*this); 
  }

  ProbabilisticSubstitutionMapping(const ProbabilisticSubstitutionMapping& psm) :
    AbstractMapping(psm), AbstractSubstitutionMapping(psm), mapTree(psm), rootPatternLinks_(psm.rootPatternLinks_), usePatterns_(psm.usePatterns_), numberOfDistinctSites_(psm.numberOfDistinctSites_)
  {}

  ProbabilisticSubstitutionMapping& operator=(const ProbabilisticSubstitutionMapping& psm)
  {
    AbstractSubstitutionMapping::operator=(psm);
    mapTree::operator=(psm);

    rootPatternLinks_ = psm.rootPatternLinks_;
    usePatterns_ = psm.usePatterns_;
    numberOfDistinctSites_ = psm.numberOfDistinctSites_;

    return *this;
  }

  virtual ~ProbabilisticSubstitutionMapping() {}

public:
  const PhyloBranchMapping& getBranch(unsigned int branchIndex) const override
  {
    return *getEdge(branchIndex);
  }

  PhyloBranchMapping& getBranch(unsigned int branchIndex) override
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

  /**
   * @brief Retrieve the counts, with REAL site positions.
   */
  double getCount(unsigned int branchId, size_t site, size_t type) const
  {
    return getEdge(branchId)->getSiteTypeCount(getSiteIndex(site), type);
  }

  const std::vector<double>& getCounts(unsigned int branchId, size_t site) const
  {
    return getEdge(branchId)->getSiteCount(getSiteIndex(site));
  }

  void setNumberOfSitesAndTypes(size_t numberOfSites, size_t numberOfTypes);

  void setNumberOfSites(size_t numberOfSites) override;

  void setNumberOfSubstitutionTypes(size_t numberOfTypes) override;

  /**
   * @brief Direct access to substitution numbers, with COMPRESSED
   * site positions (ie site indexes)
   *
   * @warning No index checking is performed, use with care!
   */
  double& operator()(unsigned int branchId, size_t siteIndex, size_t type) override
  {
    return (*getEdge(branchId))(siteIndex, type);
  }

  /**
   * @brief Fill a VVdouble with the counts at a given site.
   *    The 1st coordinate of this VVdouble correspond to edges,
   *    ordered through their ids, & the second coordinate
   *    corresponds to the type numbers
   *
   * Branches are read in the same order as returned by
   * getAllEdgesIndexes() function.
   *
   */
  void fillMappingVectorForSite(size_t siteIndex, VVdouble& counts) const;

  /**
   * @brief Direct access to substitution numbers, with COMPRESSED
   * site positions (ie site indexes)
   *
   * @warning No index checking is performed, use with care!
   */
  virtual const double& operator()(unsigned int branchId, size_t siteIndex, size_t type) const override
  {
    return (*getEdge(branchId))(siteIndex, type);
  }

  /**
   * @brief Does it use site patterns?
   *
   */
  bool usePatterns() const
  {
    return usePatterns_;
  }

  /**
   * @brief returns the vector of site patterns
   *
   */
  const PatternType getPatterns() const
  {
    return rootPatternLinks_;
  }

  const size_t getSiteIndex(size_t site) const
  {
    return usePatterns_ ? rootPatternLinks_(Eigen::Index(site)) : site;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_PROBABILISTICSUBSTITUTIONMAPPING_H
