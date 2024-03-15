// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _LEGACY_PROBABILISTIC_SUBSTITUTION_MAPPING_H_
#define _LEGACY_PROBABILISTIC_SUBSTITUTION_MAPPING_H_

#include "SubstitutionMapping.h"
#include "../../Mapping/SubstitutionCount.h" //We use the latest version of the class here
#include "../../Tree/TreeExceptions.h"

#include <Bpp/Text/TextTools.h>

// From the STL:
#include <vector>

namespace bpp
{

/**
 * @brief Legacy data storage class for probabilistic substitution mappings.
 *
 * A 'probabilistic' mapping contains an single value for each branch and each site.
 * This number can be an average number of substitutions, optionally waited, or a probability of observing a certain number of substitutions.
 * Probabilistic was coined there by opposition to the'stochastic' mapping, where a path (number of susbstitutions + there position along the branch)
 * is available for each branch and site. The probabilistic mapping can however be extended to contain a matrix will all types of substitutions, instead of their total number.
 */
class LegacyProbabilisticSubstitutionMapping:
  public LegacyAbstractSubstitutionMapping
{
  private:
    std::shared_ptr<const SubstitutionCountInterface> substitutionCount_;
    
    /**
     * @brief Substitution numbers storage.
     *
     * Numbers are stored by sites.
     */
    std::vector<std::vector<std::vector<double>>> mapping_;
  
  public:
    
    /**
     * @brief Build a new LegacyProbabilisticSubstitutionMapping object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     * @param sc A pointer toward the substitution count object that has been used for the mapping, if any.
     * This object allows to get the substitution types description, if there are several. If set to 0, then
     * the mapping will be considered as having only one type of substitution mapped.
     * @param numberOfSites The number of sites to map.
     */
    LegacyProbabilisticSubstitutionMapping(const Tree& tree, std::shared_ptr<const SubstitutionCountInterface> sc, size_t numberOfSites) :
      LegacyAbstractSubstitutionMapping(tree),
      substitutionCount_(sc),
      mapping_(0)
    {
      setNumberOfSites(numberOfSites);
    }

    /**
     * @brief Build a new ProbabilisticSubstitutionMapping object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     */
    LegacyProbabilisticSubstitutionMapping(const Tree& tree) :
      LegacyAbstractSubstitutionMapping(tree),
      substitutionCount_(nullptr),
      mapping_(0)
    {}
    

    LegacyProbabilisticSubstitutionMapping* clone() const override
    {
      return new LegacyProbabilisticSubstitutionMapping(*this);
    }

    LegacyProbabilisticSubstitutionMapping(const LegacyProbabilisticSubstitutionMapping& psm) = default;

    LegacyProbabilisticSubstitutionMapping& operator=(const LegacyProbabilisticSubstitutionMapping& psm) = default;

    virtual ~LegacyProbabilisticSubstitutionMapping() {}

  public:

    size_t getNumberOfSubstitutionTypes() const override
    {
      if (!substitutionCount_) return 1;
      return substitutionCount_->getNumberOfSubstitutionTypes();
    }
     
    double getNumberOfSubstitutions(int nodeId, size_t siteIndex, size_t type) const
    {
      return mapping_[siteIndex][getNodeIndex(nodeId)][type];
    }
    
    virtual std::vector<double> getNumberOfSubstitutions(int nodeId, size_t siteIndex) const
    {
      return mapping_[siteIndex][getNodeIndex(nodeId)];
    }
    
    /**
     * @brief (Re)-set the phylogenetic tree associated to this mapping.
     *
     * @param tree The new tree.
     */
    virtual void setTree(const Tree& tree);

    virtual void setNumberOfSites(size_t numberOfSites) override;
    
    /**
     * @brief Direct access to substitution numbers.
     *
     * @warning No index checking is performed, use with care!
     */
    virtual double& operator()(size_t nodeIndex, size_t siteIndex, size_t type) override
    {
      return mapping_[siteIndex][nodeIndex][type];
    }

    /**
     * @brief Direct access to substitution numbers.
     *
     * @warning No index checking is performed, use with care!
     */
    virtual const double& operator()(size_t nodeIndex, size_t siteIndex, size_t type) const override
    {
      return mapping_[siteIndex][nodeIndex][type];
    }
     
    /**
     * @brief Direct access to substitution numbers.
     *
     * @warning No index checking is performed, use with care!
     */
    std::vector< std::vector<double> >& operator[](size_t siteIndex)
    {
      return mapping_[siteIndex];
    }

    /**
     * @brief Direct access to substitution numbers.
     *
     * @warning No index checking is performed, use with care!
     */
    const std::vector< std::vector<double> >& operator[](size_t siteIndex) const
    {
      return mapping_[siteIndex];
    }
};

} //end of namespace bpp.

#endif //_LEGACY_PROBABILISTIC_SUBSTITUTION_MAPPING_H_

