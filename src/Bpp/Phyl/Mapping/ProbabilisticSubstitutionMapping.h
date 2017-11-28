//
// File: ProbabilisticSubstitutionMapping.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 10:47 2006
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

#ifndef _PROBABILISTICSUBSTITUTIONMAPPING_H_
#define _PROBABILISTICSUBSTITUTIONMAPPING_H_

#include "SubstitutionMapping.h"

#include <Bpp/Text/TextTools.h>

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
 *
 */

  class ProbabilisticSubstitutionMapping:
    public AbstractSubstitutionMapping,
    public AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchMapping>
  {
  public:

    typedef AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchMapping> mapTree;

  private:
    /*
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
    
    std::vector<size_t> rootPatternLinks_;
    
    bool usePatterns_;

    size_t numberOfDistinctSites_;

  public:

    typedef mapTree::EdgeIterator EdgeIterator;
    typedef mapTree::NodeIterator NodeIterator;
    
  public:
    /**
     * @brief Build a new ProbabilisticSubstitutionMapping object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     * @param sc A pointer toward the substitution count object that has been used for the mapping, if any.
     * This object allows to get the substitution types description, if there are several. If set to 0, then
     * the mapping will be considered as having only one type of substitution mapped.
     * @param numberOfSites The number of sites to map.
     */

    ProbabilisticSubstitutionMapping(const PhyloTree& tree, size_t nbTypes, size_t numberOfSites) :
      AbstractMapping(numberOfSites), AbstractSubstitutionMapping(), mapTree(tree), rootPatternLinks_(), usePatterns_(false), numberOfDistinctSites_(numberOfSites) 
    {
      setNumberOfSubstitutionTypes(nbTypes);
      
      std::unique_ptr<mapTree::EdgeIterator> nIT=allEdgesIterator();
      for (;!nIT->end(); nIT->next())
        (**nIT)->setNumberOfSitesAndTypes(numberOfDistinctSites_, nbTypes);
    }

    /**
     * @brief the same with rootPatternLinks
     *
     */
    
    ProbabilisticSubstitutionMapping(const PhyloTree& tree, size_t nbTypes, const std::vector<size_t>& rootpatterns, size_t nbDistinctSites) :
      AbstractMapping(rootpatterns.size()), AbstractSubstitutionMapping(), mapTree(tree), rootPatternLinks_(rootpatterns), usePatterns_(true), numberOfDistinctSites_(nbDistinctSites)
    {
      setNumberOfSubstitutionTypes(nbTypes);
      std::unique_ptr<mapTree::EdgeIterator> nIT=allEdgesIterator();
      for (;!nIT->end(); nIT->next())
        (**nIT)->setNumberOfSitesAndTypes(numberOfDistinctSites_, nbTypes);
    }

    /**
     * @brief Build a new ProbabilisticSubstitutionMapping object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     */
    
    ProbabilisticSubstitutionMapping(const PhyloTree& tree) :
      AbstractMapping(), AbstractSubstitutionMapping(), mapTree(tree), rootPatternLinks_(), usePatterns_(false), numberOfDistinctSites_(0) 
    {}
    

    ProbabilisticSubstitutionMapping* clone() const { return new ProbabilisticSubstitutionMapping(*this); }

    ProbabilisticSubstitutionMapping(const ProbabilisticSubstitutionMapping& psm):
      AbstractMapping(psm), AbstractSubstitutionMapping(psm), mapTree(psm), rootPatternLinks_(psm.rootPatternLinks_), usePatterns_(psm.usePatterns_), numberOfDistinctSites_(psm.numberOfDistinctSites_) 
    {}

    ProbabilisticSubstitutionMapping& operator=(const ProbabilisticSubstitutionMapping& psm)
    {
      AbstractSubstitutionMapping::operator=(psm);
      mapTree::operator=(psm);

      rootPatternLinks_=psm.rootPatternLinks_;
      usePatterns_=psm.usePatterns_;
      numberOfDistinctSites_=psm.numberOfDistinctSites_;
      
      return *this;
    }

    virtual ~ProbabilisticSubstitutionMapping() {}

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

    /*
     * @}
     */

    size_t getNumberOfDistinctSites() const
    {
      return numberOfDistinctSites_;
    }

    /**
     * @brief Retrieve the counts, with REAL site positions.
     *
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

    void setNumberOfSites(size_t numberOfSites);

    void setNumberOfSubstitutionTypes(size_t numberOfTypes);

    /**
     * @brief Direct access to substitution numbers, with COMPRESSED
     * site positions (ie site indexes)
     *
     * @warning No index checking is performed, use with care!
     */

    double& operator()(unsigned int branchId, size_t siteIndex, size_t type)
    {
      return (*getEdge(branchId))(siteIndex, type);
    }

    /**
     * @brief Direct access to substitution numbers, with COMPRESSED
     * site positions (ie site indexes)
     *
     * @warning No index checking is performed, use with care!
     */
    
    virtual const double& operator()(unsigned int branchId, size_t siteIndex, size_t type) const
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

    const std::vector<size_t>& getPatterns() const
    {
      return rootPatternLinks_;
    }

    const size_t getSiteIndex(size_t site) const
    {
      return (usePatterns_?rootPatternLinks_[site]:site);
    }
  };

} //end of namespace bpp.

#endif //_PROBABILISTICSUBSTITUTIONMAPPING_H_

