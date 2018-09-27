//
// File: ProbabilisticSubstitutionMappingForASite.h
// Created by: Laurent 
// Created on: mercredi 26 septembre 2018, à 16h 06
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

#ifndef _PROBABILISTICSUBSTITUTIONMAPPING_FOR_A_SITE_H_
#define _PROBABILISTICSUBSTITUTIONMAPPING_FOR_A_SITE_H_

#include "SubstitutionMapping.h"
#include "PhyloBranchMappingForASite.h"
#include "../Tree/PhyloTree.h"

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

  class ProbabilisticSubstitutionMappingForASite:
    public AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchMappingForASite>
  {
  private:
    size_t numberOfTypes_;

  public:
    typedef AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchMappingForASite> mapTree;

    typedef mapTree::EdgeIterator EdgeIterator;
    typedef mapTree::NodeIterator NodeIterator;
    
  public:
    /**
     * @brief Build a new ProbabilisticSubstitutionMappingForASite object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     * @param sc A pointer toward the substitution count object that has been used for the mapping, if any.
     * This object allows to get the substitution types description, if there are several. If set to 0, then
     * the mapping will be considered as having only one type of substitution mapped.
     * @param numberOfSites The number of sites to map.
     */

    ProbabilisticSubstitutionMappingForASite(const PhyloTree& tree, size_t nbTypes) :
      mapTree(tree), numberOfTypes_(nbTypes)
    {
      std::unique_ptr<mapTree::EdgeIterator> nIT=allEdgesIterator();
      for (;!nIT->end(); nIT->next())
        (**nIT)->setNumberOfTypes(nbTypes);
    }

    /**
     * @brief Build a new ProbabilisticSubstitutionMappingForASite object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     */
    
    ProbabilisticSubstitutionMappingForASite(const PhyloTree& tree) :
      mapTree(tree), numberOfTypes_(0) 
    {}
    

    ProbabilisticSubstitutionMappingForASite* clone() const { return new ProbabilisticSubstitutionMappingForASite(*this); }

    ProbabilisticSubstitutionMappingForASite(const ProbabilisticSubstitutionMappingForASite& psm):
      mapTree(psm), numberOfTypes_(psm.numberOfTypes_)
    {}

    ProbabilisticSubstitutionMappingForASite& operator=(const ProbabilisticSubstitutionMappingForASite& psm)
    {
      mapTree::operator=(psm);
      numberOfTypes_ = psm.numberOfTypes_;
      
      return *this;
    }

    virtual ~ProbabilisticSubstitutionMappingForASite() {}

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

    /**
     * @brief Retrieve the counts
     *
     */
    
    double getCount(unsigned int branchId, size_t type) const
    {
      return getEdge(branchId)->getTypeCount(type);
    }

    size_t getNumberOfTypes() const
    {
      return numberOfTypes_;
    }
    
    void setNumberOfTypes(size_t numberOfTypes)
    {
      numberOfTypes_=numberOfTypes;
      
      std::unique_ptr<mapTree::EdgeIterator> nIT=allEdgesIterator();
      for (;!nIT->end(); nIT->next())
        (**nIT)->setNumberOfTypes(numberOfTypes);
    }
    
    /**
     * @brief Direct access to substitution numbers
     */

    double operator()(unsigned int branchId, size_t type) const
    {
      return (*getEdge(branchId))(type);
    }

  };

} //end of namespace bpp.

#endif //_PROBABILISTICSUBSTITUTIONMAPPING_FOR_A_SITE_H_

