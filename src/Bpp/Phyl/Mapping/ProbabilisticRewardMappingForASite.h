//
// File: ProbabilisticRewardMapping.h
// Created by: Laurent Guéguen
// Created on: lundi 20 novembre 2017, à 16h 55
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

#ifndef _PROBABILISTICREWARDMAPPING_FOR_A_SITE_H_
#define _PROBABILISTICREWARDMAPPING_FOR_A_SITE_H_

#include "PhyloBranchRewardForASite.h"
#include "../Tree/PhyloTree.h"

#include "Reward.h"
#include "../Tree/TreeExceptions.h"

#include <Bpp/Text/TextTools.h>

// From the STL:
#include <vector>

namespace bpp
{

/**
 * @brief Data storage class for probabilistic rewards mappings.
 *
 * A 'probabilistic' mapping contains an single value for each branch.
 * This number is an average reward.
 * Probabilistic was coined there by opposition to the'stochastic'
 * mapping, where a path (sequence of rewards along the branch) is
 * available for each branch and site.
 */
  class ProbabilisticRewardMappingForASite:
    public AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchRewardForASite>
  { 
  public:
    typedef AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchRewardForASite> mapTree;

    typedef mapTree::EdgeIterator EdgeIterator;
    typedef mapTree::NodeIterator NodeIterator;
    
  public:
    
    /**
     * @brief Build a new ProbabilisticRewardMappingForASite object.
     *
     * @param tree The tree object to use. It will be cloned for internal use.
     * @param reward A pointer toward the Reward object that has been used for the mapping, if any.
     * @param numberOfSites The number of sites to map.
     */
    ProbabilisticRewardMappingForASite(const PhyloTree& tree) :
      mapTree(tree)
    {
    }
    
    ProbabilisticRewardMappingForASite* clone() const { return new ProbabilisticRewardMappingForASite(*this); }

    ProbabilisticRewardMappingForASite(const ProbabilisticRewardMappingForASite& prm):
      mapTree(prm)
    {}

    ProbabilisticRewardMappingForASite& operator=(const ProbabilisticRewardMappingForASite& prm)
    {
      mapTree::operator=(prm);

      return *this;
    }

    virtual ~ProbabilisticRewardMappingForASite() {}

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
     * @brief Retrieve the rewards, with compressed site positions.
     *
     */
    
    virtual double getReward(int branchId) const
    {
      return getEdge(branchId)->getReward();
    }
  };

} //end of namespace bpp.

#endif //_PROBABILISTICREWARDMAPPING_FOR_A_SITE_H_

