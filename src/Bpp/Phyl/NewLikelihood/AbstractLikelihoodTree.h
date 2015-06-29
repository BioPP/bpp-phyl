//
// File: AbstractLikelihoodTree.h
// Created by: Laurent Guéguen
// Created on: mardi 23 juin 2015, à 09h 19
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

#ifndef _ABSTRACT_LIKELIHOOD_TREE_H_
#define _ABSTRACT_LIKELIHOOD_TREE_H_

#include "SubstitutionProcess.h"

#include "../Tree/TreeTemplate.h"
#include "AbstractLikelihoodNode.h"
#include "LikelihoodTree.h"
#include "ComputingTree.h"


namespace bpp
{

/**
 * @brief Tree Organization of Computing Nodes
 *
 * Stores computation tools for all nodes, for all classes.
 *
 * This object has the parameters of the Tree and the rates
 * distribution, since it manages the ComputingNodes.
 *
 */

  class AbstractLikelihoodTree :
    public virtual LikelihoodTree
  {
  protected:
    /**
     * @brief Links between sites and patterns.
     * 
     * The size of this vector is equal to the number of sites in the container,
     * each element corresponds to a site in the container and points to the
     * corresponding column in the likelihood array of the root node.
     * If the container contains no repeated site, there will be a strict
     * equivalence between each site and the likelihood array of the root node.
     * However, if this is not the case, some pointers may point toward the same
     * element in the likelihood array.
     */
    std::vector<size_t> rootPatternLinks_;

    /**
     * @brief The frequency of each site.
     */
    std::vector<unsigned int> rootWeights_;

    const Alphabet* alphabet_;

    std::auto_ptr<SiteContainer> shrunkData_;
    size_t nbSites_; 
    size_t nbStates_;
    size_t nbClasses_;
    size_t nbDistinctSites_; 

  public:
    /*
     * @brief construction of an empty AbstractLikelihoodTree from a tree and a
     * discretedistribution belonging to a collection
     *
     * @param subPro the SubstitutionProcess 
     *
     */
     
    AbstractLikelihoodTree(const SubstitutionProcess& subPro);
  
    AbstractLikelihoodTree(const AbstractLikelihoodTree& tree);
 
    AbstractLikelihoodTree& operator=(const AbstractLikelihoodTree& tree);

    ~AbstractLikelihoodTree();

  public:

    /*
     * @brief reset the Likelihood arrays:
     *     resize to nbSites X nbStates
     *     set to 1 and 0
     *
     */
    
    void resetLikelihoods(int nodeId, size_t nbSites, size_t nbStates)
    {
      for (size_t i=0; i<nbClasses_; i++)
        getNodeData(nodeId, i).resetLikelihoods(nbSites, nbStates);
    }

    /*
     * @brief reset the Likelihood arrays:
     *     set to 1 and 0
     */

    void resetLikelihoods(int nodeId, unsigned char DX)
    {
      for (size_t i=0; i<nbClasses_; i++)
        getNodeData(nodeId, i).resetLikelihoods(DX);
    }

  public:
    std::vector<size_t>& getRootArrayPositions() { return rootPatternLinks_; }
      
    const std::vector<size_t>& getRootArrayPositions() const { return rootPatternLinks_; }

    virtual AbstractLikelihoodNode& getNodeData(int nodeId, size_t nClass) = 0;

    virtual const AbstractLikelihoodNode& getNodeData(int nodeId, size_t nClass) const = 0;
    
    size_t getRootArrayPosition(size_t site) const
    {
      return rootPatternLinks_[site];
    }

    unsigned int getWeight(size_t pos) const
    {
      return rootWeights_[pos];
    }

    const std::vector<unsigned int>& getWeights() const
    { 
      return rootWeights_;
    }

    const Alphabet* getAlphabet() const { return alphabet_; }

    size_t getNumberOfDistinctSites() const { return nbDistinctSites_; }
    size_t getNumberOfSites() const { return nbSites_; }
    size_t getNumberOfStates() const { return nbStates_; }
    size_t getNumberOfClasses() const { return nbClasses_; }

    const SiteContainer* getShrunkData() const {
      return shrunkData_.get();
    }

    void computeUpwardPartialLikelihoods(const ComputingTree& lTree, int nodeId, unsigned char DX, Vint* vbr=NULL)
    {
      for (size_t c = 0; c < nbClasses_; ++c)
        getNodeData(nodeId, c).computeUpwardPartialLikelihoods(*lTree[c]->getNode(nodeId), DX, vbr);
    }
    
    void computeUpwardPartialLikelihoods(const ComputingTree& lTree, int nodeId, const std::vector<const std::vector<size_t>* >& vPatterns, unsigned char DX, Vint* vbr=NULL)
    {
      for (size_t c = 0; c < nbClasses_; ++c)
        getNodeData(nodeId, c).computeUpwardPartialLikelihoods(*lTree[c]->getNode(nodeId), vPatterns, DX, vbr);
    }

  };
  
} //end of namespace bpp.

#endif //_ABSTRACT_LIKELIHOOD_TREE_H_

