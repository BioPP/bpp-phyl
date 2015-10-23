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
 * @brief Interface LikelihoodTree data structure.
 *
 * Stores all the inner computations:
 * - conditional likelihoods for each node,
 * - correspondance between sites in the dataset and array indices.
 *
 * The structure is initiated according to a tree topology, and 
 * data can be retrieved through node ids.
 *
 * @see LikelihoodNode
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

    /**
     * @brief The probabilities of the classes.
     */
    
    Vdouble vProbClass_;

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
    
    void resetLikelihoods(int nodeId, size_t nbSites, size_t nbStates, unsigned char DX)
    {
      for (size_t i=0; i<nbClasses_; i++)
        getNodeData(nodeId, i).resetLikelihoods(nbSites, nbStates, DX);
    }

    /*
     * @brief recursively reset Likelihoods on all the Inner nodes.
     *
     * This is done given the useLog
     */
    
    void resetDownwardLikelihoods(size_t nbSites, size_t nbStates, unsigned char DX)
    {
      for (size_t i=0; i<nbClasses_; i++)
        getRootData(i).resetDownwardLikelihoods(nbSites, nbStates, DX);
    }

    /*
     * @brief recursively set patterns on all the nodes.
     *
     */
    
    void setPatterns(const std::map<int, std::map<int, std::vector<size_t> > >& patterns)
    {
      for (size_t i=0; i<nbClasses_; i++)
        getRootData(i).setPatterns(patterns);
    }
    
public:

    /*
     * @brief the relations between real position and shrunked data
     * positions.
     *
     */
     
    std::vector<size_t>& getRootArrayPositions() { return rootPatternLinks_; }
      
    size_t getRootArrayPosition(size_t currentPosition) const
    {
      return rootPatternLinks_[currentPosition];
    }
    
    const std::vector<size_t>& getRootArrayPositions() const { return rootPatternLinks_; }

    const SiteContainer* getShrunkData() const {
      return shrunkData_.get();
    }

    size_t getNumberOfDistinctSites() const { return nbDistinctSites_; }

    size_t getNumberOfSites() const { return nbSites_; }

    /*
     * @brief the weights of the positions of the shrunked data
     *
     */
     
    unsigned int getWeight(size_t pos) const
    {
      return rootWeights_[pos];
    }

    const std::vector<unsigned int>& getWeights() const
    { 
      return rootWeights_;
    }

    /*
     * @brief the alphabet.
     *
     */
    
    const Alphabet* getAlphabet() const { return alphabet_; }

    size_t getNumberOfStates() const { return nbStates_; }

    size_t getNumberOfClasses() const { return nbClasses_; }

    /*
     * @brief the Node Data
     *
     */

    virtual AbstractLikelihoodNode& getNodeData(int nodeId, size_t nClass) = 0;

    virtual const AbstractLikelihoodNode& getNodeData(int nodeId, size_t nClass) const = 0;

    virtual AbstractLikelihoodNode& getRootData(size_t nClass) = 0;

    virtual const AbstractLikelihoodNode& getRootData(size_t nClass) const = 0;

    /*
     * @brief the DXLikehoodArrays
     *
     */
    
    VVdouble& getLikelihoodArray(int nodeId, size_t nClass, unsigned char DX)
    {
      return getNodeData(nodeId, nClass).getLikelihoodArray(DX);
    }

    /**
     * @brief returns if given node at given class uses log in arrays.
     *
     */
    
    bool usesLog(int nodeId, size_t nClass) const
    {
      return getNodeData(nodeId, nClass).usesLog();
    }
    
    /**
     * @brief returns if root at given class uses log in arrays.
     *
     */
    
    bool usesLogAtRoot(size_t nClass) const
    {
      return getRootData(nClass).usesLog();
    }

    /**
     * @brief sets using log in all likelihood arrays.
     *
     */
    
    void setAllUseLog(bool useLog)
    {
      for (size_t i=0; i<nbClasses_;i++)
        getRootData(i).setUseLogDownward(useLog);
    }
    

    /**
     * @brief Compute the posterior probabilities for each state and
     * each class of each distinct site.
     *
     * @param nodeId The id of the node at which probabilities must be
     * computed.
     *
     * @return A 3-dimensional array, with, for each site, joint
     * probabilities of all classes and all states (sum on each
     * site=1).
     *
     */
    
    VVVdouble getPosteriorProbabilitiesForEachStateForEachClass(int nodeId);

    /**
     * @brief Compute the posterior probabilities for each state for a
     * given node.
     *
     * This method calls the
     * getPosteriorProbabilitiesForEachStateForEachClass function and
     * average the probabilities over all sites and classes,
     * resulting in a one-dimensionnal frequency array, with one
     * frequency per model state.
     *
     * @param nodeId The id of the node at which probabilities must be
     * computed.
     * @return vector of double with state frequencies for the given
     * node.
     */
    
    Vdouble getPosteriorStateFrequencies(int nodeId);

  };
  
} //end of namespace bpp.

#endif //_ABSTRACT_LIKELIHOOD_TREE_H_

