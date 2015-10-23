//
// File: RecursiveLikelihoodTree.h
// Created by: Julien Dutheil, Laurent Guéguen
// Created on: Sat Dec 30 14:20 2006
// From file HomogeneousLikelihoodTree.h
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

#ifndef _RECURSIVE_LIKELIHOOD_TREE_H_
#define _RECURSIVE_LIKELIHOOD_TREE_H_

#include "AbstractLikelihoodTree.h"
#include "SubstitutionProcess.h"
#include "../SitePatterns.h"
#include "LikelihoodNode.h"
#include "RecursiveLikelihoodNode.h"

#include <Bpp/Numeric/VectorTools.h>

// From the STL:
#include <map>
using namespace std;

namespace bpp
{
/**
 * @brief Likelihood data structure suporting simple recursion.
 */
  class RecursiveLikelihoodTree :
    public AbstractLikelihoodTree
  {
  private:
    /*
     * a vector of trees of computing nodes
     *
     */
    
    std::vector<TreeTemplate<RecursiveLikelihoodNode>* > vTree_;


/**
 * @brief This map defines the pattern network.
 *
 * Let n1 be the id of a node in the tree, and n11 and n12 the ids of its sons.
 * Providing the likelihood array is known for nodes n11 and n12,
 * the likelihood array for node n1 and site <i>i</i> (_likelihood[n1][i]) must be computed  
 * using arrays patternLinks_[n1][n11][i] and patternLinks_[n1][n12][i].
 * This network is intialized once for all in the constructor of this class.
 *
 * The double map contains the position of the site to use (second dimension)
 * of the likelihoods array.
 */

    // TO DO: convert to a tree structure
    
    mutable std::map<int, std::map<int, std::vector<size_t> > > patternLinks_;

    bool usePatterns_;

    /*
     * @brief check in the Above Likelihoods have been initialized
     *
     */
    
    bool initializedAboveLikelihoods_;

  public:
    
    RecursiveLikelihoodTree(const SubstitutionProcess& process, bool usePatterns);

    RecursiveLikelihoodTree(const RecursiveLikelihoodTree& data);
    
    RecursiveLikelihoodTree& operator=(const RecursiveLikelihoodTree & data);

    virtual ~RecursiveLikelihoodTree();

    RecursiveLikelihoodTree* clone() const { return new RecursiveLikelihoodTree(*this); }

  public:

    /*
     * @brief operator to get numbered TreeTemplate<ComputingNode>
     *
     */
    
    TreeTemplate<RecursiveLikelihoodNode>& operator[](size_t ntree) { return *vTree_[ntree];}

    const TreeTemplate<RecursiveLikelihoodNode>& operator[](size_t ntree) const { return *vTree_[ntree];}

    /*
     * @brief the Node Data
     *
     */
    
    AbstractLikelihoodNode& getNodeData(int nodeId, size_t nClass)
    { 
      return *(*this)[nClass].getNode(nodeId);
    }

    const AbstractLikelihoodNode& getNodeData(int nodeId, size_t nClass) const
    { 
      return *(*this)[nClass].getNode(nodeId);
    }

    AbstractLikelihoodNode& getRootData(size_t nClass)
    { 
      return *(*this)[nClass].getRootNode();
    }

    const AbstractLikelihoodNode& getRootData(size_t nClass) const
    { 
      return *(*this)[nClass].getRootNode();
    }

    /*
     * @brief the recursive pattern relation between positions.
     *
     */
    
    bool usePatterns()
    {
      return usePatterns_;
    }

    /*
     * @brief resets the likelihood arrays.
     *
     */
    
    void resetBelowLikelihoods(int nodeId, size_t nbSites, size_t nbStates, unsigned char DX) const
    {
      for (size_t c = 0; c < vTree_.size(); ++c)
        vTree_[c]->getNode(nodeId)->resetBelowLikelihoods(nbSites, nbStates, DX);
    }

    void resetAboveLikelihoods(int nodeId, size_t nbSites, size_t nbStates) const
    {
      for (size_t c = 0; c < vTree_.size(); ++c)
        vTree_[c]->getNode(nodeId)->resetAboveLikelihoods(nbSites, nbStates);
    }

    void setAboveLikelihoods(int nodeId, const Vdouble& freq) {
      for (size_t c = 0; c < vTree_.size(); ++c)
        vTree_[c]->getNode(nodeId)->setAboveLikelihoods(freq);
    }

    /*
     * @brief set if log-likelihood should be computed at node.
     *
     */
    
    void setUseLog(int nodeId, bool useLog) const
    {
      for (size_t c = 0; c < vTree_.size(); ++c)
        vTree_[c]->getNode(nodeId)->setUseLog(useLog);
    }

    /*
     * @brief reset the Above Likelihood arrays all the inner nodes:
     * resize to nbSites_ X nbStates_ set to 1
     *
     */
    
    void resetInnerAboveLikelihoods();


    bool isAboveLikelihoodsInitialized() const
    {
      return initializedAboveLikelihoods_;
    }
    
    /*
     * @brief initialize the likelihoods from the data & the process.
     *
     */
    
    void initLikelihoods(const SiteContainer& sites, const SubstitutionProcess& process) throw (Exception);

    /*
     * @brief compute full likelihoods at a given node
     *
     */
    
    void computeLikelihoodsAtNode(const ComputingTree& lTree, int nodeId)
    {
      for (size_t c = 0; c < vTree_.size(); ++c)
        vTree_[c]->getNode(nodeId)->computeLikelihoods(*(lTree[c]->getNode(nodeId)), ComputingNode::D0);
    }

    /*
     * @brief compute full DXlikelihoods. In case of derivation, brId
     * is a pointer on a vector of ids of derivated branches.
     *
     */

    void computeLikelihoods(const ComputingTree& lTree, unsigned char DX, Vint* brId = NULL)
    {
      int rId=lTree[0]->getRootId();
      
      for (size_t c = 0; c < vTree_.size(); ++c)
        vTree_[c]->getNode(rId)->computeLikelihoods(*(lTree[c]->getNode(rId)), DX, brId);
    }

  protected:
    
    /**
     * @brief This method initializes the leaves according to a sequence file.
     * likelihood is set to 1 for the state corresponding to the sequence site,
     * otherwise it is set to 0.
     *
     * All likelihood arrays at each nodes are initialized according to alphabet
     * size and sequences length, and filled with 1.
     *
     * NB: This method is recursive.
     *
     * @param node      The node defining the subtree to analyse.
     * @param sequences The data to be used for initialization.
     * @param process   The substitution process to use.
     */
    
    virtual void initLikelihoodsWithoutPatterns_(const Node* node, const SiteContainer& sequences, const SubstitutionProcess& process) throw (Exception);
    
    /**
     * @brief This method initializes the leaves according to a sequence file.
     *
     * likelihood is set to 1 for the state corresponding to the sequence site,
     * otherwise it is set to 0.
     *
     * All likelihood arrays at each nodes are initialized according to alphabet
     * size and sequences length, and filled with 1.
     *
     * NB: This method is recursive.
     *
     * @param node      The node defining the subtree to analyse.
     * @param sequences The data to be used for initialization.
     * @param process   The substitution process to use.
     * @return The shrunk sub-dataset + indices for the subtree defined by <i>node</i>.
     */
    
    virtual SitePatterns* initLikelihoodsWithPatterns_(const Node* node, const SiteContainer& sequences, const SubstitutionProcess& process) throw (Exception);
    
  };
  
} //end of namespace bpp.

#endif //_RECURSIVE_LIKELIHOOD_TREE_H_

