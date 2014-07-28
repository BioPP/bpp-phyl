//
// File: ComputingTree.h
// Created by: Laurent Guéguen
// Created on: Sat Dec 30 12:48 2006
// From file AbstractTreeLikelihood.h
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

#ifndef _COMPUTINGTREE_H_
#define _COMPUTINGTREE_H_

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

#include "../Tree/Node.h"
#include "../Tree/TreeTemplate.h"

#include "ComputingNode.h"
#include "ParametrizableTree.h"

namespace bpp
{

  class SubstitutionProcessCollection;

/**
 * @brief Tree Organization of Computing Nodes
 *
 * Stores computation tools for all nodes, for all classes.
 *
 * This object has the parameters of the Tree and the rates
 * distribution, since it manages the ComputingNodes.
 *
 */

  class ComputingTree :
    public AbstractParametrizable
  {
  private:
    /*
     * A pointer towards a Parametrizable Tree
     *
     */

    const ParametrizableTree* parTree_;

    /*
     * A pointer towards a Discrete Distribution
     *
     */

    const DiscreteDistribution* pDist_;

    /*
     * a vector of trees of computing nodes
     *
     */
    
    std::vector<TreeTemplate<ComputingNode>* > vTree_;

    /*
     * boolean to say if the ComputingTree can be used for
     * computation, ie if every ComputingNode has a Model.
     *
     */

    bool isReadyToCompute_;
    
  public:
    /*
     * @brief construction of an empty ComputingTree from a tree and a
     * discretedistribution belonging to a collection
     *
     * @param pSubProColl the SubstitutionProcessCollection
     * @param nTree The number of the tree.
     * @param nDist the number of the rate distribution.
     *
     */
     
    ComputingTree(const SubstitutionProcessCollection* pSubProColl, size_t nTree, size_t nDist);
  
    /*
     * @brief construction of an empty ComputingTree with Constant Distribution
     *
     * @param tree The tree.
     *
     */
     
    ComputingTree(const ParametrizableTree& ptree);

    /*
     * @brief construction of an empty ComputingTree.
     *
     * Warning : the links towards the models do not follow the copy,
     *    so the computing tree has to be filled with models before
     *    computation performs.
     *
     * @param tree The tree.
     * @param dist the rate distribution.
     *
     */
     
    ComputingTree(const ParametrizableTree& ptree, const DiscreteDistribution& dist);

    ComputingTree(const ComputingTree& tree);
 
    ComputingTree& operator=(const ComputingTree& tree);

    ComputingTree* clone() const { return new ComputingTree(*this);}
      
    ~ComputingTree();

    /*
     * @brief construction of a complete ComputingTree.
     *
     * @param pSubMod a  pointer of SubstitutionModel.
     * @param vBr a vector of attribution of the model on the
     * branches of the tree.
     *
     */
     
    void addModel(const SubstitutionModel* pSubMod, std::vector<int> vBr);

    /*
     * @brief construction of an homogeneous ComputingTree.
     *
     * @param pSubMod a  pointer of SubstitutionModel.
     *
     */
     
    void addModel(const SubstitutionModel* pSubMod);
    
    size_t getNumberOfClasses() const { return vTree_.size();}

  private:
    
    void clearAllModels_();

  public:

    /*
     * @brief Checks if every ComputingNode has a Model
     *
     */
    
    void checkModelOnEachNode();

    /*
     * @brief operator to get numbered TreeTemplate<ComputingNode>*
     *
     */
    
    TreeTemplate<ComputingNode>* operator[](size_t ntree) { return vTree_[ntree];}

    const TreeTemplate<ComputingNode>* operator[](size_t ntree) const { return vTree_[ntree];}

    /*
     *@brief update Distribution parameters and says to the
     * ComputingNodes to be ready to update if the Branch lengths are
     * changed.
     *
     */
    
    void fireParameterChanged(const ParameterList& pl);
        
    /*
     * @brief Says to specific nodes to be ready for update
     *
     */

    void update(std::vector<int>& vId);

    /*
     * @brief Says to all nodes to be ready for update
     *
     */
    
    void updateAll();

    /**
     * @brief Methods for computing of the partial likelihoods.
     *
     */

    /**
     *@brief multiplies a partial likelihood using the partial
     * likelihoods of the sons.
     *
     * @param likelihoods_node a pointer to the partial likelihood
     * of this node [in, out].
     * @param vLikelihoods_sons a vector of the partial likelihoods of
     * the sons.
     * @param nodeId the Id of the node which partial likelihood is
     * computed.
     * @param vPatterns a vector of the corresponding positions
     * from this node to the sons.
     **/
    
    
    void multiplyPartialLikelihoods(VVVdouble* likelihoods_node, const std::vector<const VVVdouble*>& vLikelihoods_sons, int nodeId, const std::vector<const std::vector<size_t>* >& vPatterns) const
    {
      size_t nbSons=vLikelihoods_sons.size();
      for (size_t c = 0; c < vTree_.size(); ++c)
      {
        VVdouble* likelihoods_node_c = &(*likelihoods_node)[c];

        std::vector<const VVdouble*> vLikelihoods_sons_c;
        for (size_t i=0;i<nbSons;i++)
          vLikelihoods_sons_c.push_back(&(*vLikelihoods_sons[i])[c]);

        vTree_[c]->getNode(nodeId)->multiplyPartialLikelihoods(likelihoods_node_c, vLikelihoods_sons_c, vPatterns);
      }
    }

    
    /**
     *@brief multiplies a partial likelihood using the partial
     * likelihoods of the sons.
     *
     * @param likelihoods_node a pointer to the partial likelihood
     * of this node [in, out].
     * @param vLikelihoods_sons a vector of the partial likelihoods of
     * the sons.
     * @param nodeId the Id of the node which partial likelihood is
     * computed.
     **/

    
    void multiplyPartialLikelihoods(VVVdouble* likelihoods_node, const std::vector<const VVVdouble*>& vLikelihoods_sons, int nodeId) const
    {
      size_t nbSons=vLikelihoods_sons.size();
      for (size_t c = 0; c < vTree_.size(); ++c)
      {
        VVdouble* likelihoods_node_c = &(*likelihoods_node)[c];

        std::vector<const VVdouble*> vLikelihoods_sons_c;
        for (size_t i=0;i<nbSons;i++)
          vLikelihoods_sons_c.push_back((vLikelihoods_sons[i]!=0)?&(*vLikelihoods_sons[i])[c]:0);

        vTree_[c]->getNode(nodeId)->multiplyPartialLikelihoods(likelihoods_node_c, vLikelihoods_sons_c);
      }
    }

    /**
     *@brief multiplies the partial likelihood using the partial
     * likelihoods of a son.
     *
     * @param likelihoods_node a pointer to the partial likelihood
     * of this node [in, out].
     * @param likelihoods_son  the partial likelihood of
     * the used son.
     * @param nodeId the Id of the node which partial likelihood is
     * computed.
     * @param sonId the used Computing son Id.
     **/

    void multiplyPartialLikelihoods(VVVdouble* likelihoods_node, const VVVdouble* likelihoods_son, int nodeId, int sonId) const
    {
      for (size_t c = 0; c < vTree_.size(); ++c)
      {
        VVdouble* likelihoods_node_c = &(*likelihoods_node)[c];

        const VVdouble* likelihoods_sons_c = &(*likelihoods_son)[c];

        vTree_[c]->getNode(nodeId)->multiplyPartialLikelihoods(likelihoods_node_c, likelihoods_sons_c, vTree_[c]->getNode(sonId));
      }
    }


    /**
     *@brief multiplies the derivate of the partial likelihood using
     * the partial likelihoods of the sons of a given node.
     *
     * @param brId is the Id of the branch on which the derivate is
     * computed.
     * @param dXLikelihoods_node a pointer to the derivate of the
     * partial likelihood of this node [in, out].
     * @param vLikelihoods_sons a vector of the partial likelihoods of
     * the sons.
     * @param nodeId the Id of the node which the partial likelihood
     * is computed.
     * @param vPatterns a vector of the corresponding positions
     * from this node to the sons.
     * @param DX tells which matrix should be used as used for
     * transition factors, either D0 for transition probabilities, D1
     * for their first derivate, D2 for their second.
     **/

    void multiplyPartialDXLikelihoods(int brId, VVVdouble* dXLikelihoods_node, const std::vector<const VVVdouble*>& vLikelihoods_sons, int nodeId, const std::vector<const std::vector<size_t>* >& vPatterns, unsigned char DX) const
    {
      size_t nbSons=vLikelihoods_sons.size();
      
      for (size_t c = 0; c < vTree_.size(); ++c)
      {
        VVdouble* dXLikelihoods_node_c = &(*dXLikelihoods_node)[c];
        std::vector<const VVdouble*> vLikelihoods_sons_c;
        for (size_t i=0;i<nbSons;i++)
          vLikelihoods_sons_c.push_back(&(*vLikelihoods_sons[i])[c]);
        
        vTree_[c]->getNode(nodeId)->multiplyPartialDXLikelihoods(brId, dXLikelihoods_node_c, vLikelihoods_sons_c, vPatterns, DX);
      }
    }

    /**
     *@brief computes the derivate of the likelihood by the length of
     * a branch, using the conditional likelihoods of the father of
     * this node (in case of double recursive calculation).
     *
     * @param dXlikelihoods_node a pointer to the derivate of the
     * likelihood of this node [in, out].
     * @param likelihoods_father_node a pointer to the conditional
     * father-node likelihood.
     * @param vLikelihoods_sons a vector of the likelihoods of the
     * sons.
     * @param nodeId the Id of the node which the derivate is
     * computed.
     * @param DX tells which matrix should be used as used for
     * transition factors, either D0 for transition probabilities, D1
     * for their first derivate, D2 for their second.
     **/

    void computeDXLogLikelihoods(VVdouble* dXLogLikelihoods_node, const VVVdouble* likelihoods_father_node, const VVVdouble* likelihoods_father, int nodeId, unsigned char DX) const
    {
      for (size_t c = 0; c < vTree_.size(); ++c)
      {
        Vdouble* dXLogLikelihoods_node_c = &(*dXLogLikelihoods_node)[c];
        const VVdouble* likelihoods_father_node_c = &(*likelihoods_father_node)[c];
        const VVdouble* likelihoods_father_c = &(*likelihoods_father)[c];
        
        vTree_[c]->getNode(nodeId)->computeDXLogLikelihoods(dXLogLikelihoods_node_c, likelihoods_father_node_c, likelihoods_father_c, DX);
      }
    }
    
  };
  
    
} //end of namespace bpp.

#endif //_COMPUTINGTREE_H_

