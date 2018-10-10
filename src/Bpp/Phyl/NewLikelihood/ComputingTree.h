//
// File: ComputingTree.h
// Created by: Laurent Guéguen
// Created on: Sat Dec 30 12:48 2006
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

#include "../Tree/PhyloTree.h"
#include "../Tree/PhyloBranchParam.h"
#include "../Tree/AwareNode.h"

#include "SpeciationComputingNode.h"
#include "ParametrizablePhyloTree.h"

namespace bpp
{
  class SubstitutionProcessCollection;

/**
 * @brief Tree Organization of Computing Nodes
 *
 * Stores computation tools for all nodes, for all classes.
 *
 * This object has the parameters of the Tree and the rates
 * distribution, since it manages the SpeciationComputingNodes.
 *
 */

  class ComputingTree :
    public AbstractParametrizable
  {
    typedef AssociationTreeGlobalGraphObserver<ComputingNode, PhyloBranchParam>  CompTree;

    typedef AssociationTreeGlobalGraphObserver<SpeciationComputingNode, PhyloBranchParam>  SpecCompTree;

  private:
    /*
     * A pointer towards a Parametrizable Phylo Tree
     *
     */

    const ParametrizablePhyloTree* parTree_;

    /*
     * A pointer towards a Discrete Distribution
     *
     */

    const DiscreteDistribution* pDist_;

    /*
     * a vector of trees of computing nodes
     *
     */
    
    std::vector<std::shared_ptr<CompTree> > vTree_;

    /*
     * boolean to say if the ComputingTree can be used for
     * computation, ie if every SpeciationComputingNode has a Model.
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
     * @brief construction of an empty ComputingTree with no Distribution
     *
     * @param tree The tree.
     *
     */
     
    ComputingTree(const ParametrizablePhyloTree& ptree);

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
     
    ComputingTree(const ParametrizablePhyloTree& ptree, const DiscreteDistribution& dist);

    ComputingTree(const ComputingTree& tree);
 
    ComputingTree& operator=(const ComputingTree& tree);

    ComputingTree* clone() const { return new ComputingTree(*this);}
      
    ~ComputingTree();

    /*
     * @brief construction of a complete ComputingTree.
     *
     * @param pSubMod a  pointer of TransitionModel.
     * @param vBr a vector of attribution of the model on the
     * branches of the tree that need a model.
     *
     */
     
    void addModel(const TransitionModel* pSubMod, std::vector<unsigned int> vBr);

    /*
     * @brief construction of an homogeneous ComputingTree.
     *
     * @param pSubMod a  pointer of TransitionModel.
     *
     */
     
    void addModel(const TransitionModel* pSubMod);
    
    size_t getNumberOfClasses() const { return vTree_.size();}

  private:
    
    void clearAllModels_();

  public:

    /*
     * @brief Checks if every SpeciationComputingNode has a Model
     *
     */
    
    void checkModelOnEachNode();

    /*
     * @brief operator to get numbered TreeTemplate<SpeciationComputingNode>*
     *
     */
    
    std::shared_ptr<CompTree> operator[](size_t ntree) { return vTree_[ntree];}

    const std::shared_ptr<CompTree> operator[](size_t ntree) const { return vTree_[ntree];}

    /*
     *@brief update Distribution parameters and says to the
     * SpeciationComputingNodes to be ready to update if the Branch lengths are
     * changed.
     *
     */
    
    void fireParameterChanged(const ParameterList& pl);
        
    /*
     * @brief Says to specific nodes to be ready for update
     *
     * If flag = true (default), node has to be updated (false otherwise).
     */

    void update(std::vector<unsigned int>& vId, bool flag = true);

    void update(unsigned int id, bool flag = true);

    /*
     * @brief Says to all nodes to be ready for update
     *
     */
    
    void updateAll();

    /*
     * @brief Returns the list of the updated Nodes ids, excluding the root.
     *
     */
    
    Vuint toBeUpdatedNodes() const;

  private:

    /*
     * @brief Initializes the Computing tree in homogeneous constructors.
     *
     */
    
    void init_();
    
  };
  
} //end of namespace bpp.

#endif //_COMPUTINGTREE_H_

