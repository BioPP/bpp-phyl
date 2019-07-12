//
// File: ForwardLikelihoodTree.cpp
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

#ifndef _FORWARD_LIKELIHOOD_TREE_H_
#define _FORWARD_LIKELIHOOD_TREE_H_

#include "Bpp/NewPhyl/ProcessTree.h"
#include <Bpp/Seq/Container/AlignedValuesContainer.h>
#include <Bpp/Graph/AssociationDAGraphImplObserver.h>

namespace bpp
{
  namespace dataflow
  {

    inline MatrixDimension transitionMatrixDimension (std::size_t nbState) {
      return {Eigen::Index (nbState), Eigen::Index (nbState)};
    }

    inline MatrixDimension equilibriumFrequenciesDimension (std::size_t nbState) {
      return rowVectorDimension (Eigen::Index (nbState));
    }

    inline MatrixDimension conditionalLikelihoodDimension (std::size_t nbState, std::size_t nbSite) {
    return {Eigen::Index (nbState), Eigen::Index (nbSite)};
    }
    
/** conditionalLikelihood = f(forwardLikelihood[children[i]] for i).
     * conditionalLikelihood: Matrix(state, site).
     * forwardLikelihood[i]: Matrix(state, site).
     *
     * c(state, site) = prod_i f_i(state, site).
     * Using member wise multiply: c = prod_member_i f_i.
     */
    
    using SpeciationFromChildrenForward = CWiseMul<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

    /** conditionalLikelihood = f(forwardLikelihood[children[i]] for i, prop).
     * conditionalLikelihood: Matrix(state, site).
     * forwardLikelihood[i]: Matrix(state, site).
     * prop: Vector([i])
     *
     * c(state, site) = sum_i f_i(state, site) * prop_i
     * Using member wise addition: c = sum_member_i (f_i * prop_i)
     */
    
    using MixtureFromChildrenForward = CWiseMean<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>, ReductionOf<double>>;


    /** @brief forwardLikelihood = f(transitionMatrix, conditionalLikelihood).
     * - forwardLikelihood: Matrix(state, site).
     * - transitionMatrix: Matrix (fromState, toState)
     * - conditionalLikelihood: Matrix(state, site).
     *
     * f(toState, site) = sum_fromState P(fromState, toState) * c(fromState, site).
     * Using matrix multiply with transposition: f = transposed(transitionMatrix) * c.
     */
    
    using ForwardLikelihoodFromConditional =
      MatrixProduct<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>;

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

    using ConditionalLikelihoodForward = Value<Eigen::MatrixXd>;
    using ConditionalLikelihoodForwardRef = ValueRef<Eigen::MatrixXd>;

    using ForwardLikelihoodBelow = Value<Eigen::MatrixXd>;
    using ForwardLikelihoodBelowRef = ValueRef<Eigen::MatrixXd>;

    class ForwardLikelihoodTree : public AssociationDAGlobalGraphObserver<ConditionalLikelihoodForward,ForwardLikelihoodBelow>
    {
      using DAClass = AssociationDAGlobalGraphObserver<ConditionalLikelihoodForward,ForwardLikelihoodBelow>;

    private:
      Context& context_;
      std::shared_ptr<ProcessTree> processTree_;
      MatrixDimension likelihoodMatrixDim_;
      const StateMap& statemap_;
      std::size_t nbState_;
      std::size_t nbSites_;

      /* Maps of edges and nodes between ForwardLikelihoodTree and ProcessTree */
      
      std::map<ForwardLikelihoodBelowRef, ProcessEdgeRef> mapEdge_;
      std::map<ConditionalLikelihoodForwardRef, ProcessNodeRef> mapNode_;

      /* Map of the indexes between species tree and likelihood tree */

      using DAGindexes = std::vector<uint>;
      using Speciesindex = uint;
      
      std::map<Speciesindex, DAGindexes> mapIndexes_;
      
    public:
      
      ForwardLikelihoodTree(Context& c, 
                            std::shared_ptr<ProcessTree> tree,
                            const StateMap& statemap) :
        DAClass(),
        context_(c), processTree_(tree), likelihoodMatrixDim_(), statemap_(statemap), nbState_(statemap.getNumberOfModelStates()), nbSites_(0), mapEdge_(), mapNode_()
      {
      }

      void initialize(const AlignedValuesContainer& sites)
      {
        nbSites_ = sites.getNumberOfSites (); 
        likelihoodMatrixDim_ = conditionalLikelihoodDimension (nbState_, nbSites_);
        ConditionalLikelihoodForwardRef bidonRoot=ConstantZero<Eigen::MatrixXd>::create(context_, MatrixDimension(1,1));
        createNode(bidonRoot);
        rootAt(bidonRoot); // for construction, temporary top node for new edges
        auto n=makeConditionalLikelihoodNode (processTree_->getRoot(), sites);
        rootAt(n);
        outputToDot("forwardRoot_"+TextTools::toString(getNodeIndex(n))+".dot","forwardTree");
        deleteNode(bidonRoot);

        // Now map the species indexes and the likelihood DAG indexes

        setSpeciesMapIndexes_();
      }

    private:

      ConditionalLikelihoodForwardRef makeInitialConditionalLikelihood (const std::string & sequenceName, const AlignedValuesContainer & sites);

      ForwardLikelihoodBelowRef makeForwardLikelihoodEdge (std::shared_ptr<ProcessEdge> edge, const AlignedValuesContainer & sites);
      
      ConditionalLikelihoodForwardRef makeConditionalLikelihoodNode (std::shared_ptr<ProcessNode> node, const AlignedValuesContainer & sites);

      /*
       * @brief Map the species indexes and the likelihood DAG
       * indexes.
       *
       */
      
      void setSpeciesMapIndexes_();
      
    public:
      
      /*
       * @brief Links with processTree edges & nodes
       *
       */

      ProcessNodeRef getProcessNode(const ConditionalLikelihoodForwardRef& node) const
      {
        return mapNode_.at(node);
      }

      ProcessEdgeRef getProcessEdge(const ForwardLikelihoodBelowRef& edge) const
      {
        return mapEdge_.at(edge);
      }

      /*
       * @brief get the forward likehoodarray at a given node in the DAG,
       * which number may not be the number in the tree.
       *
       * @param nodeId : index of the node in the likelihood DAG.
       */
    
      const ValueRef<Eigen::MatrixXd> getForwardLikelihoodArray(int nodeId) const
      {
        return getNode(nodeId);
      }

      const ValueRef<Eigen::MatrixXd> getForwardLikelihoodArrayAtRoot() const
      {
        return getRoot();
      }


      friend class LikelihoodCalculationSingleProcess;

    };
  }
  
} //end of namespace bpp.

#endif //_FORWARD_LIKELIHOOD_TREE_H_

