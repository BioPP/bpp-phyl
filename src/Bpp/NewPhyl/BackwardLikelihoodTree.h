//
// File: BackwardLikelihoodTree.h
// Created by: Laurent Guéguen
// Created on: mercredi 12 décembre 2018, à 17h 00
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

#ifndef _BACKWARD_LIKELIHOOD_TREE_H_
#define _BACKWARD_LIKELIHOOD_TREE_H_

#include "Bpp/NewPhyl/ProcessTree.h"
#include "Bpp/NewPhyl/Model.h"

#include "Bpp/NewPhyl/ForwardLikelihoodTree.h"

namespace bpp
{
  namespace dataflow
  {
    
    /** @brief : At the top of each edge at speciation node
     *
     * conditionalLikelihood = f(backwardLikelihood[father[i]], forwardlikelihood[brothers[i]]).
     * conditionalLikelihood: Matrix(state, site).
     * backwardLikelihood, forwardlikelihood: Matrix(state, site).
     *
     * c(state, site) = prod_i f_i(state, site).
     * Using member wise multiply: c = prod_member_i f_i.
     */
    
    using ConditionalLikelihoodFromUpper = CWiseMul<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

    /** @brief : Above each node : bottom of each edge
     *
     *  backwardLikelihood = f(transitionMatrix, conditionalLikelihood).
     *  
     * - backwardLikelihood: Matrix(state, site).
     * - transitionMatrix: Matrix (FromState, toState)
     * - conditionalLikelihood: Matrix(state, site).
     *
     * f(fromState, site) = sum_toState P(fromState, toState) * c(toState, site).
     * Using matrix multiply: f = transitionMatrix * c.
     */
    
    using BackwardLikelihoodFromConditional =
      MatrixProduct<Eigen::MatrixXd, Transposed<Eigen::MatrixXd>, Eigen::MatrixXd>;

    /** @brief : Above each node : in case of mixture of above edges
     *
     * mixturefromupper = f(backwardLikelihood[edgestofathers[i]] for i, prop[i]).
     * edgestofathers: Matrix(state, site).
     * prob[i]: probability of edge leading to father[i]
     *
     * c(state, site) = sum_i f_i(state, site) * prop_i
     * Using member wise addition: c = sum_member_i (f_i * prop_i)
     */
    
    using MixtureFromUpperBackward = CWiseMean<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>, ReductionOf<double>>;

    // Upper Likelihood in nodes
    using ConditionalLikelihood = Value<Eigen::MatrixXd>;
    using ConditionalLikelihoodRef = ValueRef<Eigen::MatrixXd>;

    // Upper Likelihood at top of edges
    using BackwardLikelihoodAbove = Value<Eigen::MatrixXd>;
    using BackwardLikelihoodAboveRef = ValueRef<Eigen::MatrixXd>;

//    using FullLikelihood = Value<Patterned<Eigen::MatrixXd>>;

    /** Tree structure for all the forward computations **/
    /* All the computations are set in a DataFlow context */

    class BackwardLikelihoodTree : public AssociationDAGlobalGraphObserver<ConditionalLikelihood,BackwardLikelihoodAbove>
    {
      using DAClass = AssociationDAGlobalGraphObserver<ConditionalLikelihood,BackwardLikelihoodAbove>;

      /** For a given rate catagory, stores ProcessTree,
       * ForwardlikelihoodTree and BackwardLikelihoodTree
       **/
      
    private:
      dataflow::Context& context_;
      std::size_t nbState_;
      std::size_t nbSite_;
      std::shared_ptr<dataflow::ForwardLikelihoodTree> forwardTree_;
      std::shared_ptr<dataflow::ProcessTree> processTree_;
      ValueRef<Eigen::RowVectorXd> rFreqs_;
      MatrixDimension likelihoodMatrixDim_;
      const StateMap& statemap_;

    public:

      BackwardLikelihoodTree(dataflow::Context& c, 
                             std::shared_ptr<dataflow::ForwardLikelihoodTree> forwardTree,
                             std::shared_ptr<dataflow::ProcessTree> tree,
                             ValueRef<Eigen::RowVectorXd> rFreqs,
                             const StateMap& statemap,
                             std::size_t nbSite) :
        DAClass(forwardTree->getGraph()),
        context_(c), nbState_(statemap.getNumberOfModelStates()), nbSite_(nbSite), forwardTree_(forwardTree), processTree_(tree), rFreqs_(rFreqs), likelihoodMatrixDim_(conditionalLikelihoodDimension (nbState_, nbSite_)), statemap_(statemap)
      {
      }

      ConditionalLikelihoodRef setRootFrequencies(const dataflow::ValueRef<Eigen::RowVectorXd> rootFreqs)
      {
        auto r2=bpp::dataflow::CWiseFill<Eigen::MatrixXd, Eigen::RowVectorXd>::create(context_, {rootFreqs}, likelihoodMatrixDim_);

        associateNode(r2, forwardTree_->getRootIndex());
        setNodeIndex(r2, forwardTree_->getRootIndex());
        return r2;
      }

    private:
      BackwardLikelihoodAboveRef makeBackwardAboveLikelihoodEdge (PhyloTree::EdgeIndex index);
      
      ConditionalLikelihoodRef makeConditionalAboveLikelihoodNode (PhyloTree::NodeIndex index);

      /*
       * @brief the LikehoodArrays
       *
       * Beware: nodeIds are in the DAG, not the ids of the PhyloTree.
       *
       * Set in private to avoid bad usage, access through
       * LikelihoodCalculationSingleProcess.
       */
      
      ConditionalLikelihoodRef getBackwardLikelihoodArray(int nodeId)
      {
        if (!hasNode(nodeId))
          makeConditionalAboveLikelihoodNode(nodeId);
        
        return getNode(nodeId);
      }

      friend class LikelihoodCalculationSingleProcess;

    };
  }
  
} //end of namespace bpp.

#endif //_BACKWARD_LIKELIHOOD_TREE_H_

