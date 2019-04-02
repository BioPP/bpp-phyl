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

#include "Bpp/NewPhyl/PhyloTree_BrRef.h"
#include "Bpp/NewPhyl/Model.h"
//#include "Bpp/NewPhyl/FrequenciesSet.h"

#include "Bpp/NewPhyl/ForwardLikelihoodTree.h"

namespace bpp
{
  namespace dataflow
  {
    
    /** @brief : At the top of edges:
     *
     * conditionalLikelihood = f(backwardLikelihood[father[i], brothers[i]]).
     * conditionalLikelihood: Matrix(state, site).
     * forwardLikelihood[i]: Matrix(state, site).
     *
     * c(state, site) = prod_i f_i(state, site).
     * Using member wise multiply: c = prod_member_i f_i.
     */
    
    using ConditionalLikelihoodFromBrothersBackward = CWiseMul<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

    /** @brief : Above each node:
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
      MatrixProduct<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>;
    
    using ConditionalLikelihood = Value<Eigen::MatrixXd>;

    using BackwardLikelihoodAbove = Value<Eigen::MatrixXd>;

    /** Tree structure for all the forward computations **/
    /* All the computations are set in a DataFlow context */
    
    class BackwardLikelihoodTree : public AssociationTreeGlobalGraphObserver<ConditionalLikelihood,BackwardLikelihoodAbove>
    {
      using TreeClass = AssociationTreeGlobalGraphObserver<ConditionalLikelihood,BackwardLikelihoodAbove>;
      
    private:
      dataflow::Context& context_;
      std::shared_ptr<dataflow::ForwardLikelihoodTree> forwardTree_;
      std::shared_ptr<dataflow::PhyloTree_BrRef> tree_;
      MatrixDimension likelihoodMatrixDim_;
      const StateMap& statemap_;
      std::size_t nbState_;
      std::size_t nbSite_;

    public:

      BackwardLikelihoodTree(dataflow::Context& c, 
                             std::shared_ptr<dataflow::ForwardLikelihoodTree> forwardTree,
                             std::shared_ptr<dataflow::PhyloTree_BrRef> tree,
                             const StateMap& statemap,
                             std::size_t nbSite) :
        TreeClass(tree->getGraph()),
        context_(c), forwardTree_(forwardTree), tree_(tree), likelihoodMatrixDim_(conditionalLikelihoodDimension (nbState_, nbSite_)), statemap_(statemap), nbState_(statemap.getNumberOfModelStates()), nbSite_(nbSite)
      {}

      void setForwardTree(std::shared_ptr<dataflow::ForwardLikelihoodTree> forwardTree)
      {
        if (forwardTree_)
          throw Exception("BackwardLikelihoodTree::setForwardTree : forwardTree_ already set.");
        
        forwardTree_=forwardTree;
      }

      void setRootFrequencies(const dataflow::ValueRef<Eigen::RowVectorXd> rootFreqs)
      {
        auto r2=convertRef<bpp::dataflow::Value<Eigen::MatrixXd>, bpp::dataflow::Value<Eigen::RowVectorXd>>(rootFreqs);
        
        associateNode(r2, tree_->getRootIndex());
        setNodeIndex(r2, tree_->getRootIndex());
      }
        
      dataflow::ValueRef<Eigen::MatrixXd> makeBackwardAboveLikelihoodEdge (PhyloTree::EdgeIndex index) {

        auto fatherIndex = getFatherOfEdge(index);
        auto backNode = getNode(fatherIndex);
        
        if (!backNode)
          backNode = makeConditionalAboveLikelihoodNode(fatherIndex);

        // get forward likelihoods of brothers
        if (!forwardTree_)
          throw Exception("BackwardLikelihoodTree::makeBackwardAboveLikelihoodEdge: forwardTree_ is missing.");
        
        auto edgeIds = forwardTree_->getBranches(fatherIndex);
        dataflow::NodeRefVec deps;
        for (auto eId : edgeIds)
        {
          if (eId!=index)
            deps.push_back(forwardTree_->getEdge(eId));
        }
        deps.push_back(backNode);
        
        auto r= dataflow::ConditionalLikelihoodFromBrothersBackward::create (context_, std::move (deps),
                                                                             likelihoodMatrixDim_);
        
        associateEdge(r,index);
        setEdgeIndex(r, index);
        return r;
      }

      dataflow::ValueRef<Eigen::MatrixXd> makeConditionalAboveLikelihoodNode (PhyloTree::NodeIndex index) {
        //!!!! must be initialized before
        if (index==tree_->getRootIndex())
        {
          auto root=getNode(index);
          if (!root)
            throw Exception("BackwardLikelihoodTree::makeConditionalAboveLikelihoodNode : root is missing.");
          return root;
        }
        
        // get Edge with model
        const auto edgeToFather = tree_->getEdgeToFather(index);
        // get/check if edge with backward likelihood exists
        auto backEdge = getEdgeToFather(index);

        if (!backEdge)
          backEdge=makeBackwardAboveLikelihoodEdge(index);
        
        const auto brlen = edgeToFather->getBrLen();
        const auto model= edgeToFather->getModel();
        // useless, the transitionMatrix already exists, but is not available directly through a tree
        auto transitionMatrix =
          dataflow::ConfiguredParametrizable::createMatrix<dataflow::ConfiguredModel, dataflow::TransitionMatrixFromModel> (context_, {model, brlen}, transitionMatrixDimension (nbState_));
        auto r= dataflow::BackwardLikelihoodFromConditional::create (
          context_, {transitionMatrix, backEdge}, likelihoodMatrixDim_);

        associateNode(r, index);
        setNodeIndex(r, index);
        return r;
      }

    /*
     * @brief the LikehoodArrays
     *
     */
    
      const ValueRef<Eigen::MatrixXd> getBackwardLikelihoodArray(int nodeId) const
      {
        return getNode(nodeId);
      }

      /**
       * NbSite
       *
       */

      std::size_t getNumberOfSites() const 
      {
        return nbSite_;
      }
      
    };
  }
  
} //end of namespace bpp.

#endif //_BACKWARD_LIKELIHOOD_TREE_H_

