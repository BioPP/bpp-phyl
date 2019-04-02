//
// File: ForwardLikelihoodTree.h
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

#include "Bpp/NewPhyl/PhyloTree_BrRef.h"
#include "Bpp/NewPhyl/Model.h"
#include "Bpp/NewPhyl/DiscreteDistribution.h"
#include "Bpp/NewPhyl/FrequenciesSet.h"

namespace bpp
{
  namespace dataflow
  {

    /** conditionalLikelihood = f(forwardLikelihood[children[i]] for i).
     * conditionalLikelihood: Matrix(state, site).
     * forwardLikelihood[i]: Matrix(state, site).
     *
     * c(state, site) = prod_i f_i(state, site).
     * Using member wise multiply: c = prod_member_i f_i.
     */
    
    using ConditionalLikelihoodFromChildrenForward = CWiseMul<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

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

    using ForwardLikelihoodBelow = Value<Eigen::MatrixXd>;

    class ForwardLikelihoodTree : public AssociationTreeGlobalGraphObserver<ConditionalLikelihoodForward,ForwardLikelihoodBelow>
    {
      using TreeClass = AssociationTreeGlobalGraphObserver<ConditionalLikelihoodForward,ForwardLikelihoodBelow>;
      
    private:
      dataflow::Context& context_;
      std::shared_ptr<dataflow::PhyloTree_BrRef> tree_;
      MatrixDimension likelihoodMatrixDim_;
      const StateMap& statemap_;
      std::size_t nbState_;
      std::size_t nbSite_;

    public:

      ForwardLikelihoodTree(dataflow::Context& c, 
                            std::shared_ptr<dataflow::PhyloTree_BrRef> tree,
                            const StateMap& statemap) :
        TreeClass(tree->getGraph()),
        context_(c), tree_(tree), likelihoodMatrixDim_(), statemap_(statemap), nbState_(statemap.getNumberOfModelStates()), nbSite_(0)
      {
      }

      void initialize(const AlignedValuesContainer& sites)
      {
        nbSite_ = sites.getNumberOfSites (); 
        likelihoodMatrixDim_ = conditionalLikelihoodDimension (nbState_, nbSite_);
        makeConditionalLikelihoodNode (tree_->getRootIndex (), sites);
      }
        
      dataflow::ValueRef<Eigen::MatrixXd> makeInitialConditionalLikelihood (const std::string & sequenceName, const AlignedValuesContainer & sites)
      {
        const auto sequenceIndex = sites.getSequencePosition (sequenceName);
        Eigen::MatrixXd initCondLik (nbState_, nbSite_);
        for (std::size_t site = 0; site < nbSite_; ++site) {
          for (std::size_t state = 0; state < nbState_; ++state) {
            initCondLik (Eigen::Index (state), Eigen::Index (site)) =
              sites.getStateValueAt (site, sequenceIndex, statemap_.getAlphabetStateAsInt(state));
          }
        }
        auto r=dataflow::NumericConstant<Eigen::MatrixXd>::create (context_, std::move (initCondLik));
        
        return r;
      }

      dataflow::ValueRef<Eigen::MatrixXd> makeForwardLikelihoodNode (PhyloTree::EdgeIndex index, const AlignedValuesContainer & sites) {
        const auto brlen = tree_->getEdge(index)->getBrLen();
        const auto model= tree_->getEdge(index)->getModel();
        auto childConditionalLikelihood = makeConditionalLikelihoodNode (tree_->getSon(index), sites);
        auto transitionMatrix =
          dataflow::ConfiguredParametrizable::createMatrix<dataflow::ConfiguredModel, dataflow::TransitionMatrixFromModel> (context_, {model, brlen}, transitionMatrixDimension (nbState_));
        auto r= dataflow::ForwardLikelihoodFromConditional::create (
          context_, {transitionMatrix, childConditionalLikelihood}, likelihoodMatrixDim_);

        associateEdge(r,index);
        setEdgeIndex(r, index);
        return r;
      }

      dataflow::ValueRef<Eigen::MatrixXd> makeConditionalLikelihoodNode (PhyloTree::NodeIndex index, const AlignedValuesContainer & sites) {
        const auto childBranchIndexes = tree_->getBranches (index);
        if (childBranchIndexes.empty ()) {
          auto r=makeInitialConditionalLikelihood (tree_->getNode (index)->getName (), sites);
          associateNode(r, index);
          setNodeIndex(r, index);
          return(r);
        }
        else {
          dataflow::NodeRefVec deps (childBranchIndexes.size ());
          for (std::size_t i = 0; i < childBranchIndexes.size (); ++i) {
            deps[i] = makeForwardLikelihoodNode (childBranchIndexes[i], sites);
          }
          auto r= dataflow::ConditionalLikelihoodFromChildrenForward::create (context_, std::move (deps),
                                                                              likelihoodMatrixDim_);
          associateNode(r, index);
          setNodeIndex(r, index);
          return r;
        }
      }

    /*
     * @brief the LikehoodArrays
     *
     */
    
      const ValueRef<Eigen::MatrixXd> getForwardLikelihoodArray(int nodeId) const
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

#endif //_FORWARD_LIKELIHOOD_TREE_H_

