// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_BACKWARDLIKELIHOODTREE_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_BACKWARDLIKELIHOODTREE_H


#include "Bpp/Phyl/Likelihood/DataFlow/ForwardLikelihoodTree.h"
#include "Bpp/Phyl/Likelihood/DataFlow/Model.h"
#include "Bpp/Phyl/Likelihood/DataFlow/ProcessTree.h"

namespace bpp
{
// using RowLik = Eigen::Matrix<double, 1, Eigen::Dynamic>;
// using MatrixLik = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;


/** @brief : At the top of each edge below a speciation node
 *
 * conditionalLikelihood = f(backwardLikelihood[father[i]], forwardlikelihood[brothers[i]]).
 * conditionalLikelihood: Matrix(state, site).
 * backwardLikelihood, forwardlikelihood: Matrix(state, site).
 *
 * Using member wise multiplication: c(state, site) = prod_i f_i(state, site).
 */

using SpeciationBackward = CWiseMul<MatrixLik, ReductionOf<MatrixLik>>;

/** @brief : At the top of each edge below a mixture node
 *
 * backwardLikelihood = f(backwardLikelihood[father[i]] for i)
 * backwardLikelihood: Matrix(state, site).
 * backwardLikelihood[i]: Matrix(state, site).
 *
 * c(state, site) = sum_i f_i(state, site)
 * Using member wise addition: c = sum_member_i f_i
 * Using member wise multiplication: c(state, site) = prod_i f_i(state, site).
 */

using MixtureBackward = CWiseAdd<MatrixLik, ReductionOf<MatrixLik>>;

/** @brief : Above each node : bottom of an edge in case of transition from upper
 *
 *  backwardLikelihood = f(transitionMatrix, conditionalLikelihood).
 *
 * - backwardLikelihood: Matrix(state, site).
 * - transitionMatrix: Matrix (FromState, toState)
 * - conditionalLikelihood: Matrix(state, site).
 *
 * f(fromState, site) = sum_toState P(fromState, toState) * c(toState, site).
 * Using matrix multiply: f = transposed(transitionMatrix) * c.
 */

using BackwardTransition =
    MatrixProduct<MatrixLik, Transposed<Eigen::MatrixXd>, MatrixLik>;

/** @brief : Above each node : in case of mixture of above edges
 *
 * - backwardlikelihood: Matrix(state, site).
 * - proportion: Double
 * - conditionalLikelihood: Matrix(state, site).
 *
 * f(State, site) = c(fromState, site) * prop
 */

using BackwardProportion = CWiseMul<MatrixLik, std::tuple<double, MatrixLik>>;

// Upper Likelihood in nodes
using ConditionalLikelihood = Value<MatrixLik>;
using ConditionalLikelihoodRef = ValueRef<MatrixLik>;

// Upper Likelihood at top of edges
using BackwardLikelihoodAbove = Value<MatrixLik>;
using BackwardLikelihoodAboveRef = ValueRef<MatrixLik>;

//    using FullLikelihood = Value<Patterned<MatrixLik>>;

/** Tree structure for all the forward computations **/
/* All the computations are set in a DataFlow context */

class BackwardLikelihoodTree : public AssociationDAGlobalGraphObserver<ConditionalLikelihood, BackwardLikelihoodAbove>
{
  using DAClass = AssociationDAGlobalGraphObserver<ConditionalLikelihood, BackwardLikelihoodAbove>;

  /** For a given rate catagory, stores ProcessTree,
   * ForwardlikelihoodTree and BackwardLikelihoodTree
   **/

private:
  Context& context_;
  Eigen::Index nbState_;
  Eigen::Index nbSite_;
  std::shared_ptr<ForwardLikelihoodTree> forwardTree_;
  std::shared_ptr<ProcessTree> processTree_;
  ValueRef<Eigen::RowVectorXd> rFreqs_;
  MatrixDimension likelihoodMatrixDim_;
  const StateMapInterface& statemap_;

public:
  BackwardLikelihoodTree(Context& c,
      std::shared_ptr<ForwardLikelihoodTree> forwardTree,
      std::shared_ptr<ProcessTree> tree,
      ValueRef<Eigen::RowVectorXd> rFreqs,
      const StateMapInterface& statemap,
      Eigen::Index nbSite) :
    DAClass(forwardTree->getGraph()),
    context_(c), nbState_(Eigen::Index(statemap.getNumberOfModelStates())), nbSite_(nbSite), forwardTree_(forwardTree), processTree_(tree), rFreqs_(rFreqs), likelihoodMatrixDim_(conditionalLikelihoodDimension (nbState_, nbSite_)), statemap_(statemap)
  {}

  ConditionalLikelihoodRef setRootFrequencies(const ValueRef<Eigen::RowVectorXd> rootFreqs)
  {
    auto rootFreqsEf = Convert<ExtendedFloatRowVectorXd, Eigen::RowVectorXd>::create(context_, {rootFreqs}, RowVectorDimension (Eigen::Index (nbState_)));
    auto r2 = CWiseFill<MatrixLik, RowLik>::create(context_, {rootFreqsEf}, likelihoodMatrixDim_);

    associateNode(r2, forwardTree_->getNodeGraphid(forwardTree_->getRoot()));
    setNodeIndex(r2, forwardTree_->getRootIndex());
    return r2;
  }

private:
  /*
   * @brief Compute joined likelihood BEFORE reading edge on the
   * backward process (ie at top of the edge).
   *
   */

  BackwardLikelihoodAboveRef makeBackwardLikelihoodAtEdge (PhyloTree::EdgeIndex index);

  /*
   * @brief Compute joined likelihood at node on the backward
   * process (ie with the bottom of the edges getting to the node
   * from upward).
   *
   */

  ConditionalLikelihoodRef makeBackwardLikelihoodAtNode (PhyloTree::NodeIndex index);

  /*
   * @brief the LikehoodArrays
   *
   * Beware: nodeIds are in the DAG, not the ids of the PhyloTree.
   *
   * Set in private to avoid bad usage, access through
   * LikelihoodCalculationSingleProcess.
   */
  ConditionalLikelihoodRef getBackwardLikelihoodArray(PhyloTree::NodeIndex nodeId)
  {
    if (!hasNode(nodeId))
      makeBackwardLikelihoodAtNode(nodeId);

    return getNode(nodeId);
  }

  friend class LikelihoodCalculationSingleProcess;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_BACKWARDLIKELIHOODTREE_H
