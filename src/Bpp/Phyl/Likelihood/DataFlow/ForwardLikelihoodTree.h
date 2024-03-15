// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_FORWARDLIKELIHOODTREE_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_FORWARDLIKELIHOODTREE_H

#include <Bpp/Graph/AssociationDAGraphImplObserver.h>
#include <Bpp/Seq/Container/AlignmentData.h>

#include "Bpp/Phyl/Likelihood/DataFlow/ProcessTree.h"
#include "Definitions.h"

namespace bpp
{
// using RowLik = Eigen::Matrix<double, 1, Eigen::Dynamic>;
// using MatrixLik = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;


inline Dimension<TransitionFunction> transitionFunctionDimension (Eigen::Index nbState)
{
  return Dimension<TransitionFunction>(nbState);
}

inline MatrixDimension conditionalLikelihoodDimension (Eigen::Index nbState, Eigen::Index nbSite)
{
  return {nbState, nbSite};
}

/** conditionalLikelihood = f(forwardLikelihood[children[i]] for i).
 * conditionalLikelihood: Matrix(state, site).
 * forwardLikelihood[i]: Matrix(state, site).
 *
 * c(state, site) = prod_i f_i(state, site).
 * Using member wise multiply: c = prod_member_i f_i.
 */

using SpeciationForward = CWiseMul<MatrixLik, ReductionOf<MatrixLik> >;

/** conditionalLikelihood = f(forwardLikelihood[children[i]] for i).
 * conditionalLikelihood: Matrix(state, site).
 * forwardLikelihood[i]: Matrix(state, site).
 *
 * c(state, site) = sum_i f_i(state, site)
 * Using member wise addition: c = sum_member_i f_i
 */

using MixtureForward = CWiseAdd<MatrixLik, ReductionOf<MatrixLik> >;

/** @brief forwardLikelihood = f(transitionMatrix, conditionalLikelihood).
 * - forwardLikelihood: Matrix(state, site).
 * - transitionMatrix: Matrix (fromState, toState)
 * - conditionalLikelihood: Matrix(state, site).
 *
 * f(toState, site) = sum_fromState P(fromState, toState) * c(fromState, site).
 */

using ForwardTransition =
  MatrixProduct<MatrixLik, Eigen::MatrixXd, MatrixLik>;

using ForwardTransitionFunction =
  CWiseApply<MatrixLik, MatrixLik, TransitionFunction>;

/** @brief forwardLikelihood = f(transitionMatrix, proportion).
 * - forwardLikelihood: Matrix(state, site).
 * - proportion: Double
 * - conditionalLikelihood: Matrix(state, site).
 *
 * f(State, site) = c(fromState, site) * prop
 */

using ForwardProportion =
  CWiseMul<MatrixLik, std::tuple<double, MatrixLik> >;

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

// Lower Conditional Likelihood under nodes
using ConditionalLikelihoodForward = Value<MatrixLik>;
using ConditionalLikelihoodForwardRef = ValueRef<MatrixLik>;

// Lower Likelihood at top of edges
using ForwardLikelihoodBelow = Value<MatrixLik>;
using ForwardLikelihoodBelowRef = ValueRef<MatrixLik>;

using DAGindexes = std::vector<uint>;
using Speciesindex = uint;


class ForwardLikelihoodTree : public AssociationDAGlobalGraphObserver<ConditionalLikelihoodForward, ForwardLikelihoodBelow>
{
  using DAClass = AssociationDAGlobalGraphObserver<ConditionalLikelihoodForward, ForwardLikelihoodBelow>;

private:
  Context& context_;
  std::shared_ptr<ProcessTree> processTree_;
  MatrixDimension likelihoodMatrixDim_;
  const StateMapInterface& statemap_;
  Eigen::Index nbState_;
  Eigen::Index nbSites_;

  /* Map of the indexes of nodes between species tree and
   * likelihood tree */

  std::map<Speciesindex, DAGindexes> mapNodesIndexes_; // For nodes that bring
  // information (ie not the empty ones)

  std::map<Speciesindex, DAGindexes> mapEdgesIndexes_; // For edges that bring
  // information (ie not the empty ones)

public:
  ForwardLikelihoodTree(Context& c,
                        std::shared_ptr<ProcessTree> tree,
                        const StateMapInterface& statemap) :
    DAClass(),
    context_(c), processTree_(tree), likelihoodMatrixDim_(), statemap_(statemap), nbState_(Eigen::Index(statemap.getNumberOfModelStates())), nbSites_(0)
  {}

  void initialize(const AlignmentDataInterface& sites)
  {
    nbSites_ = Eigen::Index(sites.getNumberOfSites ());
    likelihoodMatrixDim_ = conditionalLikelihoodDimension (nbState_, nbSites_);
    ConditionalLikelihoodForwardRef bidonRoot = ConstantZero<MatrixLik>::create(context_, MatrixDimension(1, 1));
    createNode(bidonRoot);
    /* Not sure it is necessary:

       setNodeIndex(bidonRoot,processTree_->getNodeIndex(processTree_->getRoot()));
     */

    rootAt(bidonRoot); // for construction, temporary top node for new edges
    auto n = makeForwardLikelihoodAtNode (processTree_->getRoot(), sites);
    rootAt(n);
    deleteNode(bidonRoot);
  }

private:
  /**
   * @brief Compute ConditionalLikelihood after reading edge on
   * the forward proces (ie at top of the edge).
   */
  ForwardLikelihoodBelowRef makeForwardLikelihoodAtEdge(
      std::shared_ptr<ProcessEdge> edge,
      const AlignmentDataInterface& sites);

  /*
   * @brief Compute ConditionalLikelihood after reading node on
   * the forward proces (ie just above node).
   */
  ConditionalLikelihoodForwardRef makeForwardLikelihoodAtNode(
      std::shared_ptr<ProcessNode> node,
      const AlignmentDataInterface& sites);

  /**
   * @brief Compute ConditionalLikelihood for leaf.
   */
  ConditionalLikelihoodForwardRef makeInitialConditionalLikelihood(
      const std::string& sequenceName,
      const AlignmentDataInterface& sites);

  /**
   * @brief Map the species indexes and the likelihood DAG
   * indexes.
   */
  void setSpeciesMapIndexes_();

protected:
  Context& getContext()
  {
    return context_;
  }

public:
  /*
   * @brief Links with processTree edges & nodes
   *
   */
  std::shared_ptr<const ProcessTree>  getProcessTree() const
  {
    return processTree_;
  }

  /*
   * @brief Get the nodes indexes of the DAG that correspond to
   * the species Index (of the Process tree).
   */
  const DAGindexes& getDAGNodesIndexes(const Speciesindex speciesIndex) const
  {
    return mapNodesIndexes_.at(speciesIndex);
  }

  /*
   * @brief Get the edges indexes of the DAG that correspond to
   * the species Index (of the Process tree).
   */
  
  const DAGindexes& getDAGEdgesIndexes(const Speciesindex speciesIndex) const
  {
    return mapEdgesIndexes_.at(speciesIndex);
  }

  /*
   * @brief get the forward likehoodarray at a given node in the DAG,
   * which number may not be the number in the tree.
   *
   * @param nodeId : index of the node in the likelihood DAG.
   */
  const ValueRef<MatrixLik> getForwardLikelihoodArray(uint nodeId) const
  {
    return getNode(nodeId);
  }

  const ValueRef<MatrixLik> getForwardLikelihoodArrayAtRoot() const
  {
    return getRoot();
  }

  friend class LikelihoodCalculationSingleProcess;
  friend class ProbabilityDAG;
};

/**
 *@brief DAG with the same shape as ForwardLikelihoodTree with
 * computations of probabilities of nodes & branches.
 *
 */


using Proba = Value<double>;
using ProbaRef = std::shared_ptr<Value<double> >;

using DAProb = AssociationDAGlobalGraphObserver<Proba, Proba>;

using ProbaMul = CWiseMul<double, std::tuple<double, double> >;

using ProbaSum = CWiseAdd<double, ReductionOf<double> >;

class ProbabilityDAG :
  public DAProb
{
private:
  Context& context_;

public:
  ProbabilityDAG(std::shared_ptr<ForwardLikelihoodTree> tree);

public:
  double getProbaAtNode(PhyloTree::NodeIndex nodeId)
  {
    return getNode(nodeId)->targetValue();
  }

  double getProbaAtEdge(PhyloTree::EdgeIndex edgeId)
  {
    return getEdge(edgeId)->targetValue();
  }

private:
  /**
   * @brief computation of the probabilities with the same approach
   * as for BackwardLikelihoodTree.
   */
  ProbaRef makeProbaAtEdge_(PhyloTree::EdgeIndex edgeIndex, std::shared_ptr<ForwardLikelihoodTree> tree);

  ProbaRef makeProbaAtNode_(PhyloTree::EdgeIndex edgeIndex, std::shared_ptr<ForwardLikelihoodTree> tree);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_FORWARDLIKELIHOODTREE_H
