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

#include "Bpp/Phyl/NewLikelihood/DataFlow/ProcessTree.h"
#include <Bpp/Seq/Container/AlignedValuesContainer.h>
#include <Bpp/Graph/AssociationDAGraphImplObserver.h>

namespace bpp
{
  using RowLik = Eigen::Matrix<double, 1, Eigen::Dynamic>;
  using MatrixLik = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  

  inline Dimension<TransitionFunction> transitionFunctionDimension (Eigen::Index nbState) {
    return Dimension<TransitionFunction>(nbState, nbState);
  }

  inline MatrixDimension conditionalLikelihoodDimension (Eigen::Index nbState, Eigen::Index nbSite) {
    return {nbState, nbSite};
  }
    
  /** conditionalLikelihood = f(forwardLikelihood[children[i]] for i).
   * conditionalLikelihood: Matrix(state, site).
   * forwardLikelihood[i]: Matrix(state, site).
   *
   * c(state, site) = prod_i f_i(state, site).
   * Using member wise multiply: c = prod_member_i f_i.
   */
    
  using SpeciationForward = CWiseMul<MatrixLik, ReductionOf<MatrixLik>>;

  /** conditionalLikelihood = f(forwardLikelihood[children[i]] for i).
   * conditionalLikelihood: Matrix(state, site).
   * forwardLikelihood[i]: Matrix(state, site).
   *
   * c(state, site) = sum_i f_i(state, site)
   * Using member wise addition: c = sum_member_i f_i 
   */
    
  using MixtureForward = CWiseAdd<MatrixLik, ReductionOf<MatrixLik>>;

  /** @brief forwardLikelihood = f(transitionMatrix, conditionalLikelihood).
   * - forwardLikelihood: Matrix(state, site).
   * - transitionMatrix: Matrix (fromState, toState)
   * - conditionalLikelihood: Matrix(state, site).
   *
   * f(toState, site) = sum_fromState P(fromState, toState) * c(fromState, site).
   */
    
  using ForwardTransition =
    MatrixProduct<MatrixLik, MatrixLik, MatrixLik>;

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
    CWiseMul<MatrixLik, std::tuple<double, MatrixLik>>;

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

  using ConditionalLikelihoodForward = Value<MatrixLik>;
  using ConditionalLikelihoodForwardRef = ValueRef<MatrixLik>;

  using ForwardLikelihoodBelow = Value<MatrixLik>;
  using ForwardLikelihoodBelowRef = ValueRef<MatrixLik>;

  using DAGindexes = std::vector<uint>;
  using Speciesindex = uint;
      

  class ForwardLikelihoodTree : public AssociationDAGlobalGraphObserver<ConditionalLikelihoodForward,ForwardLikelihoodBelow>
  {
    using DAClass = AssociationDAGlobalGraphObserver<ConditionalLikelihoodForward,ForwardLikelihoodBelow>;

  private:
    Context& context_;
    std::shared_ptr<ProcessTree> processTree_;
    MatrixDimension likelihoodMatrixDim_;
    const StateMap& statemap_;
    Eigen::Index nbState_;
    Eigen::Index nbSites_;

    /* Map of the indexes of nodes between species tree and
     * likelihood tree */

    std::map<Speciesindex, DAGindexes> mapNodesIndexes_; // For nodes that bring
    //information (ie not the empty ones)

    std::map<Speciesindex, DAGindexes> mapEdgesIndexes_; // For edges that bring
    //information (ie not the empty ones)

  public:
      
    ForwardLikelihoodTree(Context& c, 
                          std::shared_ptr<ProcessTree> tree,
                          const StateMap& statemap) :
      DAClass(),
      context_(c), processTree_(tree), likelihoodMatrixDim_(), statemap_(statemap), nbState_(Eigen::Index(statemap.getNumberOfModelStates())), nbSites_(0)
    {
    }

    void initialize(const AlignedValuesContainer& sites)
    {
      nbSites_ = Eigen::Index(sites.getNumberOfSites ()); 
      likelihoodMatrixDim_ = conditionalLikelihoodDimension (nbState_, nbSites_);
      ConditionalLikelihoodForwardRef bidonRoot=ConstantZero<MatrixLik>::create(context_, MatrixDimension(1,1));
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

    /*
     * @brief Compute ConditionalLikelihood after reading edge on
     * the forward proces (ie at top of the edge).
     *
     */
      
    ForwardLikelihoodBelowRef makeForwardLikelihoodAtEdge (std::shared_ptr<ProcessEdge> edge, const AlignedValuesContainer & sites);

    /*
     * @brief Compute ConditionalLikelihood after reading node on
     * the forward proces (ie just above node).
     *
     */

    ConditionalLikelihoodForwardRef makeForwardLikelihoodAtNode (std::shared_ptr<ProcessNode> node, const AlignedValuesContainer & sites);

    /*
     * @brief Compute ConditionalLikelihood for leaf.
     *
     */

    ConditionalLikelihoodForwardRef makeInitialConditionalLikelihood (const std::string & sequenceName, const AlignedValuesContainer & sites);
      
    /*
     * @brief Map the species indexes and the likelihood DAG
     * indexes.
     *
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
   *computations of probabilities of nodes & branches.
   *
   */

  
  using Proba = Value<double>;
  using ProbaRef = std::shared_ptr<Value<double> >;
  
  using DAProb = AssociationDAGlobalGraphObserver<Proba, Proba>;

  using ProbaMul = CWiseMul<double, std::tuple<double, double> >;
  
  using ProbaSum = CWiseAdd<double, ReductionOf<double>>;
  
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
      return getNode(nodeId)->getTargetValue();
    }
    
    double getProbaAtEdge(PhyloTree::EdgeIndex edgeId)
    {
      return getEdge(edgeId)->getTargetValue();
    }

  private:
    /*
     *@brief computation of the probabilities with the same approach
     *as for BackwardLikelihoodTree.
     *
     */
      
    ProbaRef makeProbaAtEdge_(PhyloTree::EdgeIndex edgeIndex, std::shared_ptr<ForwardLikelihoodTree> tree);
    
    ProbaRef makeProbaAtNode_(PhyloTree::EdgeIndex edgeIndex, std::shared_ptr<ForwardLikelihoodTree> tree);


  };
  

  
} //end of namespace bpp.

#endif //_FORWARD_LIKELIHOOD_TREE_H_

