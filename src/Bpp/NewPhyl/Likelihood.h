//
// File: Likelihood.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-03
// Last modified: 2017-05-03
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#pragma once
#ifndef BPP_NEWPHYL_LIKELIHOOD_H
#define BPP_NEWPHYL_LIKELIHOOD_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h> // For dimension types, TODO remove for local forward declaration ?
#include <functional>
#include <unordered_map>

namespace bpp {
  /* Conditional likelihoods are stored in a matrix of sizes (nbState, nbSite).
   * Rows represents states (nucleotides, proteins or codon).
   * Columns represents sites (one site for each column).
   * Conditional likelihood is thus accessed by m(state,site) for an eigen matrix.
   *
   * A Transition matrix is a (nbState,nbState) matrix.
   * tm(toState, fromState) = probability of going to toState from fromState.
   *
   * Equilibrium frequencies are stored as a RowVector : matrix with 1 row and n columns.
   * This choice allows to reuse the MatrixProduct numeric node directly.
   */
  inline MatrixDimension conditionalLikelihoodDimension (std::size_t nbState, std::size_t nbSite) {
    return {Eigen::Index (nbState), Eigen::Index (nbSite)};
  }
  inline MatrixDimension transitionMatrixDimension (std::size_t nbState) {
    return {Eigen::Index (nbState), Eigen::Index (nbState)};
  }
  inline MatrixDimension equilibriumFrequenciesDimension (std::size_t nbState) {
    return rowVectorDimension (Eigen::Index (nbState));
  }

  /* Initial conditional likelihood for leaves (sequences on the tree) are computed separately.
   * The values should be used in the dataflow tree by using matrix constant nodes.
   */
  class Sequence; // TODO

  // Dataflow nodes for likelihood computation.
  namespace dataflow {
    /** conditionalLikelihood = product_i forwardLikelihood[children[i]].
     * conditionalLikelihood: Matrix.
     * forwardLikelihood[i]: Matrix.
     * c(state, site) = prod_i f_i(state, site)
     */
    using ConditionalLikelihoodFromChildrenForward = CWiseMul<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

    /** forwardLikelihood = transitionMatrix * conditionalLikelihood.
     * forwardLikelihood: Matrix.
     * conditionalLikelihood: Matrix.
     *
     * f(toState, site) = sum_fromState P(toState, fromState) * c(fromState, site).
     * Matrix multiply provides this exact computation efficiently.
     */
    using ForwardLikelihoodFromConditional = MatrixProduct<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>;

    /** likelihood = equilibriumFrequencies * rootConditionalLikelihood.
     * likelihood: RowVector.
     * equilibriumFrequencies: RowVector.
     * rootConditionalLikelihood: Matrix.
     *
     * likelihood(site) = sum_state equFreqs(state) * rootConditionalLikelihood(state, site).
     * By using RowVector, this computation is done by matrix multiply.
     */
    using LikelihoodFromRootConditional =
      MatrixProduct<Eigen::RowVectorXd, Eigen::RowVectorXd, Eigen::MatrixXd>;

    /** totalLogLikelihood = sum_site log(likelihood(site)).
     * likelihood: F (matrix-like type).
     * totalLogLikelihood: double.
     */
    using TotalLogLikelihood = SumOfLogarithms<Eigen::RowVectorXd>;
  } // namespace dataflow

  /* Likelihood transition model.
   */
  class TransitionModel;

  namespace dataflow {
    /** Helper: create a map with mutable dataflow nodes for each model parameter.
     * The map is indexed by model non namespaced names.
     */
    std::unordered_map<std::string, std::shared_ptr<NumericMutable<double>>>
    createParameterMapForModel (Context & c, const TransitionModel & model);

    /** Create a dependency vector suitable for a Model class constructor.
     * The vector is built from model parameter names, and an opaque accessor function.
     * For each named parameter in the model, getParameter(name) should return a valid node.
     * Only non-namespaced names are tried.
     * If no node is found (NodeRef was null), an exception is thrown.
     * Returned nodes must be Value<double> nodes.
     */
    NodeRefVec createDependencyVector (const TransitionModel & model,
                                       const std::function<NodeRef (const std::string &)> & getParameter);

    /** Data flow node representing a Model configured with parameter values.
     * This class wraps a bpp::TransitionModel as a data flow node.
     * It depends on Value<double> nodes (one for each parameter declared in the model).
     * It provides a dummy value representing the "model configured by its parameters".
     * This dummy value is then used by other node types to compute equilibrium frequencies,
     * transition matrices and their derivatives.
     *
     * The dummy value is implemented as a pointer to the internal model for simplicity.
     */
    class ConfiguredModel : public Value<const TransitionModel *> {
    public:
      using Self = ConfiguredModel;

      /** Create a new model node from a dependency vector.
       * Model parameters are given by a dependency vector of Value<double> nodes.
       * The number and order of parameters is given by the TransitionModel internal ParameterList.
       */
      static std::shared_ptr<ConfiguredModel> create (Context & c, NodeRefVec && deps,
                                                      std::unique_ptr<TransitionModel> && model);

      ConfiguredModel (NodeRefVec && deps, std::unique_ptr<TransitionModel> && model);
      ~ConfiguredModel ();

      /// Return the index of parameter with the given non namespaced name (or throw).
      std::size_t getParameterIndex (const std::string & name);
      /// Return the non namespaced name for parameter at the given index.
      const std::string & getParameterName (std::size_t index);

      std::string description () const final;
      std::string debugInfo () const final;

      /// Configuration for numerical derivation of computation nodes using this Model.
      NumericalDerivativeConfiguration config;

      // Recreate model node with other deps (TODO restore rebuild stuff ?)
      NodeRef recreate (Context & c, NodeRefVec && deps) const;

    private:
      void compute () final;

      std::unique_ptr<TransitionModel> model_;
    };

    /** equilibriumFrequencies = f(model).
     * equilibriumFrequencies: RowVector.
     * model: ConfiguredModel.
     */
    class EquilibriumFrequenciesFromModel : public Value<Eigen::RowVectorXd> {
    public:
      using Self = EquilibriumFrequenciesFromModel;
      using T = Eigen::RowVectorXd;

      static ValueRef<T> create (Context & c, NodeRefVec && deps, const Dimension<T> & dim);
      EquilibriumFrequenciesFromModel (NodeRefVec && deps, const Dimension<T> & dim);

      std::string debugInfo () const final;

      NodeRef derive (Context & c, const Node & node) final;

    private:
      void compute () final;

      Dimension<T> targetDimension_;
    };

    /** transitionMatrix = f(model, branchLen).
     * transitionMatrix: Matrix.
     * model: ConfiguredModel.
     * branchLen: double.
     */
    class TransitionMatrixFromModel : public Value<Eigen::MatrixXd> {
    public:
      using Self = TransitionMatrixFromModel;
      using T = Eigen::MatrixXd;

      static ValueRef<T> create (Context & c, NodeRefVec && deps, const Dimension<T> & dim);
      TransitionMatrixFromModel (NodeRefVec && deps, const Dimension<T> & dim);

      std::string debugInfo () const final;

      NodeRef derive (Context & c, const Node & node) final;

    private:
      void compute () final;

      Dimension<T> targetDimension_;
    };

  } // namespace dataflow
} // namespace bpp

#endif // BPP_NEWPHYL_LIKELIHOOD_H
