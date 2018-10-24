//
// File: Model.h
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
#ifndef BPP_NEWPHYL_MODEL_H
#define BPP_NEWPHYL_MODEL_H 

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/Exceptions.h>
#include <functional>
#include <unordered_map>

namespace bpp {
  /* Likelihood transition model.
   */
  class TransitionModel;

  namespace dataflow {

    /** @brief Data flow node representing a Model configured with parameter values.
     *
     * This class wraps a bpp::TransitionModel as a data flow node.
     * It depends on Value<double> nodes (one for each parameter declared in the model).
     * It provides a dummy value representing the "model configured by its parameters".
     * This dummy value is then used by other node types to compute equilibrium frequencies,
     * transition matrices and their derivatives.
     *
     * The dummy value is implemented as a pointer to the internal model for simplicity.
     */
    class ConfiguredModel : public Value<const TransitionModel*>// ,
                            // public ConfiguredParametrizable
    {
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

      std::string description () const final;
      std::string debugInfo () const final;

      /// Return the index of parameter with the given non namespaced name (or throw).
      std::size_t getParameterIndex (const std::string & name);
      /// Return the non namespaced name for parameter at the given index.
      const std::string & getParameterName (std::size_t index);

      bool compareAdditionalArguments (const Node & other) const;
      
      std::size_t hashAdditionalArguments () const;
      
      /// Configuration for numerical derivation of computation nodes using this Model.
      NumericalDerivativeConfiguration config;

      NodeRef recreate (Context & c, NodeRefVec && deps) final;
      
    private:
      void compute ();

      std::unique_ptr<TransitionModel> model_;
      
    };

    /** equilibriumFrequencies = f(model).
     * equilibriumFrequencies: RowVector(nbState).
     * model: ConfiguredModel.
     *
     * Node construction should be done with the create static method.
     */

    class EquilibriumFrequenciesFromModel : public Value<Eigen::RowVectorXd> {
    public:
      using Self = EquilibriumFrequenciesFromModel;
      using T = Eigen::RowVectorXd;

      static ValueRef<T> create (Context & c, NodeRefVec && deps, const Dimension<T> & dim);
      EquilibriumFrequenciesFromModel (NodeRefVec && deps, const Dimension<T> & dim);

      std::string debugInfo () const final;

      bool compareAdditionalArguments (const Node & other) const final;

      NodeRef derive (Context & c, const Node & node) final;
      NodeRef recreate (Context & c, NodeRefVec && deps) final;

    private:
      void compute () final;

      Dimension<T> targetDimension_;
    };

    /** transitionMatrix = f(model, branchLen).
     * transitionMatrix: Matrix(fromState, toState).
     * model: ConfiguredModel.
     * branchLen: double.
     *
     * Node construction should be done with the create static method.
     */

    class TransitionMatrixFromModel : public Value<Eigen::MatrixXd> {
    public:
      using Self = TransitionMatrixFromModel;
      using T = Eigen::MatrixXd;

      /// Build a new TransitionMatrixFromModel node with the given output dimensions.
      static ValueRef<T> create (Context & c, NodeRefVec && deps, const Dimension<T> & dim);
      TransitionMatrixFromModel (NodeRefVec && deps, const Dimension<T> & dim);

      std::string debugInfo () const final;

      bool compareAdditionalArguments (const Node & other) const final;

      NodeRef derive (Context & c, const Node & node) final;
      NodeRef recreate (Context & c, NodeRefVec && deps) final;

    private:
      void compute () final;

      Dimension<T> targetDimension_;
    };

    /** dtransitionMatrix/dbrlen = f(model, branchLen).
     * dtransitionMatrix/dbrlen: Matrix(fromState, toState).
     * model: ConfiguredModel.
     * branchLen: double.
     *
     * Node construction should be done with the create static method.
     */

    class TransitionMatrixFromModelFirstBrlenDerivative : public Value<Eigen::MatrixXd> {
    public:
      using Self = TransitionMatrixFromModelFirstBrlenDerivative;
      using T = Eigen::MatrixXd;

      static ValueRef<T> create (Context & c, NodeRefVec && deps, const Dimension<T> & dim);
      TransitionMatrixFromModelFirstBrlenDerivative (NodeRefVec && deps, const Dimension<T> & dim);

      std::string debugInfo () const final;

      bool compareAdditionalArguments (const Node & other) const final;

      NodeRef derive (Context & c, const Node & node) final;
      NodeRef recreate (Context & c, NodeRefVec && deps) final;

    private:
      void compute () final;

      Dimension<T> targetDimension_;
    };

    /** d2transitionMatrix/dbrlen2 = f(model, branchLen).
     * d2transitionMatrix/dbrlen2: Matrix(fromState, toState).
     * model: ConfiguredModel.
     * branchLen: double.
     *
     * Node construction should be done with the create static method.
     */

    class TransitionMatrixFromModelSecondBrlenDerivative : public Value<Eigen::MatrixXd> {
    public:
      using Self = TransitionMatrixFromModelSecondBrlenDerivative;
      using T = Eigen::MatrixXd;

      static ValueRef<T> create (Context & c, NodeRefVec && deps, const Dimension<T> & dim);
      TransitionMatrixFromModelSecondBrlenDerivative (NodeRefVec && deps, const Dimension<T> & dim);

      std::string debugInfo () const final;

      bool compareAdditionalArguments (const Node & other) const final;

      NodeRef derive (Context & c, const Node & node) final;
      NodeRef recreate (Context & c, NodeRefVec && deps) final;

    private:
      void compute () final;

      Dimension<T> targetDimension_;
    };
  } // namespace dataflow
} // namespace bpp

#endif // BPP_NEWPHYL_MODEL_H
