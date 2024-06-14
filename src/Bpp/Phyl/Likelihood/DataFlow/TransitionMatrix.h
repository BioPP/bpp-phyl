// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_TRANSITIONMATRIX_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_TRANSITIONMATRIX_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Hmm/HmmTransitionMatrix.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowCWiseComputing.h>
#include <functional>
#include <unordered_map>

#include "Definitions.h"

namespace bpp
{
/** @brief Data flow node representing a TransitionMatrix
 * configured with parameter values.
 *
 * It depends on Value<double> nodes (one for each parameter).
 * It provides a dummy value representing the "transitionmatrix configured by its parameters".
 */

class ConfiguredTransitionMatrix : public Value<const HmmTransitionMatrix*>,
  public AbstractParametrizable
{
  // private:
  //   Context& context_;

public:
  using Self = ConfiguredTransitionMatrix;
  using Target = HmmTransitionMatrix;

  ConfiguredTransitionMatrix (Context& context, NodeRefVec&& deps, std::unique_ptr<HmmTransitionMatrix>&& hmm);
  ~ConfiguredTransitionMatrix ();

  ConfiguredTransitionMatrix* clone() const
  {
    throw bpp::Exception("ConfiguredTransitionMatrix clone should not happen.");
  }

  std::string description () const final;
  std::string debugInfo () const final;

  std::string color () const final
  {
    return "green";
  }

  bool compareAdditionalArguments (const Node_DF& other) const;

  std::size_t hashAdditionalArguments () const;

  /// Configuration for numerical derivation of computation nodes using this TransitionMatrix.
  NumericalDerivativeConfiguration config;

  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

  const ConfiguredParameter& getConfiguredParameter(const std::string& name) const
  {
    return static_cast<const ConfiguredParameter&>(parameter(name));
  }

private:
  void compute ()
  {
    hmm_->matchParametersValues(getParameters());
  }

  std::unique_ptr<HmmTransitionMatrix> hmm_;
};

/** equilibriumFrequencies = f(model).
 * equilibriumFrequencies: Vector(nbState).
 * model: ConfiguredTransitionMatrix.
 *
 * Node construction should be done with the create static method.
 */

class EquilibriumFrequenciesFromTransitionMatrix : public Value<Eigen::VectorXd>
{
public:
  using Self = EquilibriumFrequenciesFromTransitionMatrix;
  using Dep = ConfiguredTransitionMatrix;
  using T = Eigen::VectorXd;

  EquilibriumFrequenciesFromTransitionMatrix (NodeRefVec&& deps, const Dimension<T>& dim);

  std::string debugInfo () const final;

  bool compareAdditionalArguments (const Node_DF& other) const;

  NodeRef derive (Context& c, const Node_DF& node) final;
  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

  std::string color () const final
  {
    return "#ffff66";
  }

private:
  void compute () final;

  Dimension<T> targetDimension_;
};

/** transitionMatrix = f(HmmTransitionMatrix).
 * transitionMatrix: Matrix(fromState, toState).
 *
 * Node construction should be done with the create static method.
 */

class TransitionMatrixFromTransitionMatrix : public Value<Eigen::MatrixXd>
{
public:
  using Self = TransitionMatrixFromTransitionMatrix;
  using Dep = ConfiguredTransitionMatrix;
  using T = Eigen::MatrixXd;

private:
  Dimension<T> targetDimension_;

public:
/// Build a new TransitionMatrixFromTransitionMatrix node with the given output dimensions.
  TransitionMatrixFromTransitionMatrix (NodeRefVec&& deps, const Dimension<T>& dim);

  std::string debugInfo () const final;

  bool compareAdditionalArguments (const Node_DF& other) const;

  NodeRef derive (Context& c, const Node_DF& node) final;
  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

  std::string color () const final
  {
    return "#aaff00";
  }

  std::string description () const final
  {
    return "TransitionMatrix";
  }

  std::string shape() const
  {
    return "octagon";
  }

private:
  void compute () final;
};
} // namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_TRANSITIONMATRIX_H
