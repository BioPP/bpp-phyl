// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_MODEL_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_MODEL_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowCWiseComputing.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <functional>
#include <unordered_map>

#include "Definitions.h"

namespace bpp
{
inline MatrixDimension transitionMatrixDimension (std::size_t nbState)
{
  return {Eigen::Index (nbState), Eigen::Index (nbState)};
}

inline RowVectorDimension equilibriumFrequenciesDimension (std::size_t nbState)
{
  return RowVectorDimension (Eigen::Index (nbState));
}


/**
 * @brief Likelihood transition model.
 */
//class TransitionModelInterface;

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
class ConfiguredModel :
  public Value<std::shared_ptr<BranchModelInterface>>,
  public AbstractParametrizable
{
  // private:
  //   Context& context_;

public:
  using Self = ConfiguredModel;
  using Target = BranchModelInterface;

  ConfiguredModel(Context& context, NodeRefVec&& deps, std::shared_ptr<BranchModelInterface>&& model);
  virtual ~ConfiguredModel();

  ConfiguredModel* clone() const override
  {
    throw bpp::Exception("ConfiguredModel clone should not happen.");
  }

  std::string description() const final;
  std::string debugInfo() const final;

  std::string color() const final
  {
    return "red";
  }

  bool compareAdditionalArguments(const Node_DF& other) const override;

  std::size_t hashAdditionalArguments () const override;

  /// Configuration for numerical derivation of computation nodes using this Model.
  NumericalDerivativeConfiguration config;

  NodeRef recreate(Context& c, NodeRefVec&& deps) final;

  const ConfiguredParameter& getConfiguredParameter(const std::string& name) const
  {
    return static_cast<const ConfiguredParameter&>(parameter(name));
  }

private:
  void compute() override
  {
    model_->matchParametersValues(getParameters());
  }

  std::shared_ptr<BranchModelInterface> model_;
};

/** 
 * equilibriumFrequencies = f(model).
 * equilibriumFrequencies: RowVector(nbState).
 * model: ConfiguredModel.
 *
 * Node construction should be done with the create static method.
 */

class EquilibriumFrequenciesFromModel : public Value<Eigen::RowVectorXd>
{
public:
  using Self = EquilibriumFrequenciesFromModel;
  using Dep = ConfiguredModel;
  using T = Eigen::RowVectorXd;

  EquilibriumFrequenciesFromModel (NodeRefVec&& deps, const Dimension<T>& dim);

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

/** transitionMatrix = f(model, branchLen, nDeriv, nMod, mult).
 * transitionMatrix: Matrix(fromState, toState).
 * model: ConfiguredModel.
 * branchLen: double.
 * nDeriv: degree of derivate (default: 0)
 * nMod: in case of mixture model, takes the number of submodel
 *       where the generator comes from (optional, or null).
 *
 * Node construction should be done with the create static method.
 */

class TransitionMatrixFromModel : public Value<Eigen::MatrixXd>
{
public:
  using Self = TransitionMatrixFromModel;
  using Dep = ConfiguredModel;
  using T = Eigen::MatrixXd;

private:
  Dimension<T> targetDimension_;

public:
/// Build a new TransitionMatrixFromModel node with the given output dimensions.
  TransitionMatrixFromModel (NodeRefVec&& deps, const Dimension<T>& dim);

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

/** transitionProbability = f(model, branchLen, nDeriv).
 * transitionProbability: f(fromState, vector) -> probability
 *
 * model: ConfiguredModel.
 * branchLen: double.
 * nDeriv: derivation level
 *
 * Node construction should be done with the create static method.
 */

class TransitionFunctionFromModel : public Value<bpp::TransitionFunction>
{
public:
  using Self = TransitionFunctionFromModel;
  using Dep = ConfiguredModel;
  using T = bpp::TransitionFunction;

private:
  Dimension<T> targetDimension_;

public:
  /// Build a new TransitionMatrixFromModel node with the given output dimensions.
  TransitionFunctionFromModel(NodeRefVec&& deps, const Dimension<T>& dim);

  std::string debugInfo() const final;

  bool compareAdditionalArguments(const Node_DF& other) const;

  NodeRef derive(Context& c, const Node_DF& node) final;
  NodeRef recreate(Context& c, NodeRefVec&& deps) final;

  std::string color () const final
  {
    return "#aaff00";
  }

  std::string description () const final
  {
    return "TransitionFunction";
  }

  std::string shape() const
  {
    return "octagon";
  }

private:
  void compute() final;

public:
  static std::shared_ptr<Self> create(Context& c, NodeRefVec&& deps, const Dimension<T>& dim);
};

/** 
 * dtransitionMatrix/dbrlen = f(model, branchLen).
 * dtransitionMatrix/dbrlen: Matrix(fromState, toState).
 * model: ConfiguredModel.
 * branchLen: double.
 *
 * Node construction should be done with the create static method.
 */
class ProbabilitiesFromMixedModel : public Value<Eigen::RowVectorXd>
{
public:
  using Self = ProbabilitiesFromMixedModel;
  using Dep = ConfiguredModel;
  using T = Eigen::RowVectorXd;

  ProbabilitiesFromMixedModel (NodeRefVec&& deps, const Dimension<T>& dim);

  std::string debugInfo () const final;

  std::string color() const final
  {
    return "blue";
  }

  bool compareAdditionalArguments (const Node_DF& other) const final;

  NodeRef derive (Context& c, const Node_DF& node) final;
  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

private:
  void compute () final;

  Dimension<T> nbClass_;

public:
  static std::shared_ptr<Self> create (Context& c, NodeRefVec&& deps);
};

/** Proba = f(MixedModel, Category).
 * Proba: Double.
 * MixedModel: ConfiguredModel
 * Category: number of the category
 *
 * Node construction should be done with the create static method.
 */

class ProbabilityFromMixedModel : public Value<double>
{
private:
  size_t nCat_;

public:
  using Self = ProbabilityFromMixedModel;
  using Dep = ConfiguredModel;
  using T = double;

  ProbabilityFromMixedModel (NodeRefVec&& deps, size_t nCat);

  std::string debugInfo () const final;

  bool compareAdditionalArguments (const Node_DF& other) const final;

  NodeRef derive (Context& c, const Node_DF& node) final;
  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

  std::string color() const final
  {
    return "blue";
  }

private:
  void compute () final;

public:
  static std::shared_ptr<Self> create (Context& c, NodeRefVec&& deps, size_t nCat);
};
} // namespace bpp
#endif// BPP_PHYL_LIKELIHOOD_DATAFLOW_MODEL_H
