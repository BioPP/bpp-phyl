// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_DISCRETEDISTRIBUTION_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_DISCRETEDISTRIBUTION_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlow.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Parametrizable.h>
#include <functional>
#include <unordered_map>

#include "Definitions.h"

namespace bpp
{
/* Likelihood discrete distribution.
 */

/** @brief Data flow node representing a DiscreteDistribution
 * configured with parameter values.
 *
 * This class wraps a bpp::DiscreteDistribution as a data flow node.
 * It depends on Value<double> nodes (one for each parameter
 * declared in the distribution).
 * It provides a dummy value representing the "distribution
 * configured by its parameters".
 *
 * The dummy value is implemented as a pointer to the internal
 * distribution for simplicity.
 */

class ConfiguredDistribution : public Value<const DiscreteDistributionInterface*>,
  public AbstractParametrizable
{
  // private:
  //   Context& context_;

public:
  using Self = ConfiguredDistribution;
  using Target = DiscreteDistributionInterface;

  ConfiguredDistribution (Context& context, NodeRefVec&& deps, std::unique_ptr<DiscreteDistributionInterface>&& distrib);
  ~ConfiguredDistribution ();

  ConfiguredDistribution* clone() const
  {
    throw bpp::Exception("ConfiguredDistribution clone should not happen.");
  }

  std::string description () const final;
  std::string debugInfo () const final;
  std::string color() const final
  {
    return "blue";
  }

  bool compareAdditionalArguments (const Node_DF& other) const;

  std::size_t hashAdditionalArguments () const;

  /// Configuration for numerical derivation of computation nodes using this Model.
  NumericalDerivativeConfiguration config;

  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

  const ConfiguredParameter& getConfiguredParameter(const std::string& name)
  {
    return static_cast<const ConfiguredParameter&>(parameter(name));
  }

private:
  void compute ()
  {
    distrib_->matchParametersValues(getParameters());
  }

  std::unique_ptr<DiscreteDistributionInterface> distrib_;
};

/** Probabilities = f(DiscreteDistribution).
 * Probabilities: RowVector(nbClass).
 * DiscreteDistribution: ConfiguredDistribution.
 *
 * Node construction should be done with the create static method.
 */

class ProbabilitiesFromDiscreteDistribution : public Value<Eigen::RowVectorXd>
{
public:
  using Self = ProbabilitiesFromDiscreteDistribution;
  using Dep = ConfiguredDistribution;
  using T = Eigen::RowVectorXd;

  ProbabilitiesFromDiscreteDistribution (NodeRefVec&& deps, const Dimension<T>& dim);

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

/** Proba = f(DiscreteDistribution, Category).
 * Proba: Double.
 * DiscreteDistribution: ConfiguredDistribution.
 * Category: number of the category
 *
 * Node construction should be done with the create static method.
 */

class ProbabilityFromDiscreteDistribution : public Value<double>
{
private:
  unsigned int nCat_;

public:
  using Self = ProbabilityFromDiscreteDistribution;
  using Dep = ConfiguredDistribution;
  using T = double;

  ProbabilityFromDiscreteDistribution (NodeRefVec&& deps, unsigned int nCat_);

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
  static std::shared_ptr<Self> create (Context& c, NodeRefVec&& deps, unsigned int nCat);
};

/** Rate = f(DiscreteDistribution, Category).
 * Rate: Double.
 * DiscreteDistribution: ConfiguredDistribution.
 * Category: number of the category
 *
 * Node construction should be done with the create static method.
 */

class CategoryFromDiscreteDistribution : public Value<double>
{
private:
  unsigned int nCat_;

public:
  using Self = CategoryFromDiscreteDistribution;
  using Dep = ConfiguredDistribution;
  using T = double;

  CategoryFromDiscreteDistribution (NodeRefVec&& deps, unsigned int nCat_);

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
  static std::shared_ptr<Self> create (Context& c, NodeRefVec&& deps, unsigned int nCat);
};
} // namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_DISCRETEDISTRIBUTION_H
