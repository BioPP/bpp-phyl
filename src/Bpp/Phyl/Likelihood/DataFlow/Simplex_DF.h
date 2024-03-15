// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_SIMPLEX_DF_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_SIMPLEX_DF_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Prob/Simplex.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlow.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowCWiseComputing.h>
#include <functional>
#include <unordered_map>

#include "Definitions.h"

namespace bpp
{
class Simplex;

/** @brief Data flow node representing a Frequencies Set
 * configured with parameter values.
 *
 * This class wraps a bpp::Simplex as a data flow node.
 *
 * It depends on Value<double> nodes (one for each parameter
 * declared in the freq set).
 * It provides a dummy value representing the "frequencies set
 * configured by its parameters".
 *
 * The dummy value is implemented as a pointer to the internal
 * frequencies set for simplicity.
 */

class ConfiguredSimplex : public Value<const Simplex*>,
  public AbstractParametrizable
{
  // private:

  //   const Context& context_;

public:
  using Self = ConfiguredSimplex;
  using Target = Simplex;

  ConfiguredSimplex (const Context& context, NodeRefVec&& deps, std::unique_ptr<Simplex>&& simplex);
  ~ConfiguredSimplex ();

  ConfiguredSimplex* clone() const
  {
    throw bpp::Exception("ConfiguredSimplex clone should not happen.");
  }

  std::string description () const final;
  std::string debugInfo () const final;
  std::string color() const final
  {
    return "blue";
  }

  bool compareAdditionalArguments (const Node_DF& other) const;

  std::size_t hashAdditionalArguments () const;

  /// Configuration for numerical derivation of computation nodes using this Simplex.
  NumericalDerivativeConfiguration config;

  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

  const ConfiguredParameter& getConfiguredParameter(const std::string& name)
  {
    return static_cast<const ConfiguredParameter&>(parameter(name));
  }

private:
  void compute ()
  {
    simplex_->matchParametersValues(getParameters());
  }


  std::unique_ptr<Simplex> simplex_;
};

/** Frequencies = f(Simplex).
 * Frequencies: RowVector(nbState).
 * Frequenciesset: ConfiguredSimplex.
 *
 * Node construction should be done with the create static method.
 */

class FrequenciesFromSimplex : public Value<Eigen::RowVectorXd>
{
public:
  using Self = FrequenciesFromSimplex;
  using T = Eigen::RowVectorXd;

  // static ValueRef<T> create (Context & c, NodeRefVec && deps, const Dimension<T> & dim);
  FrequenciesFromSimplex (NodeRefVec&& deps, const Dimension<T>& dim);

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

  Dimension<T> targetDimension_;

public:
  static std::shared_ptr<Self> create (Context& c, NodeRefVec&& deps);
};
} // namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_SIMPLEX_DF_H
