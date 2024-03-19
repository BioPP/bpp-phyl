// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_FREQUENCYSET_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_FREQUENCYSET_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlow.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowCWiseComputing.h>
#include <Bpp/Phyl/Model/FrequencySet/FrequencySet.h>
#include <functional>
#include <unordered_map>

#include "Definitions.h"

namespace bpp
{
class FrequencySetInterface;

/**
 * @brief Data flow node representing a Frequencies Set
 * configured with parameter values.
 *
 * This class wraps a bpp::FrequencySet as a data flow node.
 *
 * It depends on Value<double> nodes (one for each parameter
 * declared in the freq set).
 * It provides a dummy value representing the "frequencies set
 * configured by its parameters".
 *
 * The dummy value is implemented as a pointer to the internal
 * frequencies set for simplicity.
 */
class ConfiguredFrequencySet :
  public Value<const FrequencySetInterface*>,
  public AbstractParametrizable
{
  // private:

  //   const Context& context_;

public:
  using Self = ConfiguredFrequencySet;
  using Target = FrequencySetInterface;

  ConfiguredFrequencySet (const Context& context, NodeRefVec&& deps, std::unique_ptr<FrequencySetInterface>&& freqset);

  virtual ~ConfiguredFrequencySet ();

  ConfiguredFrequencySet* clone() const
  {
    throw bpp::Exception("ConfiguredFrequencySet clone should not happen.");
  }

  std::string description () const final;
  std::string debugInfo () const final;
  std::string color() const final
  {
    return "#ffff00";
  }

  bool compareAdditionalArguments (const Node_DF& other) const;

  std::size_t hashAdditionalArguments () const;

  /// Configuration for numerical derivation of computation nodes using this FrequencySet.
  NumericalDerivativeConfiguration config;

  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

  const ConfiguredParameter& getConfiguredParameter(const std::string& name)
  {
    return static_cast<const ConfiguredParameter&>(parameter(name));
  }

private:
  void compute ()
  {
    freqset_->matchParametersValues(getParameters());
  }


  std::unique_ptr<FrequencySetInterface> freqset_;
};

/** Frequencies = f(FrequencySet).
 * Frequencies: RowVector(nbState).
 * Frequenciesset: ConfiguredFrequencySet.
 *
 * Node construction should be done with the create static method.
 */

class FrequenciesFromFrequencySet : public Value<Eigen::RowVectorXd>
{
public:
  using Self = FrequenciesFromFrequencySet;
  using T = Eigen::RowVectorXd;

  // static ValueRef<T> create (Context & c, NodeRefVec && deps, const Dimension<T> & dim);
  FrequenciesFromFrequencySet (NodeRefVec&& deps, const Dimension<T>& dim);

  std::string debugInfo () const final;

  bool compareAdditionalArguments (const Node_DF& other) const final;

  NodeRef derive (Context& c, const Node_DF& node) final;
  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

  std::string color() const final
  {
    return "#ffff66";
  }

private:
  void compute () final;

  Dimension<T> targetDimension_;
};
} // namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_FREQUENCYSET_H
