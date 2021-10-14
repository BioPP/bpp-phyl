//
// File: Parameter.h
// Authors:
//   Laurent Gueguen
// Created: mercredi 23 janvier 2019, ÃÂ  06h 02
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_PARAMETER_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_PARAMETER_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlow.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowNumeric.h>
#include <functional>
#include <unordered_map>


namespace bpp
{
/** @brief Data flow node representing a parameter.
 *
 * This class wraps a bpp::Parameter as a data flow node.
 *
 * It depends on a Value<double> node.
 */

class ConfiguredParameter : public Parameter,
  public Value<Parameter*>
{
private:
  const Context& context_;

public:
  using Self = ConfiguredParameter;
  using Target = Parameter;

  /// Build a new ConfiguredParameter node.
  static std::shared_ptr<Self> create (Context& c, NodeRefVec&& deps, const Parameter& param)
  {
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 1);
    checkNthDependencyIsValue<double>(typeid (Self), deps, 0);
    return cachedAs<Self>(c, std::make_shared<Self>(c, std::move(deps), param));
  }

  static std::shared_ptr<Self> create (Context& c, const Parameter& param)
  {
    NodeRefVec deps({NumericMutable<double>::create(c, param.getValue())});

    return cachedAs<Self>(c, std::make_shared<Self>(c, std::move(deps), param));
  }

  static std::shared_ptr<Self> resetDependencies(Context& c, std::shared_ptr<Self> self, NodeRefVec&& deps)
  {
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 1);
    checkNthDependencyIsValue<double>(typeid (Self), deps, 0);

    self->resetDependencies_(std::move(deps));

    return cachedAs<Self>(c, self);
  }


  ConfiguredParameter (const ConfiguredParameter& param)
    : Parameter (param), Value<Parameter*>({std::shared_ptr<NumericMutable<double> >(new NumericMutable<double>(param.getValue()))}, this), context_(param.context_)
  {}

  ConfiguredParameter (const Context& context, NodeRefVec&& deps, const Parameter& parameter);

  ConfiguredParameter (const Context& context, NodeRefVec&& deps, Parameter&& parameter);

  ConfiguredParameter* clone() const override
  {
    return new ConfiguredParameter(*this);
  }

  ~ConfiguredParameter ();

  std::string description () const override;
  std::string color () const override;
  std::string debugInfo () const override;

  /*
   * @brief setValue is transfered to the double dependency,
   * through the parameter constraints test.
   *
   */
  void setValue(double v) override
  {
    Parameter::setValue (v);                  // Will apply possible constraints
    static_cast<NumericMutable<double>&>(*dependency(0)).setValue(Parameter::getValue());
  }

  /*
   * @brief computation of double dependency is done and
   * transfered to parameter value, with constraint.
   *
   */
  double getValue() const override
  {
    double x = static_cast<NumericMutable<double>&>(*dependency(0)).getTargetValue();
    accessValueConst()->Parameter::setValue(x);
    return Parameter::getValue();
  }

  bool compareAdditionalArguments (const Node_DF& other) const override;

  std::size_t hashAdditionalArguments () const override;

  NodeRef derive (Context& c, const Node_DF& node) override;

  NodeRef recreate (Context& c, NodeRefVec&& deps) override;

  NodeRef recreate (Context& c);

private:
  void compute () override;
};


/** @brief shift param value = n * delta + x.
 * - x, delta: double.
 * - n: constant int.
 * - order of dependencies: (x, delta).
 *
 * Adds n * delta to value of x.
 * Used to generate x +/- delta values for numerical derivation.
 * Node construction should be done with the create static method.
 */

class ShiftParameter : public ConfiguredParameter
{
public:
  using Self = ShiftParameter;

  /// Build a new ShiftDelta node with the given output dimensions and shift number.
  static std::shared_ptr<ConfiguredParameter> create(Context& c, NodeRefVec&& deps, Parameter& param, const int n)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 2);
    checkDependencyRangeIsValue<double>(typeid (Self), deps, 0, 2);
    // Detect if we have a chain of ShiftParameter with the same delta.

    auto& delta = deps[1];
    auto* paramAsShiftParameter = dynamic_cast<ShiftParameter*>(&param);

    if (paramAsShiftParameter != nullptr && paramAsShiftParameter->dependency (1) == delta)
    {
      // Merge with ShiftParameter dependency by summing the n.
      Parameter p2(param);
      return Self::create (c, std::move(deps), p2, n + paramAsShiftParameter->getN ());
    }
    // Not a merge, select node implementation.
    if (n == 0 || delta->hasNumericalProperty (NumericalProperty::ConstantZero))
    {
      auto* paramAsConf = dynamic_cast<ConfiguredParameter*>(&param);
      if (paramAsConf)
        return std::dynamic_pointer_cast<ConfiguredParameter>(paramAsConf->shared_from_this());
      else
        return ConfiguredParameter::create(c, {deps[0]}, param);
    }
    else
    {
      return cachedAs<ConfiguredParameter>(c, std::make_shared<Self>(c, std::move(deps), param, n));
    }
  }

  ShiftParameter (const Context& context, NodeRefVec&& deps, const Parameter& parameter, int n)
    : ConfiguredParameter (context, std::move (deps), parameter), n_ (n)
  {}

  ShiftParameter (const Context& context, NodeRefVec&& deps, Parameter&& parameter, int n)
    : ConfiguredParameter (context, std::move (deps), std::move(parameter)), n_ (n)
  {}

  std::string description () const override { return "Shift" + ConfiguredParameter::description();}

  std::string color () const override
  {
    auto& name = getName();
    if (name.substr(0, 5) == "BrLen")
      return "#aa99aa";
    else
      return "#ffcc44";
  }


  std::string debugInfo () const override
  {
    return ConfiguredParameter::debugInfo() +
           " n=" + std::to_string (n_);
  }

  // ShiftDelta additional arguments = (n_).
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    const auto* derived = dynamic_cast<const Self*>(&other);
    return derived != nullptr && n_ == derived->n_;
  }

  /*
   * @brief setValue is not possible here, since computation is
   * done only through ompute method.
   *
   */
  void setValue(double v) override
  {
    throw Exception("ShiftParameter setValue should not be called");
  }


  /** @brief Raw value access (const).
   *
   * Value is not guaranteed to be valid (no recomputation).
   */
  double getValue() const override
  {
    return Parameter::getValue();
  }

  std::size_t hashAdditionalArguments () const final
  {
    std::size_t seed = 0;
    combineHash (seed, n_);
    return seed;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<Parameter>::create (c, Dimension<Parameter>());
    }
    auto dx = this->dependency (0)->derive (c, node);
    auto ddelta = this->dependency (1)->derive (c, node);
    if (dx->hasNumericalProperty (NumericalProperty::ConstantZero))
      return ConstantZero<Parameter>::create (c, Dimension<Parameter>());

    if (dx->hasNumericalProperty (NumericalProperty::ConstantOne))
      return ConstantOne<Parameter>::create (c, Dimension<Parameter>());

    return Self::create (c, {dx, ddelta}, *this, n_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), *this, n_);
  }

  int getN () const { return n_; }

private:
  void compute () final
  {
    const auto&  x = accessValueConstCast<double>(*this->dependency (0));
    const auto& delta = accessValueConstCast<double>(*this->dependency (1));
    double r = n_ * delta + x;
    // Boundary mgmt not so clean!
    this->accessValueMutable()->Parameter::setValue(hasConstraint() ? (getConstraint()->isCorrect(r) ? r : getConstraint()->getAcceptedLimit(r)) : r);
  }

  int n_;
};

/** Value = parameter.getValue().
 * parameter: ConfiguredParameter.
 *
 * Node construction should be done with the create static method.
 */


class ValueFromConfiguredParameter : public Value<double>
{
public:
  using Self = ValueFromConfiguredParameter;
  using Dep = ConfiguredParameter;
  using T = double;

  ValueFromConfiguredParameter (NodeRefVec&& deps);

  std::string debugInfo () const final;

  bool compareAdditionalArguments (const Node_DF& other) const final;

  NodeRef derive (Context& c, const Node_DF& node) final;
  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

  std::string description () const
  {
    return "ValueFrom(" + dependency(0)->description() + ")";
  }

  std::string color () const;

private:
  void compute () final;

public:
  static std::shared_ptr<Self> create (Context& c, NodeRefVec&& deps);
};
} // namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_PARAMETER_H
