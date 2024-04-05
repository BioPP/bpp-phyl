// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Parameter.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Parametrizable.h>


using namespace std;

namespace bpp
{
// Parameter node

ConfiguredParameter::ConfiguredParameter (const Context& context, NodeRefVec&& deps, const Parameter& parameter)
  : Parameter(parameter), Value<Parameter*>(deps, this), context_(context)
{}

ConfiguredParameter::ConfiguredParameter (const Context& context, NodeRefVec&& deps, Parameter&& parameter)
  : Parameter(std::move(parameter)), Value<Parameter*>(deps, this), context_(context)
{}

ConfiguredParameter::~ConfiguredParameter () = default;

std::string ConfiguredParameter::description () const
{
  return "Parameter(" + getName () + ")\nvalue=" + std::to_string(getValue());
}

std::string ConfiguredParameter::debugInfo () const
{
  return "";
}

std::string ConfiguredParameter::color () const
{
  auto& name = getName();
  if (name.substr(0, 5) == "BrLen")
    return "#00ff00";
  else
    return "#ff8800";
}

// Parameter node additional arguments = (type of bpp::Parameter).
// Everything else is determined by the node dependencies.

bool ConfiguredParameter::compareAdditionalArguments (const Node_DF& other) const
{
  const auto* derived = dynamic_cast<const Self*>(&other);
  if (derived == nullptr)
  {
    return false;
  }
  else
  {
    return this->getName() == derived->getName();
  }
}

std::size_t ConfiguredParameter::hashAdditionalArguments () const
{
  const auto& bppFS = *this;
  return typeid (bppFS).hash_code ();
}

NodeRef ConfiguredParameter::derive (Context& c, const Node_DF& node)
{
  if (&node == this)
  {
    return ConstantOne<double>::create (c, Dimension<double>());
  }
  else
    return this->dependency(0)->derive(c, node);
}

NodeRef ConfiguredParameter::recreate (Context& c)
{
  return ConfiguredParameter::create (c, NodeRefVec{NumericMutable<double>::create(c, accessValueConstCast<double>(*dependency(0)))}, *this);
}

NodeRef ConfiguredParameter::recreate (Context& c, NodeRefVec&& deps)
{
  return ConfiguredParameter::create (c, std::move (deps), *this);
}

void ConfiguredParameter::compute ()
{
  // Update with the numerical dependency
  auto& v = accessValueConstCast<double>(*this->dependency (0));
  if (Parameter::getValue () != v)
  {
    Parameter::setValue(v);
  }
}

////////////////////////////////////////////////////
// ValueFromConfiguredParameter

ValueFromConfiguredParameter::ValueFromConfiguredParameter (
    NodeRefVec&& deps)
  : Value<double>(std::move (deps)) {}

std::string ValueFromConfiguredParameter::debugInfo () const
{
  return dependency(0)->debugInfo();
}

// ValueFromConfiguredParameter additional arguments = ().
bool ValueFromConfiguredParameter::compareAdditionalArguments (const Node_DF& other) const
{
  const auto* derived = dynamic_cast<const Self*>(&other);
  return derived != nullptr;
}

std::shared_ptr<ValueFromConfiguredParameter> ValueFromConfiguredParameter::create (Context& c, NodeRefVec&& deps)
{
  checkDependenciesNotNull (typeid (Self), deps);
  checkDependencyVectorSize (typeid (Self), deps, 1);
  checkNthDependencyIs<ConfiguredParameter>(typeid (Self), deps, 0);
  return cachedAs<ValueFromConfiguredParameter>(c, std::make_shared<ValueFromConfiguredParameter>(std::move (deps)));
}

NodeRef ValueFromConfiguredParameter::derive (Context& c, const Node_DF& node)
{
  if (&node == this)
  {
    return ConstantOne<double>::create (c, Dimension<double>());
  }
  else
    return this->dependency(0)->derive(c, node);
}

std::string ValueFromConfiguredParameter::color () const
{
  auto& name = accessValueConstCast<const ConfiguredParameter*>(*this->dependency (0))->getName();
  if (name.substr(0, 5) == "BrLen")
    return "#66ff66";
  else
    return "#ffcc00";
}

NodeRef ValueFromConfiguredParameter::recreate (Context& c, NodeRefVec&& deps)
{
  return ValueFromConfiguredParameter::create (c, {std::move (deps)});
}

void ValueFromConfiguredParameter::compute ()
{
  const auto* param = accessValueConstCast<const ConfiguredParameter*>(*this->dependency (0));
  this->accessValueMutable () = param->getValue();
}
} // namespace bpp
