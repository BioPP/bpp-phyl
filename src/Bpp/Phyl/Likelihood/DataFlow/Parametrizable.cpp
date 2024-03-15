// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/ParameterAliasable.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Parameter.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Parametrizable.h>


using namespace std;

namespace bpp
{
std::unordered_map<std::string, std::shared_ptr<ConfiguredParameter> >
createParameterMap (Context& c, const ParameterAliasable& parametrizable)
{
  const auto& parameters = parametrizable.getIndependentParameters ();
  const auto nbParameters = parameters.size ();
  std::unordered_map<std::string, std::shared_ptr<ConfiguredParameter> > map;
  for (std::size_t i = 0; i < nbParameters; ++i)
  {
    const auto& param = parameters[i];
    auto value = NumericMutable<double>::create (c, param.getValue ());
    map.emplace (param.getName (),
                 ConfiguredParameter::create (c, {std::move(value)}, param));
  }
  return map;
}

std::unordered_map<std::string, std::shared_ptr<ConfiguredParameter> >
createParameterMap (Context& c, const Parametrizable& parametrizable)
{
  const auto& parameters = parametrizable.getParameters ();
  const auto nbParameters = parameters.size ();
  std::unordered_map<std::string, std::shared_ptr<ConfiguredParameter> > map;
  for (std::size_t i = 0; i < nbParameters; ++i)
  {
    const auto& param = parameters[i];
    auto value = NumericMutable<double>::create (c, param.getValue ());
    map.emplace (param.getName (),
                 ConfiguredParameter::create (c, {std::move(value)}, param));
  }
  return map;
}

NodeRefVec createDependencyVector (const Parametrizable& parametrizable,
                                   const std::function<NodeRef (const std::string&)>& parameter)
{
  const auto& parameters = parametrizable.getParameters ();
  const auto nbParameters = parameters.size ();
  NodeRefVec deps (nbParameters);
  for (std::size_t i = 0; i < nbParameters; ++i)
  {
    auto dep = parameter(parameters[i].getName ());
    if (!dep)
    {
      throw Exception ("createDependencyVector (Parametrizable): parameter not found: " + parameters[i].getName ());
    }
    deps[i] = std::move (dep);
  }
  return deps;
}

NodeRefVec createDependencyVector (const ParameterAliasable& parametrizable,
                                   const std::function<NodeRef (const std::string&)>& parameter)
{
  const auto& parameters = parametrizable.getIndependentParameters ();
  const auto nbParameters = parameters.size ();
  NodeRefVec deps (nbParameters);
  for (std::size_t i = 0; i < nbParameters; ++i)
  {
    auto dep = parameter(parameters[i].getName ());
    if (!dep)
    {
      throw Exception ("createDependencyVector (Parametrizable): parameter not found: " + parameters[i].getName ());
    }
    deps[i] = std::move (dep);
  }
  return deps;
}
} // namespace bpp
