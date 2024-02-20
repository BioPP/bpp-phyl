//
// File: Parametrizable.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-03 00:00:00
// Last modified: 2017-05-03 00:00:00
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

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_PARAMETRIZABLE_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_PARAMETRIZABLE_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowCWiseComputing.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Parameter.h>
#include <functional>
#include <iostream>
#include <unordered_map>

#include "Definitions.h"

namespace bpp
{
/** 
 * Helper: create a map with mutable dataflow nodes for each
 *   parameter of the parametrizable.
 * The map is indexed by parameter names.
 */
std::unordered_map<std::string, std::shared_ptr<ConfiguredParameter>>
createParameterMap(Context& c, const Parametrizable& parametrizable);

std::unordered_map<std::string, std::shared_ptr<ConfiguredParameter>>
createParameterMap(Context& c, const ParameterAliasable& parametrizable);


/**
 * @brief Create a dependency vector suitable for a parametrizable class constructor.
 * The vector is built from parameter names, and an opaque accessor function.
 * For each named parameter, parameter(name) should return a valid node.
 * If no node is found (NodeRef was null), an exception is thrown.
 * Returned nodes must be Value<double> nodes.
 */
NodeRefVec createDependencyVector(const Parametrizable& parametrizable,
                                  const std::function<NodeRef(const std::string&)>& getParameter);

NodeRefVec createDependencyVector(const ParameterAliasable& parametrizable,
                                  const std::function<NodeRef(const std::string&)>& getParameter);


class ConfiguredParametrizable
{
public:
  /**
   * @brief Create a new node to an Aliased Parameterized object from
   * a dependency vector of Independent Parameters (which are already
   * Configure).
   *
   * Object parameters are given by a dependency vector of Value<double> nodes.
   * The number and order of parameters is given by the Object
   * internal ParameterList.
   */
  
  template<typename Object, typename Self>
  static std::shared_ptr<Self> createConfigured (Context& c, NodeRefVec&& deps,
                                                 std::unique_ptr<Object>&& object,
                                                 typename std::enable_if<std::is_base_of<ParameterAliasable, Object>::value>::type* = 0)
  {
    if (!object)
    {
      throw Exception ("createConfigured(): nullptr object");
    }
    // Check dependencies
    const auto nbParameters = object->getIndependentParameters ().size ();
    checkDependenciesNotNull(typeid (Self), deps);
    checkDependencyVectorSize(typeid (Self), deps, nbParameters);
    checkDependencyRangeIsValue<Parameter*>(typeid (Self), deps, 0, nbParameters);
    return cachedAs<Self>(c, std::make_shared<Self>(c, std::move (deps), std::move (object)));
  }


  /**
   * @brief Create a new node to a Parameterized object from a
   * dependency vector of Parameters (which are already Configured).
   *
   * Object parameters are given by a dependency vector of Value<double> nodes.
   * The number and order of parameters is given by the Object
   * internal ParameterList.
   */
  
  template<typename Object, typename Self>
  static std::shared_ptr<Self> createConfigured (Context& c, NodeRefVec&& deps,
                                                 std::unique_ptr<Object>&& object,
                                                 typename std::enable_if<!std::is_base_of<ParameterAliasable, Object>::value>::type* = 0)
  {
    if (!object)
    {
      throw Exception ("createConfigured(): nullptr object");
    }
    // Check dependencies
    const auto nbParameters = object->getParameters ().size ();
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, nbParameters);
    checkDependencyRangeIsValue<Parameter*>(typeid (Self), deps, 0, nbParameters);
    return cachedAs<Self>(c, std::make_shared<Self>(c, std::move (deps), std::move (object)));
  }

  
  /**
   * @brief Create a new node to a Parameterized object
   * ConfiguredParameters are built from its set of independent Parameters.
   *
   */
  
  template<typename Object, typename Self>
  static std::shared_ptr<Self> createConfigured (Context& context, std::unique_ptr<Object>&& object)
  {
    auto objectParameters = createParameterMap(context, *object);

    auto depvecObject = createDependencyVector(
      *object, [&objectParameters](const std::string& paramName) { return objectParameters[paramName]; });

    auto objectNode = ConfiguredParametrizable::createConfigured<Object, Self>(
      context,
      std::move(depvecObject),
      std::move(object));


    return objectNode;
  }

  /**
   * @brief Create a new node from a clonable object. Object
   * parameters are get from ConfiguredParameters PREVIOUSLY built
   * and stored in a ParameterList.
   *
   * @param context
   * @param object : Object which will be cloned to the node.
   * @param parList: list of the ConfiguredParameters previously built
   * @param suff  optional suffix for the parameter name, in the context of a full DataFlow
   */
  template<typename Object, typename Self>
  static std::shared_ptr<Self> createConfigured(Context& context, const Object& object, ParameterList& parList, const std::string& suff = "")
  {
    auto nObject = std::unique_ptr<Object>(dynamic_cast<Object*>(object.clone()));

    auto pa = dynamic_cast<const ParameterAliasable*>(&object);

    const ParameterList& lParams = pa ? pa->getIndependentParameters() : object.getParameters();

    std::vector<NodeRef> dep;
    for (size_t i = 0; i < lParams.size(); i++)
    {
      std::string name = lParams[i].getName() + suff;
      if (!parList.hasParameter(name) && suff == "")
      {
        if (!parList.hasParameter(lParams[i].getName() + "_1"))
          throw Exception("createConfigured: unknow ConfiguredParameter " + name);
        else name = lParams[i].getName() + "_1";
      }
      auto confPar = dynamic_cast<ConfiguredParameter*>(parList.getParameter(name).get());
      if (!confPar)
        throw Exception("createConfigured: unknown ConfiguredParameter " + name);

      dep.push_back(ConfiguredParameter::create(context, {confPar->dependency(0)}, lParams[i]));
    }
    return ConfiguredParametrizable::createConfigured<Object, Self>(context, std::move(dep), std::move(nObject));
  }


  /** Create a new vector of dependencies node to a double,
   * through a inheriting class (aka Self), from a
   * ConfiguredObject
   *
   */

  template<typename ConfiguredObject, typename Self>
  static ValueRef<double>
  createDouble (Context& c, NodeRefVec&& deps)
  {
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 1);
    checkNthDependencyIs<ConfiguredObject>(typeid (Self), deps, 0);
    return cachedAs<Value<double> >(c, std::make_shared<Self>(std::move (deps)));
  }

  /** Create a new vector of dependencies node to a RowVectorXd,
   * through a inheriting class (aka Self), from a
   * ConfiguredObject
   *
   * Additional dependencies are allowed
   *
   */

  template<typename ConfiguredObject, typename Self, typename Row>
  static ValueRef<Row>
  createRowVector (Context& c, NodeRefVec&& deps,
                   const Dimension<Row>& dim)
  {
    checkDependencyVectorMinSize (typeid (Self), deps, 1);
    checkNthDependencyNotNull (typeid (Self), deps, 0);
    checkNthDependencyIs<ConfiguredObject>(typeid (Self), deps, 0);
    return cachedAs<Value<Row> >(c, std::make_shared<Self>(std::move (deps), dim));
  }

  template<typename ConfiguredObject, typename Self, typename Col>
  static ValueRef<Col>
  createVector (Context& c, NodeRefVec&& deps,
                const Dimension<Col>& dim)
  {
    checkDependencyVectorMinSize (typeid (Self), deps, 1);
    checkNthDependencyNotNull (typeid (Self), deps, 0);
    checkNthDependencyIs<ConfiguredObject>(typeid (Self), deps, 0);
    return cachedAs<Value<Col> >(c, std::make_shared<Self>(std::move (deps), dim));
  }

  /** Create a new vector of dependencies node to a MatrixXd,
   * through a inheriting class (aka Self), from a
   * ConfiguredObject and an optional ConfiguredParameter (ie branch length)
   *
   */

  template<typename ConfiguredObject, typename Self, typename Matrix>
  static ValueRef<Matrix>
  createMatrix (Context& c, NodeRefVec&& deps,
                const Dimension<Matrix>& dim)
  {
    checkDependencyVectorMinSize (typeid (Self), deps, 1);
    checkNthDependencyNotNull (typeid (Self), deps, 0);
    if (deps.size() > 1)
      checkNthDependencyNotNull (typeid (Self), deps, 1);
    checkNthDependencyIs<ConfiguredObject>(typeid (Self), deps, 0);
    if (deps.size() > 1)
      checkNthDependencyIs<ConfiguredParameter>(typeid (Self), deps, 1);

    return cachedAs<Value<Matrix> >(c, std::make_shared<Self>(std::move (deps), dim));
  }


  /* Helper function for generating numerical derivatives of
   * computation nodes. Assuming we have a v = f(object, stuff)
   * node, with v of type T. df/dn = sum_i df/dx_i * dx_i/dn +
   * df/dstuff * dstuff/dn.
   * This function returns a NodeRefVec containings the nodes:
   * {df/dx_i * dx_i/dn} for i in order.
   *
   * buildFWithFreqSet(newFreqSet) should create the f(newFS,
   * stuff) node.
   */

  template<typename ConfiguredObject, typename T, typename B>
  static NodeRefVec generateDerivativeSumDepsForComputations (
    Context& c, ConfiguredObject& object, const Node_DF& derivationNode, const Dimension<T>& targetDimension, B buildFWithNewObject)
  {
    NodeRefVec derivativeSumDeps;

    for (std::size_t i = 0; i < object.nbDependencies (); ++i)
    {
      // First compute dxi_dn. If this maps to a constant 0, do not compute df_dxi at all (costly).
      if (!object.dependency(i))
        continue;

      auto dxi_dn = object.dependency (i)->derive (c, derivationNode);
      if (!dxi_dn->hasNumericalProperty (NumericalProperty::ConstantZero))
      {
        auto buildFWithNewXi = [&c, i, &object, &buildFWithNewObject](std::shared_ptr<ConfiguredParameter> newDep) {
                                 // The sub-graph that will be replicated with shifted inputs is: f(freqset(x_i), stuff)
                                 NodeRefVec newObjectDeps = object.dependencies ();
                                 newObjectDeps[i] = std::move (newDep);
                                 auto newObject = object.recreate (c, std::move (newObjectDeps));
                                 return buildFWithNewObject (std::move (newObject));
                               };

        auto df_dxi = generateNumericalDerivative<T>(
          c, object.config, object.dependency (i), targetDimension, buildFWithNewXi);

        if (dxi_dn->hasNumericalProperty (NumericalProperty::ConstantOne))
          derivativeSumDeps.emplace_back (std::move (df_dxi));
        else
          derivativeSumDeps.emplace_back (CWiseMul<T, std::tuple<double, T> >::create (
                                            c, {std::move (dxi_dn), std::move (df_dxi)}, targetDimension));
      }
    }
    return derivativeSumDeps;
  }
};
} // namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_PARAMETRIZABLE_H
