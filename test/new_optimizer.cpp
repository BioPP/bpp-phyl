//
// File: new_optimizer.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-07-07 00:00:00
// Last modified: 2017-07-07
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

  This software is a computer program whose purpose is to provide classes
  for numerical calculus. This file is part of the Bio++ project.

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

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <Bpp/NewPhyl/DataFlowTemplates.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/Optimizer.h>

#include <algorithm> // remove_if for Add optimisation
#include <cassert>

using namespace bpp::DF;

// Addition
template<typename T>
struct AdditionOp : public OperationBase<AdditionOp<T>>
{
  using ArgumentType = T;
  using ResultType = T;
  static void reset(ResultType& r) { r = 0.; }
  static void reduce(ResultType& r, const ArgumentType& x) { r += x; }
  static NodeRef derive(Node& self, const Node& variable);
  static std::string description() { return "+"; }
};
template<typename T>
using AdditionNode = GenericReductionComputation<AdditionOp<T>>;
template<typename T>
NodeRef makeAdditionNode(NodeRefVec deps);

template<typename T>
NodeRef AdditionOp<T>::derive(Node& self, const Node& variable)
{
  NodeRefVec derivatives;
  for (auto& subExpr : self.dependencies())
    derivatives.emplace_back(subExpr->derive(variable));
  return makeAdditionNode<T>(std::move(derivatives));
}

// Multiplication
template<typename T>
struct MultiplicationOp : public OperationBase<MultiplicationOp<T>>
{
  using ArgumentType = T;
  using ResultType = T;
  static void reset(ResultType& r) { r = 1.; }
  static void reduce(ResultType& r, const ArgumentType& x) { r *= x; }
  static NodeRef derive(Node& self, const Node& variable);
  static std::string description() { return "*"; }
};
template<typename T>
using MultiplicationNode = GenericReductionComputation<MultiplicationOp<T>>;
template<typename T>
NodeRef makeMultiplicationNode(NodeRefVec deps);

template<typename T>
NodeRef MultiplicationOp<T>::derive(Node& self, const Node& variable)
{
  NodeRefVec additionDeps;
  for (auto i : bpp::index_range(self.dependencies()))
  {
    NodeRefVec mulDeps = self.dependencies();
    mulDeps[i] = self.dependencies()[i]->derive(variable);
    additionDeps.emplace_back(makeMultiplicationNode<T>(std::move(mulDeps)));
  }
  return makeAdditionNode<T>(std::move(additionDeps));
}

// Square
template<typename T>
struct SquareOp : public OperationBase<SquareOp<T>>
{
  using ArgumentTypes = std::tuple<T>;
  using ResultType = T;
  static void compute(ResultType& r, const T& x) { r = x * x; }
  static NodeRef derive(Node& self, const Node& variable);
  static std::string description() { return "x^2"; }
};
template<typename T>
using SquareNode = GenericFunctionComputation<SquareOp<T>>;
template<typename T>
NodeRef SquareOp<T>::derive(Node& self, const Node& variable)
{
  auto& x = self.dependencies()[0];
  return makeMultiplicationNode<T>({createNode<Constant<T>>(2.0), x, x->derive(variable)});
}

// Optimisations
template<typename T>
NodeRef makeAdditionNode(NodeRefVec deps)
{
  // FIXME crutch. need to store a ref to have at least one node to get build arguments from.
  assert(deps.size() > 0);
  auto savedRef = convertRef<Value<T>>(deps[0]);
  // Remove '0s' from deps
  deps.erase(std::remove_if(deps.begin(),
                            deps.end(),
                            [](const NodeRef& nodeRef) { return nodeRef->numericProperties().isConstantZero; }),
             deps.end());
  // Node choice
  if (deps.size() == 1)
  {
    return deps[0];
  }
  else if (deps.size() == 0)
  {
    return createNode<Constant<T>>(createZeroValue(savedRef->value()));
  }
  else
  {
    return createNode<AdditionNode<T>>(std::move(deps));
  }
}

template<typename T>
NodeRef makeMultiplicationNode(NodeRefVec deps)
{
  // Same crutch again
  assert(deps.size() > 0);
  auto savedRef = convertRef<Value<T>>(deps[0]);
  // Return 0 if any dep is 0
  if (std::any_of(
        deps.begin(), deps.end(), [](const NodeRef& nodeRef) { return nodeRef->numericProperties().isConstantZero; }))
  {
    return createNode<Constant<T>>(createZeroValue(savedRef->value()));
  }
  // Remove any 1s
  deps.erase(std::remove_if(deps.begin(),
                            deps.end(),
                            [](const NodeRef& nodeRef) { return nodeRef->numericProperties().isConstantOne; }),
             deps.end());
  // Node choice
  if (deps.size() == 1)
  {
    return deps[0];
  }
  else if (deps.size() == 0)
  {
    return createNode<Constant<T>>(createOneValue(savedRef->value()));
  }
  else
  {
    return createNode<MultiplicationNode<T>>(std::move(deps));
  }
}

///////////////// TESTS /////////////////::

// FIXME fix them later when API is stable
#if 0
TEST_CASE("derive constant")
{
  auto konst = createNode<Constant<double>>(42.0);
  CHECK(getUpToDateValue(konst) == 42.0);
  CHECK(konst->isConstant());

  auto dummy = createNode<Parameter<double>>(0);
  auto derived = convertRef<Value<double>>(konst->derive(*dummy));
  CHECK(derived->isConstant());
  CHECK(getUpToDateValue(derived) == 0);
}

TEST_CASE("derive parameter")
{
  auto x = createNode<Parameter<double>>(42.0);
  auto dummy = createNode<Parameter<double>>(3);

  auto dx_dx = convertRef<Value<double>>(x->derive(*x));
  CHECK(dx_dx->isConstant());
  CHECK(getUpToDateValue(dx_dx) == 1.0);

  auto dx_dummy = convertRef<Value<double>>(x->derive(*dummy));
  CHECK(dx_dummy->isConstant());
  CHECK(getUpToDateValue(dx_dummy) == 0.0);
}
#endif

#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Function/SimpleNewtonMultiDimensions.h>
#include <fstream>
#include <iostream>

TEST_CASE("test")
{
  // x^2 + (y - 3)^2
  bpp::DataFlowParameter xp{"x", 2.0};
  bpp::DataFlowParameter yp{"y", -3.0};

  auto& x = xp.getDataFlowParameter();
  auto& y = yp.getDataFlowParameter();
  auto x2 = createNode<SquareNode<double>>({x});
  auto konst3 = createNode<Constant<double>>(-3.);
  auto sy = makeAdditionNode<double>({y, konst3});
  auto y2 = createNode<SquareNode<double>>({sy});
  auto f = makeAdditionNode<double>({x2, y2});

#if 0
  std::cout << "x2 + y2 = " << getUpToDateValue(f) << "\n";
  auto df_dx = f->derive(*x);
  auto df2_dx2 = df_dx->derive(*x);
  std::cout << "d2(x2 + y2)/dx2 = " << getUpToDateValue(convertRef<Value<double>>(df2_dx2)) << "\n";
  std::ofstream fd("df_debug");
  debugDag(fd, df2_dx2, DebugOptions::ShowDependencyIndex | DebugOptions::FollowUpwardLinks);
#endif

  bpp::ParameterList params;
  params.addParameter(xp);
  params.addParameter(yp);
  bpp::DataFlowFunction dfFunc{convertRef<Value<double>>(f), params};

  // bpp::ConjugateGradientMultiDimensions optimizer(&dfFunc);
  bpp::SimpleNewtonMultiDimensions optimizer(&dfFunc);
  optimizer.setVerbose(1);
  optimizer.setProfiler(bpp::ApplicationTools::message);
  optimizer.setMessageHandler(bpp::ApplicationTools::message);
  optimizer.setMaximumNumberOfEvaluations(1000000);
  optimizer.getStopCondition()->setTolerance(0.000001);
  optimizer.setConstraintPolicy(bpp::AutoParameter::CONSTRAINTS_AUTO);
  optimizer.init(dfFunc.getParameters());
  optimizer.optimize();

  std::cout << "(x, y) == (" << xp.getValue() << ", " << yp.getValue() << ")\n";

  std::ofstream fd("df_debug");
  debugDag(fd, dfFunc.getAllNamedNodes("f"), DebugOptions::ShowDependencyIndex);
}
