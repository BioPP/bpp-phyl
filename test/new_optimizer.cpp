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

#include <Bpp/NewPhyl/DataFlowTemplateUtils.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/Optimizer.h>

#include <algorithm> // remove_if for Add optimisation

using namespace bpp;
using namespace bpp::DF;

// Addition
template<typename T>
NodeRef makeAdditionNode(NodeRefVec deps);

template<typename T>
struct Add : public Value<T>
{
  using Dependencies = ReductionOfValue<T>;
  Add(NodeRefVec&& deps)
    : Value<T>(std::move(deps))
  {
    checkDependencies(*this);
  }
  void compute() override final
  {
    DF::callWithValues(*this, [](T& r) { r = 0.; }, [](T& r, const T& t) { r += t; });
  }
  std::string description() const override final { return "+"; }
  NodeRef derive(const Node& variable) override final
  {
    NodeRefVec derivatives;
    for (auto& subExpr : this->dependencies())
      derivatives.emplace_back(subExpr->derive(variable));
    return makeAdditionNode<T>(std::move(derivatives));
  }
};

// Multiplication
template<typename T>
NodeRef makeMultiplicationNode(NodeRefVec deps);

template<typename T>
struct Mul : public Value<T>
{
  using Dependencies = ReductionOfValue<T>;
  Mul(NodeRefVec&& deps)
    : Value<T>(std::move(deps))
  {
    checkDependencies(*this);
  }
  void compute() override final
  {
    DF::callWithValues(*this, [](T& r) { r = 1.; }, [](T& r, const T& t) { r *= t; });
  }
  std::string description() const override final { return "*"; }
  NodeRef derive(const Node& variable) override final
  {
    NodeRefVec additionDeps;
    for (auto i : bpp::index_range(this->dependencies()))
    {
      NodeRefVec mulDeps = this->dependencies();
      mulDeps[i] = this->dependencies()[i]->derive(variable);
      additionDeps.emplace_back(makeMultiplicationNode<T>(std::move(mulDeps)));
    }
    return makeAdditionNode<T>(std::move(additionDeps));
  }
};

// Square
template<typename T>
struct Square : public Value<T>
{
  using Dependencies = FunctionOfValues<T>;
  Square(NodeRefVec&& deps)
    : Value<T>(std::move(deps))
  {
    checkDependencies(*this);
  }
  void compute() override final
  {
    DF::callWithValues(*this, [](T& r, const T& t) { r = t * t; });
  }
  std::string description() const override final { return "x^2"; }
  NodeRef derive(const Node& variable) override final
  {
    auto& x = this->dependencies()[0];
    return makeMultiplicationNode<T>({createNode<Constant<T>>(2.0), x, x->derive(variable)});
  }
};

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
    return createNode<Constant<T>>(createZeroValue(savedRef->accessValue()));
  }
  else
  {
    return createNode<Add<T>>(std::move(deps));
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
    return createNode<Constant<T>>(createZeroValue(savedRef->accessValue()));
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
    return createNode<Constant<T>>(createOneValue(savedRef->accessValue()));
  }
  else
  {
    return createNode<Mul<T>>(std::move(deps));
  }
}

///////////////// TESTS /////////////////::

// FIXME fix them later when API is stable
#if 0
TEST_CASE("derive constant")
{
  auto konst = createNode<Constant<double>>(42.0);
  CHECK(konst->getValue () == 42.0);
  CHECK(konst->isConstant());

  auto dummy = createNode<Parameter<double>>(0);
  auto derived = convertRef<Value<double>>(konst->derive(*dummy));
  CHECK(derived->isConstant());
  CHECK(derived->getValue () == 0);
}

TEST_CASE("derive parameter")
{
  auto x = createNode<Parameter<double>>(42.0);
  auto dummy = createNode<Parameter<double>>(3);

  auto dx_dx = convertRef<Value<double>>(x->derive(*x));
  CHECK(dx_dx->isConstant());
  CHECK(dx_dx->getValue () == 1.0);

  auto dx_dummy = convertRef<Value<double>>(x->derive(*dummy));
  CHECK(dx_dummy->isConstant());
  CHECK(dx_dummy->getValue () == 0.0);
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
  auto x2 = createNode<Square<double>>({x});
  auto konst3 = createNode<Constant<double>>(-3.);
  auto sy = makeAdditionNode<double>({y, konst3});
  auto y2 = createNode<Square<double>>({sy});
  auto f = makeAdditionNode<double>({x2, y2});

#if 0
  std::cout << "x2 + y2 = " << f->getValue () << "\n";
  auto df_dx = f->derive(*x);
  auto df2_dx2 = df_dx->derive(*x);
  std::cout << "d2(x2 + y2)/dx2 = " << convertRef<Value<double>>(df2_dx2)->getValue () << "\n";
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
