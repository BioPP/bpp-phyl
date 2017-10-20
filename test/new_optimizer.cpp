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

#include <Bpp/NewPhyl/DataFlowInternalTemplates.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/Optimizer.h>

using namespace bpp::DF;

// Square
struct Square : public Value<double>
{
  using Dependencies = FunctionOfValues<double>;
  Square(NodeRefVec&& deps)
    : Value<double>(std::move(deps))
  {
    checkDependencies(*this);
  }
  void compute() override final
  {
    callWithValues(*this, [](double& r, const double& t) { r = t * t; });
  }
  std::string description() const override final { return "x^2"; }
  NodeRef derive(const Node& variable) override final
  {
    auto& x = this->dependencies()[0];
    return make<MulDouble>({make<Constant<double>>(2.0), x, x->derive(variable)});
  }
};

///////////////// TESTS /////////////////

TEST_CASE("derive constant")
{
  auto konst = make<Constant<double>>(42.0);
  CHECK(konst->getValue() == 42.0);
  CHECK(konst->isConstant());

  auto dummy = make<Parameter<double>>(0);
  auto derived = convertRef<Value<double>>(konst->derive(*dummy));
  CHECK(derived->isConstant());
  CHECK(derived->getValue() == 0.);
}

TEST_CASE("derive parameter")
{
  auto x = make<Parameter<double>>(42.0);
  auto dummy = make<Parameter<double>>(3);

  auto dx_dx = convertRef<Value<double>>(x->derive(*x));
  CHECK(dx_dx->isConstant());
  CHECK(dx_dx->getValue() == 1.);

  auto dx_dummy = convertRef<Value<double>>(x->derive(*dummy));
  CHECK(dx_dummy->isConstant());
  CHECK(dx_dummy->getValue() == 0.);
}

#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Function/SimpleNewtonMultiDimensions.h>
#include <fstream>
#include <iostream>

TEST_CASE("optimizer test")
{
  // x^2 + (y - 3)^2
  bpp::DataFlowParameter xp{"x", 2.0};
  bpp::DataFlowParameter yp{"y", -3.0};

  auto& x = xp.getDataFlowParameter();
  auto& y = yp.getDataFlowParameter();
  auto x2 = make<Square>({x});
  auto konst3 = make<Constant<double>>(-3.);
  auto sy = make<AddDouble>({y, konst3});
  auto y2 = make<Square>({sy});
  auto f = make<AddDouble>({x2, y2});

  bpp::ParameterList params;
  params.addParameter(xp);
  params.addParameter(yp);
  bpp::DataFlowFunction dfFunc{f, params};

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

  CHECK(doctest::Approx(xp.getValue()) == 0.);
  CHECK(doctest::Approx(yp.getValue()) == 3.);
}
