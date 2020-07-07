//
// File: new_dataflow.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-02-23 00:00:00
// Last modified: 2017-06-08
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

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <algorithm>

#include <Bpp/Exceptions.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowCWiseComputing.h>

static bool enableDotOutput = false;
using namespace bpp;

static void dotOutput(const std::string& testName, const std::vector<const Node_DF*>& nodes)
{
  if (enableDotOutput)
  {
    using bpp::DotOptions;
    writeGraphToDot(
      "debug_" + testName + ".dot", nodes, DotOptions::DetailedNodeInfo | DotOptions::ShowDependencyIndex);
  }
}

using bpp::Dimension;
using bpp::MatrixDimension;

/******************************************************************************
 * Test dataflow basics.
 */
struct DoNothingNode : public Node_DF
{
  DoNothingNode() = default;
  DoNothingNode(NodeRefVec&& deps)
    : Node_DF(std::move(deps))
  {
  }

  using Node_DF::invalidateRecursively; // Make it public to test it

  // Flag to check recomputations
  bool computeCalled = false;
  void compute() final { computeCalled = true; }
};

TEST_CASE("dataflow_dependent_list")
{
  auto l = std::make_shared<DoNothingNode>();
  const auto count = [&l](const Node_DF* dependent) {
    const auto& dependents = l->dependentNodes();
    return std::count(dependents.begin(), dependents.end(), dependent);
  };

  CHECK(l->dependentNodes().empty());

  auto d0 = std::make_shared<DoNothingNode>(NodeRefVec{l});
  CHECK(count(d0.get()) == 1);

  auto d1 = std::make_shared<DoNothingNode>(NodeRefVec{l});
  CHECK(count(d0.get()) == 1);
  CHECK(count(d1.get()) == 1);

  d0.reset(); // Destroy d0
  CHECK(count(d0.get()) == 0);
  CHECK(count(d1.get()) == 1);

  dotOutput("dataflow_dependend_list", {d1.get()});
}

TEST_CASE("dataflow_invariants")
{
  /* Build the following DAG:
   * l0__n0__n1
   * l1_/   /
   * l2____/
   */
  auto l0 = std::make_shared<DoNothingNode>();
  auto l1 = std::make_shared<DoNothingNode>();
  auto l2 = std::make_shared<DoNothingNode>();
  auto n0 = std::make_shared<DoNothingNode>(NodeRefVec{l0, l1});
  auto n1 = std::make_shared<DoNothingNode>(NodeRefVec{n0, l2});

  // Initial state is invalid
  CHECK_FALSE(l0->isValid());
  CHECK_FALSE(l1->isValid());
  CHECK_FALSE(l2->isValid());
  CHECK_FALSE(n0->isValid());
  CHECK_FALSE(n1->isValid());

  n0->computeRecursively();

  // n0 and its dependencies must be valid, and compute must have been called
  CHECK(l0->isValid());
  CHECK(l1->isValid());
  CHECK(n0->isValid());
  CHECK(l0->computeCalled);
  CHECK(l1->computeCalled);
  CHECK(n0->computeCalled);

  l0->computeCalled = false;
  l1->computeCalled = false;
  n0->computeCalled = false;
  n1->computeRecursively();

  // Compute n1: dependencies must be valid. n0 and its deps should not be recomputed (already valid)
  CHECK(l2->isValid());
  CHECK(n1->isValid());
  CHECK(l2->computeCalled);
  CHECK(n1->computeCalled);
  CHECK_FALSE(l0->computeCalled);
  CHECK_FALSE(l1->computeCalled);
  CHECK_FALSE(n0->computeCalled);

  l2->computeCalled = false;
  n1->computeCalled = false;
  l2->invalidateRecursively();

  // l2 and n1 must be invalid, but not n0 and its deps
  CHECK(l0->isValid());
  CHECK(l1->isValid());
  CHECK(n0->isValid());
  CHECK_FALSE(l2->isValid());
  CHECK_FALSE(n1->isValid());

  dotOutput("dataflow_invariants", {n1.get()});
}

TEST_CASE("dataflow_node_basic_errors")
{
  auto doNothing = std::make_shared<DoNothingNode>();

  // Failed conversion
  auto asNodeClass = std::shared_ptr<Node_DF>(doNothing);
  CHECK_THROWS_AS(convertRef<Value<int>>(asNodeClass), Exception&);

  // By default, derive fails
  Context c;
  CHECK_THROWS_AS(asNodeClass->derive(c, *asNodeClass), Exception&);
}

/******************************************************************************
 * Test dataflow numerical nodes.
 */
TEST_CASE("ConstantZero")
{
  Context c;
  auto d = ConstantZero<double>::create(c, Dimension<double>());
  auto m = ConstantZero<Eigen::MatrixXd>::create(c, MatrixDimension(1, 2));

  // Check value
  CHECK(d->getTargetValue() == 0);

  const auto& mValue = m->getTargetValue();
  CHECK(MatrixDimension(mValue) == MatrixDimension(1, 2));
  CHECK(mValue(0, 0) == 0);
  CHECK(mValue(0, 1) == 0);

  // Check derivative
  auto dummy = std::make_shared<DoNothingNode>();
  CHECK(d->deriveAsValue(c, *dummy)->getTargetValue() == 0);
  CHECK(d->deriveAsValue(c, *d)->getTargetValue() == 1);

  dotOutput("ConstantZero", {d.get(), m.get()});
}

TEST_CASE("ConstantOne")
{
  Context c;
  auto d = ConstantOne<double>::create(c, Dimension<double>());
  auto m = ConstantOne<Eigen::MatrixXd>::create(c, MatrixDimension(1, 2));

  // Check value
  CHECK(d->getTargetValue() == 1);

  const auto& mValue = m->getTargetValue();
  CHECK(MatrixDimension(mValue) == MatrixDimension(1, 2));
  CHECK(mValue(0, 0) == 1);
  CHECK(mValue(0, 1) == 1);

  // Check derivative
  auto dummy = std::make_shared<DoNothingNode>();
  CHECK(d->deriveAsValue(c, *dummy)->getTargetValue() == 0);
  CHECK(d->deriveAsValue(c, *d)->getTargetValue() == 1);

  dotOutput("ConstantOne", {d.get(), m.get()});
}

TEST_CASE("NumericConstant")
{
  Context c;
  auto d = NumericConstant<double>::create(c, 42);
  auto m = NumericConstant<Eigen::MatrixXd>::create(c, (Eigen::MatrixXd(1, 2) << 42, -42).finished());

  // Check value
  CHECK(d->getTargetValue() == 42);

  const auto& mValue = m->getTargetValue();
  CHECK(MatrixDimension(mValue) == MatrixDimension(1, 2));
  CHECK(mValue(0, 0) == 42);
  CHECK(mValue(0, 1) == -42);

  // Check derivative
  auto dummy = std::make_shared<DoNothingNode>();
  CHECK(d->deriveAsValue(c, *dummy)->getTargetValue() == 0);
  CHECK(d->deriveAsValue(c, *d)->getTargetValue() == 1);

  dotOutput("NumericConstant", {d.get(), m.get()});
}

TEST_CASE("NumericMutable")
{
  Context c;
  auto d = NumericMutable<double>::create(c, 42);
  auto m = NumericMutable<Eigen::MatrixXd>::create(c, (Eigen::MatrixXd(1, 2) << 42, -42).finished());

  // Check value
  CHECK(d->getTargetValue() == 42);

  const auto& mValue = m->getTargetValue();
  CHECK(MatrixDimension(mValue) == MatrixDimension(1, 2));
  CHECK(mValue(0, 0) == 42);
  CHECK(mValue(0, 1) == -42);

  // Check invalidation on change
  auto dependent = std::make_shared<DoNothingNode>(NodeRefVec{d});
  dependent->computeRecursively();
  CHECK(dependent->isValid());
  d->setValue(-42);
  CHECK(d->getTargetValue() == -42);
  CHECK_FALSE(dependent->isValid());

  // Check derivative
  auto dummy = std::make_shared<DoNothingNode>();
  CHECK(d->deriveAsValue(c, *dummy)->getTargetValue() == 0);
  CHECK(d->deriveAsValue(c, *d)->getTargetValue() == 1);

  dotOutput("NumericMutable", {d.get(), m.get()});
}

TEST_CASE("Convert")
{
  Context c;
  auto d = NumericMutable<double>::create(c, 42);
  auto m = NumericMutable<Eigen::MatrixXd>::create(c, (Eigen::MatrixXd(1, 2) << 42., -42.).finished());

  // Check wrong dependency detection
  CHECK_THROWS_AS((Convert<double, float>::create(c, {}, Dimension<double>())), Exception&);
  CHECK_THROWS_AS((Convert<double, float>::create(c, {nullptr}, Dimension<double>())), Exception&);
  CHECK_THROWS_AS((Convert<double, float>::create(c, {d}, Dimension<double>())), Exception&);
  CHECK_THROWS_AS((Convert<double, double>::create(c, {d, d}, Dimension<double>())), Exception&);

  // Scalar to scalar
  auto f = Convert<float, double>::create(c, {d}, Dimension<float>());
  CHECK(f->getTargetValue() == 42.);

  // Scalar to matrix
  auto m2 = Convert<Eigen::MatrixXd, double>::create(c, {d}, MatrixDimension(1, 2));
  const auto& m2Value = m2->getTargetValue();
  CHECK(MatrixDimension(m2Value) == MatrixDimension(1, 2));
  CHECK(m2Value(0, 0) == 42);
  CHECK(m2Value(0, 1) == 42);

  // Matrix to a Fixed size type (RowVector2d).
  auto m3 = Convert<Eigen::RowVector2d, Eigen::MatrixXd>::create(c, {m}, MatrixDimension(1, 2));
  const auto& m3Value = m3->getTargetValue();
  CHECK(MatrixDimension(m3Value) == MatrixDimension(1, 2));
  CHECK(m3Value(0, 0) == 42);
  CHECK(m3Value(0, 1) == -42);

  // Matrix to matrix with transposition
  auto m4 = Convert<Eigen::MatrixXd, Transposed<Eigen::MatrixXd>>::create(c, {m}, MatrixDimension(1, 2));
  const auto& m4Value = m4->getTargetValue();
  CHECK(MatrixDimension(m4Value) == MatrixDimension(2, 1));
  CHECK(m4Value(0, 0) == 42);
  CHECK(m4Value(1, 0) == -42);

  // Check derivative
  auto dummy = std::make_shared<DoNothingNode>();
  auto df_ddummy = f->deriveAsValue(c, *dummy);
  CHECK(df_ddummy->getTargetValue() == 0.);
  auto df_df = f->deriveAsValue(c, *f);
  CHECK(df_df->getTargetValue() == 1.);
  auto df_dd = f->deriveAsValue(c, *d);
  CHECK(df_dd->getTargetValue() == 1.);

  dotOutput("Convert", {f.get(), m2.get(), m3.get(), df_df.get(), df_dd.get()});
  // Not tested: Constant simplifications
}

TEST_CASE("CWiseAdd")
{
  Context c;
  auto d = NumericMutable<double>::create(c, 42);
  auto m = NumericMutable<Eigen::MatrixXd>::create(c, (Eigen::MatrixXd(1, 2) << 42., -42.).finished());

  // Check wrong dependency detection (tuple<A,B>)
  CHECK_THROWS_AS((CWiseAdd<double, std::tuple<double, float>>::create(c, {d, d}, Dimension<double>())),
                  Exception&);
  CHECK_THROWS_AS((CWiseAdd<double, std::tuple<double, double>>::create(c, {nullptr, d}, Dimension<double>())),
                  Exception&);
  CHECK_THROWS_AS((CWiseAdd<double, std::tuple<double, double>>::create(c, {d}, Dimension<double>())), Exception&);
  CHECK_THROWS_AS((CWiseAdd<double, std::tuple<double, double>>::create(c, {d, d, d}, Dimension<double>())),
                  Exception&);
  // Check wrong dependency detection (reduction<A>)
  CHECK_THROWS_AS((CWiseAdd<double, ReductionOf<double>>::create(c, {nullptr}, Dimension<double>())), Exception&);
  CHECK_THROWS_AS((CWiseAdd<double, ReductionOf<float>>::create(c, {d}, Dimension<double>())), Exception&);

  // Scalar + scalar
  auto add_d_d = CWiseAdd<double, std::tuple<double, double>>::create(c, {d, d}, Dimension<double>());
  auto add_0_d = CWiseAdd<double, ReductionOf<double>>::create(c, {}, Dimension<double>());
  auto add_1_d = CWiseAdd<double, ReductionOf<double>>::create(c, {d}, Dimension<double>());
  auto add_2_d = CWiseAdd<double, ReductionOf<double>>::create(c, {d, d}, Dimension<double>());
  auto add_3_d = CWiseAdd<double, ReductionOf<double>>::create(c, {d, d, d}, Dimension<double>());
  CHECK(add_d_d->getTargetValue() == d->getTargetValue() * 2);
  CHECK(add_0_d->getTargetValue() == 0);
  CHECK(add_1_d->getTargetValue() == d->getTargetValue());
  CHECK(add_2_d->getTargetValue() == d->getTargetValue() * 2);
  CHECK(add_3_d->getTargetValue() == d->getTargetValue() * 3);

  // Matrix + matrix
  auto add_m_m = CWiseAdd<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>::create(c, {m, m}, MatrixDimension(1, 2));
  const auto& add_m_m_value = add_m_m->getTargetValue();
  CHECK(MatrixDimension(add_m_m_value) == MatrixDimension(1, 2));
  CHECK(add_m_m_value(0, 0) == 42 * 2);
  CHECK(add_m_m_value(0, 1) == -42 * 2);

  // Matrix + scalar
  auto add_m_d =
    CWiseAdd<Eigen::MatrixXd, std::tuple<Eigen::MatrixXd, double>>::create(c, {m, d}, MatrixDimension(1, 2));
  const auto& add_m_d_value = add_m_d->getTargetValue();
  CHECK(MatrixDimension(add_m_d_value) == MatrixDimension(1, 2));
  CHECK(add_m_d_value(0, 0) == 42 * 2);
  CHECK(add_m_d_value(0, 1) == 0);

  // Derivatives
  auto dummy = std::make_shared<DoNothingNode>();
  auto& x = add_3_d;
  auto dx_ddummy = x->deriveAsValue(c, *dummy);
  CHECK(dx_ddummy->getTargetValue() == 0);
  auto dx_dx = x->deriveAsValue(c, *x);
  CHECK(dx_dx->getTargetValue() == 1);
  auto dx_dd = x->deriveAsValue(c, *d);
  CHECK(dx_dd->getTargetValue() == 3);

  dotOutput("CWiseAdd",
            {add_d_d.get(),
             add_0_d.get(),
             add_1_d.get(),
             add_2_d.get(),
             add_3_d.get(),
             add_m_m.get(),
             add_m_d.get(),
             dx_dx.get(),
             dx_dd.get()});
  // Not tested: Constant simplifications
}

TEST_CASE("CWiseMul")
{
  Context c;
  auto d = NumericMutable<double>::create(c, 42);
  auto m = NumericMutable<Eigen::MatrixXd>::create(c, (Eigen::MatrixXd(1, 2) << 42., -42.).finished());

  // Check wrong dependency detection (tuple<A,B>)
  CHECK_THROWS_AS((CWiseMul<double, std::tuple<double, float>>::create(c, {d, d}, Dimension<double>())),
                  Exception&);
  CHECK_THROWS_AS((CWiseMul<double, std::tuple<double, double>>::create(c, {nullptr, d}, Dimension<double>())),
                  Exception&);
  CHECK_THROWS_AS((CWiseMul<double, std::tuple<double, double>>::create(c, {d}, Dimension<double>())), Exception&);
  CHECK_THROWS_AS((CWiseMul<double, std::tuple<double, double>>::create(c, {d, d, d}, Dimension<double>())),
                  Exception&);
  // Check wrong dependency detection (reduction<A>)
  CHECK_THROWS_AS((CWiseMul<double, ReductionOf<double>>::create(c, {nullptr}, Dimension<double>())), Exception&);
  CHECK_THROWS_AS((CWiseMul<double, ReductionOf<float>>::create(c, {d}, Dimension<double>())), Exception&);

  // Scalar * scalar
  auto mul_d_d = CWiseMul<double, std::tuple<double, double>>::create(c, {d, d}, Dimension<double>());
  auto mul_0_d = CWiseMul<double, ReductionOf<double>>::create(c, {}, Dimension<double>());
  auto mul_1_d = CWiseMul<double, ReductionOf<double>>::create(c, {d}, Dimension<double>());
  auto mul_2_d = CWiseMul<double, ReductionOf<double>>::create(c, {d, d}, Dimension<double>());
  auto mul_3_d = CWiseMul<double, ReductionOf<double>>::create(c, {d, d, d}, Dimension<double>());
  CHECK(mul_d_d->getTargetValue() == d->getTargetValue() * d->getTargetValue());
  CHECK(mul_0_d->getTargetValue() == 1);
  CHECK(mul_1_d->getTargetValue() == d->getTargetValue());
  CHECK(mul_2_d->getTargetValue() == d->getTargetValue() * d->getTargetValue());
  CHECK(mul_3_d->getTargetValue() == d->getTargetValue() * d->getTargetValue() * d->getTargetValue());

  // Matrix * matrix
  auto mul_m_m = CWiseMul<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>::create(c, {m, m}, MatrixDimension(1, 2));
  const auto& mul_m_m_value = mul_m_m->getTargetValue();
  CHECK(MatrixDimension(mul_m_m_value) == MatrixDimension(1, 2));
  CHECK(mul_m_m_value(0, 0) == 42 * 42);
  CHECK(mul_m_m_value(0, 1) == -42 * -42);

  // Matrix * scalar
  auto mul_m_d =
    CWiseMul<Eigen::MatrixXd, std::tuple<Eigen::MatrixXd, double>>::create(c, {m, d}, MatrixDimension(1, 2));
  const auto& mul_m_d_value = mul_m_d->getTargetValue();
  CHECK(MatrixDimension(mul_m_d_value) == MatrixDimension(1, 2));
  CHECK(mul_m_d_value(0, 0) == 42 * 42);
  CHECK(mul_m_d_value(0, 1) == -42 * 42);

  // Derivatives
  auto dummy = std::make_shared<DoNothingNode>();
  auto& x = mul_3_d;
  auto dx_ddummy = x->deriveAsValue(c, *dummy);
  CHECK(dx_ddummy->getTargetValue() == 0);
  auto dx_dx = x->deriveAsValue(c, *x);
  CHECK(dx_dx->getTargetValue() == 1);
  auto dx_dd = x->deriveAsValue(c, *d);
  CHECK(dx_dd->getTargetValue() == 3 * d->getTargetValue() * d->getTargetValue());

  dotOutput("CWiseMul",
            {mul_d_d.get(),
             mul_0_d.get(),
             mul_1_d.get(),
             mul_2_d.get(),
             mul_3_d.get(),
             mul_m_m.get(),
             mul_m_d.get(),
             dx_dx.get(),
             dx_dd.get()});
  // Not tested: Constant simplifications
}

// Test dataflow node for numerical derivation. Can serve as an example of a simple case.
struct OpaqueTestFunction : public Value<double>
{
  using Self = OpaqueTestFunction;

  NumericalDerivativeConfiguration config{};

  static std::shared_ptr<OpaqueTestFunction> create(Context& c, NodeRefVec&& deps)
  {
    checkDependenciesNotNull(typeid(Self), deps);
    checkDependencyVectorSize(typeid(Self), deps, 2);
    checkDependencyRangeIsValue<double>(typeid(Self), deps, 0, deps.size());
    return cachedAs<Self>(c, std::make_shared<Self>(std::move(deps)));
  }
  OpaqueTestFunction(NodeRefVec&& deps)
    : Value<double>(std::move(deps))
  {
  }

  // OpaqueTestFunction additional arguments = (). Entirely defined by deps.
  bool compareAdditionalArguments(const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive(Context& c, const Node_DF& node) final
  {
    // df/dn = df/dx * dx/dn + df/dy * dy/dn (derivative of multivariable func)
    auto dim = Dimension<double>{}; // Dimension is trivial for this test
    NodeRefVec derivativeSumDeps;
    for (std::size_t i = 0; i < this->nbDependencies(); ++i)
    {
      // First compute dxi_dn. If this maps to a constant 0, do not compute df_dxi at all (costly).
      auto dxi_dn = this->dependency(i)->derive(c, node);
      if (!dxi_dn->hasNumericalProperty(NumericalProperty::ConstantZero))
      {
        auto buildFWithNewXi = [this, i, &c](ValueRef<double> newDep) {
          // Build a duplicate of Self (OpaqueTestFunction) with replaced dependency.
          // The function supports a general case for sub-expressions.
          NodeRefVec newDeps = this->dependencies();
          newDeps[i] = std::move(newDep);
          auto newNode = Self::create(c, std::move(newDeps));
          newNode->config = this->config; // Duplicate config
          return newNode;
        };
        auto df_dxi =
          generateNumericalDerivative<double, double>(c, config, this->dependency(i), dim, dim, buildFWithNewXi);
        derivativeSumDeps.emplace_back(CWiseMul<double, std::tuple<double, double>>::create(
          c, {std::move(df_dxi), std::move(dxi_dn)}, Dimension<double>()));
      }
    }
    return CWiseAdd<double, ReductionOf<double>>::create(c, std::move(derivativeSumDeps), Dimension<double>());
  }

  void compute() final
  {
    auto& result = this->accessValueMutable();
    const auto& x = accessValueConstCast<double>(*this->dependency(0));
    const auto& y = accessValueConstCast<double>(*this->dependency(1));
    result = 3 * x * x + y * y; // 3x^2 + y^2
  }
};

TEST_CASE("numerical_derivation")
{
  Context c;
  auto x = NumericMutable<double>::create(c, 1);
  auto y = NumericMutable<double>::create(c, -1);
  auto f = OpaqueTestFunction::create(c, {x, y});

  auto dummy = std::make_shared<DoNothingNode>();
  auto delta = NumericMutable<double>::create(c, 0.0001);

  // Initial state. Numerical diff not configured, so it should fail.
  CHECK(f->getTargetValue() == 4);
  CHECK_THROWS_AS(f->derive(c, *x), Exception&);

  // Configure and test derivations
  f->config.delta = delta;
  f->config.type = NumericalDerivativeType::ThreePoints;
  auto df_ddummy = f->deriveAsValue(c, *dummy);
  auto df_dx = f->deriveAsValue(c, *x);
  auto df_dy = f->deriveAsValue(c, *y);
  CHECK(df_ddummy->getTargetValue() == 0);
  CHECK(df_dx->getTargetValue() == doctest::Approx(3 * 2 * x->getTargetValue()));
  CHECK(df_dy->getTargetValue() == doctest::Approx(1 * 2 * y->getTargetValue()));

  // Second order
  auto d2f_dx2 = df_dx->deriveAsValue(c, *x);
  CHECK(d2f_dx2->getTargetValue() == doctest::Approx(3 * 2));

  // Multiple dependencies
  auto g = OpaqueTestFunction::create(c, {x, x});
  g->config.delta = delta;
  g->config.type = NumericalDerivativeType::ThreePoints;
  auto dg_dx = g->deriveAsValue(c, *x);
  CHECK(dg_dx->getTargetValue() == doctest::Approx(4 * 2 * x->getTargetValue()));

  dotOutput("numerical_derivation", {f.get(), df_ddummy.get(), df_dx.get(), df_dy.get(), d2f_dx2.get()});
}

int main(int argc, char** argv)
{
  const std::string keyword = "dot_output";
  for (int i = 1; i < argc; ++i)
  {
    if (argv[i] == keyword)
    {
      enableDotOutput = true;
    }
  }
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  return context.run();
}
