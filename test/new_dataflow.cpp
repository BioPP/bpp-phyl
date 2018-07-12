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

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <algorithm>

#include <Bpp/Exceptions.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>

using namespace bpp::dataflow;
using bpp::MatrixDimension;

/******************************************************************************
 * Test dataflow basics.
 */
struct DoNothingNode : public Node
{
  DoNothingNode() = default;
  DoNothingNode(NodeRefVec&& deps)
    : Node(std::move(deps))
  {
  }

  using Node::invalidateRecursively; // Make it public to test it

  // Flag to check recomputations
  bool computeCalled = false;
  void compute() final { computeCalled = true; }
};

TEST_CASE("dataflow_dependent_list")
{
  auto l = std::make_shared<DoNothingNode>();
  const auto count = [&l](const Node* dependent) {
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
}

TEST_CASE("dataflow_node_basic_errors")
{
  auto doNothing = std::make_shared<DoNothingNode>();

  // Failed conversion
  auto asNodeClass = std::shared_ptr<Node>(doNothing);
  CHECK_THROWS_AS(convertRef<Value<int>>(asNodeClass), bpp::Exception);

  // By default, derive fails
  Context c;
  CHECK_THROWS_AS(asNodeClass->derive(c, *asNodeClass), bpp::Exception);
}

/******************************************************************************
 * Test dataflow numerical nodes.
 */
TEST_CASE("ConstantZero")
{
  Context c;

  auto d = ConstantZero<double>::create(c);
  auto m = ConstantZero<Eigen::MatrixXd>::create(c, MatrixDimension(1, 2));

  // Check value
  CHECK(d->getValue() == 0.);

  const auto& mValue = m->getValue();
  CHECK(MatrixDimension(mValue) == MatrixDimension(1, 2));
  CHECK(mValue(0, 0) == 0.);
  CHECK(mValue(0, 1) == 0.);

  // Check derivative
  auto dummy = std::make_shared<DoNothingNode>();
  CHECK(d->deriveAsValue(c, *dummy)->getValue() == 0.);
}

TEST_CASE("ConstantOne")
{
  Context c;

  auto d = ConstantOne<double>::create(c);
  auto m = ConstantOne<Eigen::MatrixXd>::create(c, MatrixDimension(1, 2));

  // Check value
  CHECK(d->getValue() == 1.);

  const auto& mValue = m->getValue();
  CHECK(MatrixDimension(mValue) == MatrixDimension(1, 2));
  CHECK(mValue(0, 0) == 1.);
  CHECK(mValue(0, 1) == 1.);

  // Check derivative
  auto dummy = std::make_shared<DoNothingNode>();
  CHECK(d->deriveAsValue(c, *dummy)->getValue() == 0.);
}

TEST_CASE("test")
{
  using namespace bpp::dataflow;

  Context context;

  auto a = ConstantZero<double>::create(context);
  auto b = ConstantZero<Eigen::MatrixXd>::create(context, bpp::MatrixDimension{42, 32});
  auto c = ConstantOne<Eigen::VectorXd>::create(context, bpp::vectorDimension(42));

  auto d = NumericConstant<Eigen::MatrixXd>::create(context, Eigen::MatrixXd::Random(42, 32));
  auto e = CWiseAdd<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>::create(context, {b, d});
  e->getValue();

  auto z = e->derive(context, *a);
  // Print DF graph
  // bpp::debugDag("df_debug", *root);
}
