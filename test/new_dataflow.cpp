// File: test_likelihood.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 23/02/2017

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Debug.h>
#include <cassert>
#include <fstream>

using namespace bpp;

struct Sum : public DF::Value<int>::Impl
{
  Sum(std::vector<DF::Node> deps)
    : DF::Value<int>::Impl(std::move(deps))
  {
    // Check deps
    this->foreachDependencyNode([](DF::Node::Impl* n) { assert(dynamic_cast<DF::Value<int>::Impl*>(n)); });
  }

  void compute() override
  {
    int a = 0;
    // Use static downcast as it was checked before (TODO make this more ergonomic)
    this->foreachDependencyNode([&a](DF::Node::Impl* n) { a += static_cast<DF::Value<int>::Impl&>(*n).getValue(); });
    this->value_ = a;
  }
};

struct NegOne : public DF::Value<int>::Impl
{
  NegOne(DF::Node dep)
    : DF::Value<int>::Impl(std::vector<DF::Node>{std::move(dep)})
  {
    // Check deps
    this->foreachDependencyNode([](DF::Node::Impl* n) { assert(dynamic_cast<DF::Value<int>::Impl*>(n)); });
  }

  void compute() override
  {
    // Use static downcast as it was checked before (TODO make this more ergonomic)
    this->value_ = -static_cast<DF::Value<int>::Impl&>(this->dependencyNodes_.front().getImpl()).getValue();
  }
};

TEST_CASE("Testing data flow system on simple int reduction tree")
{
  using namespace bpp::DF;

  /* Build the following DAG:
   * p1__n1__n2__root
   * p2_/   /
   * p3____/__n3
   * p4______/
   *
   * With p_i parameters, n_i sum nodes, root a neg node.
   */

  auto p1 = Parameter<int>::create(42);
  auto p2 = Parameter<int>::create(1);
  auto p3 = Parameter<int>::create(0);
  auto p4 = Parameter<int>::create(3);

  auto n1 = Value<int>::create<Sum>(std::vector<Node>{Node(p1), Node(p2)});
  auto n2 = Value<int>::create<Sum>(std::vector<Node>{Node(n1), Node(p3)});
  auto n3 = Value<int>::create<Sum>(std::vector<Node>{Node(p3), Node(p4)});

  auto root = Value<int>::create<NegOne>(Node(n2));

  // Initial state
  CHECK(p1.getImpl().isValid());
  CHECK_FALSE(n2.getImpl().isValid());
  CHECK_FALSE(root.getImpl().isValid());
  CHECK_FALSE(n3.getImpl().isValid());

  // Get an intermediate value
  CHECK(n2.getValue() == 43);
  CHECK(n2.getImpl().isValid());
  CHECK_FALSE(root.getImpl().isValid());
  CHECK_FALSE(n3.getImpl().isValid());

  // Get root
  CHECK(root.getValue() == -43);
  CHECK(root.getImpl().isValid());
  CHECK_FALSE(n3.getImpl().isValid());

  // Get n3
  CHECK(n3.getValue() == 3);
  CHECK(root.getImpl().isValid());
  CHECK(n3.getImpl().isValid());

  // Change p3, check invalidations
  p3.setValue(10);
  CHECK(p3.getImpl().isValid());
  CHECK_FALSE(root.getImpl().isValid());
  CHECK_FALSE(n3.getImpl().isValid());
  CHECK(n1.getImpl().isValid()); // Not dependent on p3

  // Recompute root
  CHECK(root.getValue() == -53);
  CHECK(root.getImpl().isValid());
  CHECK_FALSE(n3.getImpl().isValid());

  // Recompute n3
  CHECK(n3.getValue() == 13);
  CHECK(n3.getImpl().isValid());

  // Print DF graph
  std::ofstream fd("df_debug");
  debugDag(fd, Node(root));
}
