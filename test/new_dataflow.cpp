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

#include <Bpp/Exceptions.h>
#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowInternal.h>
#include <Bpp/NewPhyl/DataFlowTemplates.h>
#include <Bpp/NewPhyl/Debug.h>
#include <fstream>

using namespace bpp::DF;

struct AddInt : public Value<int>
{
  using Dependencies = ReductionOfValue<int>;
  AddInt(NodeRefVec&& deps)
    : Value<int>(std::move(deps))
  {
    checkDependencies(*this);
  }
  void compute() final
  {
    callWithValues(*this, [](int& r) { r = 0; }, [](int& r, int i) { r += i; });
  }
  NodeRef rebuild(NodeRefVec&& deps) const final { return makeNode<AddInt>(std::move(deps)); }
};

struct NegInt : public Value<int>
{
  using Dependencies = FunctionOfValues<int>;
  NegInt(NodeRefVec&& deps)
    : Value<int>(std::move(deps))
  {
    checkDependencies(*this);
  }
  void compute() final
  {
    callWithValues(*this, [](int& r, int i) { r = -i; });
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

  auto p1 = makeNode<Mutable<int>>(42);
  auto p2 = makeNode<Mutable<int>>(1);
  auto p3 = makeNode<Mutable<int>>(0);
  auto p4 = makeNode<Mutable<int>>(3);

  auto n1 = makeNode<AddInt>({p1, p2});
  auto n2 = makeNode<AddInt>({n1, p3});
  auto n3 = makeNode<AddInt>({p3, p4});

  auto root = makeNode<NegInt>({n2});

  // Initial state
  CHECK(p1->isValid());
  CHECK_FALSE(n2->isValid());
  CHECK_FALSE(root->isValid());
  CHECK_FALSE(n3->isValid());

  // Get an intermediate value
  CHECK(n2->getValue() == 43);
  CHECK(n2->isValid());
  CHECK_FALSE(root->isValid());
  CHECK_FALSE(n3->isValid());

  // Get root
  CHECK(root->getValue() == -43);
  CHECK(root->isValid());
  CHECK_FALSE(n3->isValid());

  // Get n3
  CHECK(n3->getValue() == 3);
  CHECK(root->isValid());
  CHECK(n3->isValid());

  // Change p3, check invalidations
  p3->setValue(10);
  CHECK(p3->isValid());
  CHECK_FALSE(root->isValid());
  CHECK_FALSE(n3->isValid());
  CHECK(n1->isValid()); // Not dependent on p3

  // Recompute root
  CHECK(root->getValue() == -53);
  CHECK(root->isValid());
  CHECK_FALSE(n3->isValid());

  // Recompute n3
  CHECK(n3->getValue() == 13);
  CHECK(n3->isValid());

  // Print DF graph
  std::ofstream fd("df_debug");
  bpp::debugDag(fd, root);
}

template<typename T>
struct MyParam : bpp::DF::Mutable<T>
{
  using bpp::DF::Mutable<T>::Mutable;
  using bpp::DF::Mutable<T>::invalidateRecursively;
  using bpp::DF::Mutable<T>::computeRecursively;
};

TEST_CASE("Exceptions")
{
  using namespace bpp::DF;

  // Check that param should crash if made invalid
  auto param = makeNode<MyParam<int>>(42);
  param->invalidateRecursively(); // Bad !
  CHECK_THROWS_AS(param->computeRecursively(), bpp::Exception);

  // Bad node dynamic downcast
  CHECK_THROWS_AS((void)convertRef<Value<bool>>(param), bpp::Exception);

  // GenericFunctionComputation: bad dep vec len
  CHECK_THROWS_AS(makeNode<NegInt>({}), bpp::Exception);

  // GenericFunctionComputation: type mismatch
  auto p = makeNode<Mutable<bool>>();
  CHECK_THROWS_AS(makeNode<NegInt>({p}), bpp::Exception);

  // GenericReductionComputation: type mismatch
  CHECK_THROWS_AS(makeNode<AddInt>({p}), bpp::Exception);
}

TEST_CASE("Properties")
{
  using namespace bpp::DF;
  auto konst = makeNode<Constant<int>>(42);
  auto param = makeNode<Mutable<int>>(42);
  auto add = makeNode<AddInt>({konst, param});

  // int not derivable
  CHECK_FALSE(param->isDerivable(*param));
  CHECK_FALSE(param->isDerivable(*konst));
  CHECK_FALSE(add->isDerivable(*param));
  CHECK_THROWS_AS(param->derive(*param), bpp::Exception);

  // isConstant
  CHECK_FALSE(param->isConstant());
  CHECK(konst->isConstant());
  CHECK_FALSE(add->isConstant());

  // isTransitivelyDependentOn
  CHECK_FALSE(konst->isTransitivelyDependentOn(*param));
  CHECK(konst->isTransitivelyDependentOn(*konst));
  CHECK(add->isTransitivelyDependentOn(*konst));
  CHECK_FALSE(konst->isTransitivelyDependentOn(*add));
}

TEST_CASE("Rebuild")
{
  using namespace bpp::DF;
  auto k1 = makeNode<Constant<int>>(42);
  auto k2 = makeNode<Constant<int>>(-21);
  auto k3 = makeNode<Constant<int>>(0);

  auto a12 = makeNode<AddInt>({k1, k2});
  auto a123 = makeNode<AddInt>({a12, k3});
  auto na123 = makeNode<NegInt>({a123});

  CHECK(na123->getValue() == -21);

  auto p = makeNode<Mutable<int>>(-21);
  std::map<const Node*, NodeRef> replace_k3_with_p = {{k3.get(), p}};

  CHECK_THROWS_AS(rebuildWithSubstitution(na123, replace_k3_with_p),
                  bpp::Exception); // NegInt has no support for rebuild

  auto rebuilt_k1 = rebuildWithSubstitution(k1, replace_k3_with_p);
  CHECK(rebuilt_k1 == k1); // Kept the same, leaves are unchanged
  auto rebuilt_a12 = rebuildWithSubstitution(a12, replace_k3_with_p);
  CHECK(rebuilt_a12 == a12); // Not rebuilt, did not depend on k3

  auto rebuilt_a123 = rebuildWithSubstitution(a123, replace_k3_with_p);
  CHECK(rebuilt_a123 != a123); // Rebuilt, did depend on k3
  CHECK(convertRef<Value<int>>(rebuilt_a123)->getValue() == 0);
}
