//
// File: new_phyl_dataflow.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-19
// Last modified: 2017-04-19
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

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Registry.h>
#include <Bpp/NewPhyl/Topology.h>

#include <fstream>
#include <iosfwd>

using bpp::DF::Node;
using bpp::DF::Value;
using bpp::DF::Parameter;

struct A : public Value<int>::Impl
{
  int i_;
  A(int i)
    : Value<int>::Impl(0)
    , i_(i)
  {
  }

  void compute() override
  {
    int a = i_;
    this->foreachDependencyNode([&a](Node::Impl* n) { a += dynamic_cast<Value<int>::Impl&>(*n).getValue(); });
    this->value_ = a;
  }

  void addDep(Node n)
  {
    n.get().registerNode(this);
    dependencyNodes_.emplace_back(std::move(n));
  }
};

TEST_CASE("test")
{
  bpp::Topology::Tree tree;
  {
    auto a = tree.createNode({}, "A");
    auto b = tree.createNode({}, "B");
    auto c = tree.createNode({a, b});
    auto d = tree.createNode({c});
    tree.rootId() = d;
    std::ofstream file("topology_debug");
    bpp::Topology::debugTree(file, tree);

    auto e = tree.nodeRef(c);
    auto e2 = bpp::Topology::Element (e.asNodeRef().fatherBranch().childNode());
    CHECK (e == e2);
    CHECK (e.hashCode () == e2.hashCode ());
  }

  auto a = Parameter<int>::create (3);
  auto b = Parameter<int>::create (42);

  CHECK (a.getValue() == 3);
  a.setValue(-42);
  CHECK (a.getValue() == -42);

 // auto a = Node::create<A>(1);
 // auto b = Node::create<A>(2);
 // auto c = Node::create<A>(3);
 // auto d = Node::create<A>(4);
 // dynamic_cast<A&>(a.get()).addDep(b);
 // dynamic_cast<A&>(a.get()).addDep(c);
 // dynamic_cast<A&>(d.get()).addDep(a);
 // Value<int> v(d);
 // CHECK(v.getValue() == 10);
 // std::ofstream file("df_debug");
 // bpp::DF::debugDag(file, a);
}
