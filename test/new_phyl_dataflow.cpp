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
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/Registry.h>
#include <Bpp/NewPhyl/Topology.h>

#include <fstream>
#include <iostream>
#include <unordered_map>

using bpp::DF::Node;
using bpp::DF::Value;
using bpp::DF::Parameter;
using bpp::DF::Registry;
using bpp::Topology::Element;

using DataSet = std::unordered_map<std::string, Node>;

struct Sum : public Value<int>::Impl
{
  Sum() = default;

  void compute() override
  {
    int a = 0;
    this->foreachDependencyNode([&a](Node::Impl* n) { a += dynamic_cast<Value<int>::Impl&>(*n).getValue(); });
    this->value_ = a;
  }

  void addDep(Node n) { this->appendDependency(std::move(n)); }

  static Node build(Registry& registry, const Element& element, const DataSet& ds)
  {
    return registry.node<Sum>(element, [&] {
      auto& phyloNode = element.asNodeRef();
      auto dfNode = Node::create<Sum>();
      auto& sumNode = static_cast<Sum&>(dfNode.getImpl());
      if (phyloNode.nbChildBranches() > 0)
      {
        // Internal node
        for (bpp::Topology::IndexType i = 0; i < phyloNode.nbChildBranches(); ++i)
          sumNode.addDep(build(registry, phyloNode.childBranch(i).childNode(), ds));
      }
      else
      {
        // Leaf
        sumNode.addDep(ds.at(phyloNode.name()));
      }
      return dfNode;
    });
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
    auto e2 = bpp::Topology::Element(e.asNodeRef().fatherBranch().childNode());
    CHECK(e == e2);
    CHECK(e.hashCode() == e2.hashCode());
  }

  auto a = Parameter<int>::create(3);
  auto b = Parameter<int>::create(42);

  DataSet ds;
  ds.emplace("A", Node(a));
  ds.emplace("B", Node(b));

  Registry registry;

  Value<int> sum{Sum::build(registry, tree.nodeRef(tree.rootId()), ds)};
  CHECK(sum.getValue() == 45);
  a.setValue(-42);
  CHECK(sum.getValue() == 0);

  Value<int> partialSum{Sum::build(registry, tree.nodeRef(0), ds)};

  std::ofstream file("df_debug");
  bpp::DF::debugRegistry(file, registry);
}
