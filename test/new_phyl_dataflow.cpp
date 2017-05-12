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
#include <Bpp/NewPhyl/DataFlowTemplates.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/Likelihood.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/NewPhyl/Registry.h>
#include <Bpp/NewPhyl/Topology.h>

#include <fstream>
#include <iostream>
#include <unordered_map>

using bpp::DF::Node;
using bpp::DF::Value;
using bpp::Topology::Element;

struct MyParam : public bpp::DF::Parameter<int>::Impl
{
  MyParam(int i)
    : bpp::DF::Parameter<int>::Impl(i)
  {
  }
};

struct Sum : public bpp::DF::Value<int>::Impl
{
  Sum() = default;

  void compute() override
  {
    int a = 0;
    this->foreachDependencyNode([&a](Node::Impl* n) { a += dynamic_cast<Value<int>::Impl&>(*n).getValue(); });
    this->value_ = a;
  }

  void addDep(Node n) { this->appendDependency(std::move(n)); }

  static std::vector<bpp::DF::NodeSpecification> computeDepencies(const bpp::DF::NodeSpecification& key)
  {
    auto& phyloNode = key.element().asNodeRef();
    std::vector<bpp::DF::NodeSpecification> deps;
    if (phyloNode.nbChildBranches() > 0)
    {
      // Internal node
      for (auto i : bpp::range(phyloNode.nbChildBranches()))
        deps.emplace_back(key.withElement(phyloNode.childBranch(i).childNode()));
    }
    else
    {
      // Leaf
      deps.emplace_back(key.withOperation<MyParam>());
    }
    return deps;
  }
  static Node buildNode(std::vector<Node> deps)
  {
    // TODO automatize ?
    auto node = Node::create<Sum>();
    auto& impl = static_cast<Sum&>(node.getImpl());
    for (auto& d : deps)
      impl.addDep(std::move(d));
    return node;
  }
};

TEST_CASE("test")
{
  bpp::DF::Builder::registerOperation<Sum>(Sum::computeDepencies, Sum::buildNode);

  bpp::Topology::Tree tree;
  auto ta = tree.createNode({}, "A");
  auto tb = tree.createNode({}, "B");
  auto tc = tree.createNode({ta, tb});
  auto td = tree.createNode({tc});
  tree.rootId() = td;
  /*
  std::ofstream ft("topology_debug");
  bpp::Topology::debugTree(ft, tree);
  */

  bpp::DF::Registry registry;

  using Key = bpp::DF::NodeSpecification;

  auto a = bpp::DF::Parameter<int>::create<MyParam>(3);
  auto b = bpp::DF::Parameter<int>::create<MyParam>(42);

  registry.setNode(Key::create<MyParam>(tree.nodeRef(ta)), Node(a));
  registry.setNode(Key::create<MyParam>(tree.nodeRef(tb)), Node(b));

  Value<int> sum{registry.instantiate(Key::create<Sum>(tree.nodeRef(tree.rootId())))};
  CHECK(sum.getValue() == 45);
  a.setValue(-42);
  CHECK(sum.getValue() == 0);

  // Value<int> partialSum{Sum::build(registry, tree.nodeRef(0), ds)};

  std::ofstream fd("df_debug");
  bpp::DF::debugRegistry(fd, registry);
}
