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
#include <Bpp/NewPhyl/NodeSpecification.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/NewPhyl/Topology.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <unordered_map>

using bpp::DF::Node;
using bpp::DF::NodeSpecification;
using bpp::DF::Value;
using bpp::Topology::NodeRef;

using NodeVec = std::vector<Node>;

class MyParamSpec
{
  // Just return the right param
public:
  MyParamSpec(NodeRef node, NodeVec& params)
    : node_(node)
    , params_(params)
  {
  }

  static std::vector<NodeSpecification> computeDependencies() { return {}; }
  Node buildNode(std::vector<Node>) const { return Node(params_[node_.nodeId()]); }
  static std::type_index nodeType() { return typeid(int); } // dummy
  std::string description() const { return "MyParam-N" + std::to_string(node_.nodeId()); }

private:
  bpp::Topology::NodeRef node_;
  NodeVec& params_;
};

struct Sum : public bpp::DF::Value<int>::Impl
{
  Sum(std::vector<Node> deps)
    : bpp::DF::Value<int>::Impl(std::move(deps))
  {
    // Check deps
    this->foreachDependencyNode([](Node::Impl* n) { assert(dynamic_cast<Value<int>::Impl*>(n)); });
  }

  void compute() override
  {
    int a = 0;
    // Use static downcast as it was checked before (TODO make this more ergonomic)
    this->foreachDependencyNode([&a](Node::Impl* n) { a += static_cast<Value<int>::Impl&>(*n).getValue(); });
    this->value_ = a;
  }

  class Spec
  {
  public:
    Spec(NodeRef node, NodeVec& params)
      : node_(node)
      , params_(params)
    {
    }

    std::vector<NodeSpecification> computeDependencies() const
    {
      std::vector<NodeSpecification> deps;
      if (node_.nbChildBranches() > 0)
      {
        // Internal node
        node_.foreachChildBranch([this, &deps](bpp::Topology::BranchRef&& branch) {
          deps.emplace_back(Spec{std::move(branch).childNode(), params_});
        });
      }
      else
      {
        // Leaf
        deps.emplace_back(MyParamSpec{node_, params_});
      }
      return deps;
    }
    static Node buildNode(std::vector<Node> deps) { return Node::create<Sum>(std::move(deps)); }
    static std::type_index nodeType() { return typeid(Sum); }
    std::string description() const { return "Sum-N" + std::to_string(node_.nodeId()); }

  private:
    bpp::Topology::NodeRef node_;
    NodeVec& params_;
  };
};

TEST_CASE("test")
{
  bpp::Topology::Tree tree;
  auto ta = tree.createNode();
  auto tb = tree.createNode();
  auto tc = tree.createNode({ta, tb});
  auto td = tree.createNode({tc});
  tree.rootId() = td;

  std::ofstream ft("topology_debug");
  bpp::Topology::debugTree(ft, tree);

  auto a = bpp::DF::Parameter<int>::create(3);
  auto b = bpp::DF::Parameter<int>::create(42);
  NodeVec params(tree.nbNodes());
  params[ta] = Node(a);
  params[tb] = Node(b);

  auto sumSpec = NodeSpecification::create<Sum::Spec>(tree.nodeRef(tree.rootId()), params);
  auto partialSumSpec = NodeSpecification::create<Sum::Spec>(tree.nodeRef(0), params);

  bpp::DF::Registry registry;

  Value<int> sum{sumSpec.instantiateWithReuse(registry)};
  CHECK(sum.getValue() == 45);
  a.setValue(-42);
  CHECK(sum.getValue() == 0);

  Value<int> partialSum{partialSumSpec.instantiateWithReuse(registry)};

  std::ofstream fd("df_debug");
  bpp::DF::debugNodeSpecInstantiationInRegistry(fd, sumSpec, registry, true);
}
