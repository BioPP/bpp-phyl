//
// File: new_nodespec_instantiation.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-19 00:00:00
// Last modified: 2017-05-29
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
#include <Bpp/NewPhyl/NodeSpecification.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/NewPhyl/TopologyAnnotation.h>
#include <cassert>
#include <fstream>

using namespace bpp;
using DF::Node;
using DF::NodeVec;
using DF::NodeSpecification;
using DF::Value;

class MyParamSpec
{
  // Just return the right param
public:
  MyParamSpec(Topology::Node node, const Topology::NodeMap<DF::Parameter<int>>& params)
    : node_(node)
    , params_(params)
  {
  }

  static std::vector<NodeSpecification> computeDependencies() { return {}; }
  Node buildNode(NodeVec) const { return params_[node_].value(); }
  static std::type_index nodeType() { return typeid(int); } // dummy
  std::string description() const { return "MyParam-N" + std::to_string(node_.nodeId()); }

private:
  Topology::Node node_;
  const Topology::NodeMap<DF::Parameter<int>>& params_;
};

struct SumOp
{
  using ResultType = int;
  using ArgumentType = int;
  static void reset(int& r) { r = 0; }
  static void reduce(int& acc, int i) { acc += i; }
};
using Sum = DF::GenericReductionComputation<SumOp>;

class SumSpec
{
public:
  SumSpec(Topology::Node node, const Topology::NodeMap<DF::Parameter<int>>& params)
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
      node_.foreachChildBranch([this, &deps](Topology::Branch&& branch) {
        deps.emplace_back(SumSpec{std::move(branch).childNode(), params_});
      });
    }
    else
    {
      // Leaf
      deps.emplace_back(MyParamSpec{node_, params_});
    }
    return deps;
  }
  static Node buildNode(NodeVec deps) { return Node::create<Sum>(std::move(deps)); }
  static std::type_index nodeType() { return typeid(Sum); }
  std::string description() const { return "Sum-N" + std::to_string(node_.nodeId()); }

private:
  Topology::Node node_;
  const Topology::NodeMap<DF::Parameter<int>>& params_;
};

TEST_CASE("test")
{
  auto buildTree = Topology::Tree::create();
  auto ta = buildTree->createNode();
  auto tb = buildTree->createNode();
  auto tc = buildTree->createNode({ta, tb});
  auto td = buildTree->createNode({tc});
  buildTree->setRootNodeId(td);
  auto tree = Topology::Tree::finalize(std::move(buildTree));

  std::ofstream ft("topology_debug");
  Topology::debugTree(ft, tree);

  Topology::NodeMap<DF::Parameter<int>> params{tree};
  params.access(ta) = DF::Parameter<int>::create(3);
  params.access(tb) = DF::Parameter<int>::create(42);

  auto sumSpec = NodeSpecification::create<SumSpec>(Topology::Node{tree, tree->rootNodeId()}, params);
  auto partialSumSpec = NodeSpecification::create<SumSpec>(Topology::Node{tree, 0}, params);

  DF::Registry registry;

  Value<int> sum{sumSpec.instantiateWithReuse(registry)};
  CHECK(sum.getValue() == 45);
  params.access (ta)->setValue(-42);
  CHECK(sum.getValue() == 0);

  Value<int> partialSum{partialSumSpec.instantiateWithReuse(registry)};

  std::ofstream fd("df_debug");
  DF::debugNodeSpecInstantiationInRegistry(fd, sumSpec, registry, true);
}
