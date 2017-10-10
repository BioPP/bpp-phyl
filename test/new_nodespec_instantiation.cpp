//
// File: new_nodespec_instantiation.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-19 00:00:00
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

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowTemplateUtils.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/FrozenPtr.h>
#include <Bpp/NewPhyl/NodeSpecification.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/NewPhyl/TopologyMap.h>
#include <cassert>
#include <fstream>

using namespace bpp;
using DF::Node;
using DF::Value;

struct Sum : public DF::Value<int>
{
  using Dependencies = DF::ReductionOfValue<int>;
  Sum(DF::NodeRefVec&& deps)
    : DF::Value<int>(std::move(deps))
  {
    DF::checkDependencies(*this);
  }
  void compute() override final
  {
    DF::callWithValues(*this, [](int& r) { r = 0; }, [](int& r, int i) { r += i; });
  }
};

struct SumSpec : DF::NodeSpecAlwaysGenerate<Sum>
{
  Topology::Node node;
  FrozenPtr<Topology::NodeValueMap<DF::ParameterRef<int>>> params;

  SumSpec(const Topology::Node& n, const FrozenPtr<Topology::NodeValueMap<DF::ParameterRef<int>>>& p)
    : node(n)
    , params(p)
  {
  }

  DF::NodeSpecificationVec computeDependencies() const
  {
    if (node.nbChildBranches() > 0)
    {
      // Internal node
      DF::NodeSpecificationVec deps;
      node.foreachChildBranch([this, &deps](Topology::Branch&& branch) {
        deps.emplace_back(SumSpec{std::move(branch).childNode(), params});
      });
      return deps;
    }
    else
    {
      // Leaf
      return DF::makeNodeSpecVec(DF::NodeSpecReturnParameter{params->access(node).value()});
    }
  }
  std::string description() const { return "Sum-N" + std::to_string(node.nodeId()); }
};

TEST_CASE("test")
{
  auto buildTree = make_freezable<Topology::Tree>();
  auto ta = buildTree->createNode();
  auto tb = buildTree->createNode();
  auto tc = buildTree->createNode({ta, tb});
  auto td = buildTree->createNode({tc});
  buildTree->setRootNodeId(td);
  auto tree = std::move(buildTree).freeze();

  std::ofstream ft("topology_debug");
  Topology::debugTree(ft, tree);

  auto buildParams = make_freezable<Topology::NodeValueMap<DF::ParameterRef<int>>>(tree);
  buildParams->access(tree->node(ta)) = DF::Parameter<int>::create(3);
  buildParams->access(tree->node(tb)) = DF::Parameter<int>::create(42);
  auto params = std::move(buildParams).freeze();

  DF::Registry registry;

  auto sumSpec = SumSpec{tree->rootNode(), params};
  auto sum = DF::convertRef<Value<int>>(DF::instantiateNodeSpecWithReuse(sumSpec, registry));
  CHECK(sum->getValue() == 45);
  params->access(tree->node(ta)).value()->setValue(-42);
  CHECK(sum->getValue() == 0);

  auto partialSum =
    DF::convertRef<Value<int>>(DF::instantiateNodeSpecWithReuse(SumSpec{tree->node(0), params}, registry));

  std::ofstream fd("df_debug");
  DF::debugNodeSpecInstantiationInRegistry(fd, sumSpec, registry, DF::DebugOptions::ShowRegistryLinks);
}
