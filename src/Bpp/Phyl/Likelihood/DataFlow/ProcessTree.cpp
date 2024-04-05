// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Graph/AssociationTreeGraphImplObserver.h>
#include <Bpp/Phyl/Likelihood/DataFlow/CollectionNodes.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Model/MixedTransitionModel.h>

#include "Model.h"
#include "Parametrizable.h"
#include "ProcessTree.h"

// From the stl:
#include <string>

using namespace bpp;
using namespace std;

ProcessTree::ProcessTree(Context& context,
    const ParametrizablePhyloTree& tree) :
  AssociationTreeGlobalGraphObserver<ProcessNode, ProcessEdge>(tree.getGraph()),
  context_(context)
{
  // Set Nodes
  auto vNodes = tree.getAllNodes();

  for (const auto& node:vNodes)
  {
    uint index = tree.getNodeIndex(node);
    auto pn = make_shared<ProcessNode>(*node, index);
    associateNode(pn, tree.getNodeGraphid(node));
    setNodeIndex(pn, index);
  }

  // Set Edges

  std::vector<std::shared_ptr<PhyloBranchParam>> vB = tree.getAllEdges();

  for (auto& branch:vB)
  {
    const auto& bp = branch->getParameters()[0]; // Get BrLen Parameter

    auto parDF = ConfiguredParameter::create(context, bp);
    auto index = tree.getEdgeIndex(branch);
    auto brref = make_shared<ProcessEdge>(index, parDF, nullptr);

    associateEdge(brref, tree.getEdgeGraphid(branch));
    setEdgeIndex(brref, index);
  }
}

ProcessTree::ProcessTree(const ProcessTree& tree,
    ValueRef<double> rate) :
  AssociationTreeGlobalGraphObserver<ProcessNode, ProcessEdge>(tree)
{
  // Adjust Edges
  auto aEit = allEdgesIterator();

  while (!aEit->end())
  {
    auto edge = **aEit;

    if (edge->getBrLen())
    {
      auto mulref = CWiseMul<double, std::tuple<double, double>>::create (context_, {edge->getBrLen()->dependency(0), rate}, Dimension<double>());
      auto confpar = std::dynamic_pointer_cast<ConfiguredParameter>(edge->getBrLen()->recreate(context_, {std::move(mulref)}));
      edge->setBrLen(confpar);
    }
    aEit->next();
  }
}

ProcessTree::ProcessTree(Context& context,
    const ParametrizablePhyloTree& tree,
    const ParameterList& parList,
    const std::string& suff) :
  AssociationTreeGlobalGraphObserver<ProcessNode, ProcessEdge>(tree.getGraph()),
  context_(context)
{
  // Set Nodes
  auto vNodes = tree.getAllNodes();

  for (const auto& node:vNodes)
  {
    uint index = tree.getNodeIndex(node);
    auto pn = make_shared<ProcessNode>(*node, index);
    associateNode(pn, tree.getNodeGraphid(node));
    setNodeIndex(pn, index);
  }

  // Set Edges

  std::vector<std::shared_ptr<PhyloBranchParam>> vB = tree.getAllEdges();

  for (auto& branch:vB)
  {
    const auto& bp = branch->getParameters()[0]; // Get BrLen Parameter

    std::string name = bp.getName() + suff;
    if (!parList.hasParameter(name) && suff == "")
    {
      if (!parList.hasParameter(bp.getName() + "_1"))
        throw Exception("makeTreeNode: unknown ConfiguredParameter " + name);
      else
        name = bp.getName() + "_1";
    }

    auto confPar = dynamic_cast<ConfiguredParameter*>(parList.getParameter(name).get());
    if (!confPar)
      throw Exception("ProcessTree::ProcessTree: unknown ConfiguredParameter " + name);

    // Share numeric dependency with this parameter
    auto parDF = ConfiguredParameter::create(context, {confPar->dependency(0)}, bp);
    auto index = tree.getEdgeIndex(branch);
    auto brref = make_shared<ProcessEdge>(index, parDF, nullptr);

    associateEdge(brref, tree.getEdgeGraphid(branch));
    setEdgeIndex(brref, index);
  }
}

ProcessTree::ProcessTree(const ProcessComputationTree& tree,
    ParametrizableCollection<ConfiguredModel>& modelColl,
    const ProcessTree& phyloTree) :
  AssociationTreeGlobalGraphObserver<ProcessNode, ProcessEdge>(tree.getGraph()), context_(phyloTree.context_)
{
  // Set Nodes
  auto vNodes = tree.getAllNodes();

  for (const auto& node:vNodes)
  {
    auto pn = make_shared<ProcessNode>(*node);
    associateNode(pn, tree.getNodeGraphid(node));
    setNodeIndex(pn, tree.getNodeIndex(node));
  }

  // Assign References on all branches

  auto vEdges = tree.getAllEdges();

  for (const auto& edge:vEdges)
  {
    std::shared_ptr<ProcessEdge> brref;

    auto id = tree.getEdgeGraphid(edge);
    uint spIndex = edge->getSpeciesIndex(); // index of the matching
    // edge in the
    // ParametrizablePhyloTree

    auto model = edge->getModel();
    if (!model) // ie empty branch
    {
      brref = make_shared<ProcessEdge>(spIndex, nullptr);
      associateEdge(brref, id);
      setEdgeIndex(brref, tree.getEdgeIndex(edge));
      continue;
    }

    if (!modelColl.hasObject(edge->getModelNumber()))
      throw Exception("ProcessTree::ProcessTree : Model unknown " + model->getName() + "  for node " + TextTools::toString(spIndex));

    std::shared_ptr<ConfiguredModel> pmodel = dynamic_pointer_cast<ConfiguredModel>(modelColl[edge->getModelNumber()]);

    auto vNb = edge->subModelNumbers();
    if (vNb.size() > 1)
      throw Exception("ProcessTree::ProcessTree : only simple submodels are used, not combinations. Ask developers");

    if (!edge->useProb()) // model transition is used
    {
      if (!phyloTree.hasEdge(spIndex)) // ie there is no branch length
        throw Exception("ProcessTree::ProcessTree missing branch length for branch " + TextTools::toString(spIndex));

      auto edge2 = phyloTree.getEdge(spIndex);
      if (vNb.size() == 0) // ie full model
        brref = make_shared<ProcessEdge>(spIndex, edge2->getBrLen(), pmodel);
      else
      {
        auto nMod = NumericConstant<size_t>::create(context_, vNb[0]); // Only 1st submodel is used
        brref = make_shared<ProcessEdge>(spIndex, edge2->getBrLen(), pmodel, nMod);
      }
    }
    else // a branch issued from a mixture.
    {
      if (vNb.size() == 0)
        brref = make_shared<ProcessEdge>(spIndex, ConstantOne<double>::create (context_, Dimension<double>()));
      else
      {
        auto nMod = NumericConstant<size_t>::create(context_, (size_t)vNb[0]); // Only 1st submodel is used
        auto brprob = ProbabilityFromMixedModel::create(context_, {pmodel}, (size_t)vNb[0]);
        brref = make_shared<ProcessEdge>(spIndex, brprob);
      }
    }

    associateEdge(brref, id);
    setEdgeIndex(brref, tree.getEdgeIndex(edge));
  }
}


shared_ptr<ProcessTree> ProcessTree::makeProcessTree(CollectionNodes& collection, size_t pNum)
{
  auto process = collection.collection().getSubstitutionProcess(pNum);

  auto& pt = *collection.getProcessTree(process->getTreeNumber());

  ProcessComputationTree tree(process);

  // then process tree with DF objects

  return std::make_shared<ProcessTree>(tree, collection.getModelCollection(), pt);
}

shared_ptr<ProcessTree> ProcessTree::makeProcessTree(
    Context& context,
    shared_ptr<const SubstitutionProcessInterface> process,
    ParameterList& parList,
    const std::string& suff)
{
  auto parTree = process->getParametrizablePhyloTree();

  if (!parTree)
    throw Exception("ProcessTree::makeProcessTree: missing Tree in process.");

  auto modelColl = makeConfiguredModelCollection(context, *process, parList);

  ProcessTree pt(context, *parTree, parList, suff); // tree with only branches

  ProcessComputationTree tree(process);

  // then process tree with DF objects

  return std::make_shared<ProcessTree>(tree, modelColl, pt);
}

DAGindexes ProcessTree::getDAGEdgesIndexes(const Speciesindex speciesIndex) const
{
  auto vB = getAllEdges();

  DAGindexes dag;
  for (auto& branch:vB)
  {
    if (branch->getSpeciesIndex() == speciesIndex)
      dag.push_back(getEdgeIndex(branch));
  }

  return dag;
}
