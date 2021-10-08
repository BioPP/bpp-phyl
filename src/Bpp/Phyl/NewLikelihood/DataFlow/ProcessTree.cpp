//
// File: ProcessTree.cpp
// Authors:
//   Laurent GuÃ©guen
// Created: mardi 11 juin 2019, Ã  09h 39
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
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

#include <Bpp/Graph/AssociationTreeGraphImplObserver.h>
#include <Bpp/Phyl/Model/MixedTransitionModel.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/CollectionNodes.h>
#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>

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

  std::vector<std::shared_ptr<PhyloBranchParam> > vB = tree.getAllEdges();

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
      auto mulref = CWiseMul<double, std::tuple<double, double> >::create (context_, {edge->getBrLen()->dependency(0), rate}, Dimension<double>());
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

  std::vector<std::shared_ptr<PhyloBranchParam> > vB = tree.getAllEdges();

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

    auto confPar = dynamic_cast<ConfiguredParameter*>(parList.getSharedParameter(name).get());
    if (!confPar)
      throw Exception("makeProcessTree: unknown ConfiguredParameter " + name);

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
      throw Exception("ProcessTree::ProcessTree : only simple submodels are used, not combinations. Ask developpers");

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


std::shared_ptr<ProcessTree> ProcessTree::makeProcessTree(CollectionNodes& collection, size_t pNum)
{
  auto& process = collection.getCollection().getSubstitutionProcess(pNum);

  auto& pt = *collection.getProcessTree(process.getTreeNumber());

  ProcessComputationTree tree(process);

  // then process tree with DF objects

  return std::make_shared<ProcessTree>(tree, collection.getModelCollection(), pt);
}

std::shared_ptr<ProcessTree> ProcessTree::makeProcessTree(Context& context, const SubstitutionProcess& process, ParameterList& parList, const std::string& suff)
{
  auto& parTree = process.getParametrizablePhyloTree();

  auto modelColl = makeConfiguredModelCollection(context, process, parList);

  ProcessTree pt(context, parTree, parList, suff); // tree with only branches

  ProcessComputationTree tree(process);

  // then process tree with DF objects

  return std::make_shared<ProcessTree>(tree, modelColl, pt);
}
