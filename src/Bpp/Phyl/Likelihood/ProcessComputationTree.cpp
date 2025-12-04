// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <numeric>

#include "ProcessComputationTree.h"

using namespace bpp;
using namespace std;

ProcessComputationTree::ProcessComputationTree(
    shared_ptr<const SubstitutionProcessInterface> process) :
  BaseTree(process->getModelScenario() ? 0 :
      (process->getParametrizablePhyloTree() ? process->getParametrizablePhyloTree()->getGraph() : 0)),
  process_(process)
{
  shared_ptr<const ParametrizablePhyloTree> ptree = process->getParametrizablePhyloTree();
  if (!ptree)
    throw Exception("ProcessComputationTree::ProcessComputationTree: missing tree.");

  // if no model scenario, copy the basic tree
  auto scenario = process->getModelScenario();

  if (!scenario)
  {
    auto itN = ptree->allNodesIterator();
    for (itN->start(); !itN->end(); itN->next())
    {
      auto index = ptree->getNodeIndex(*(*itN));
      auto nnode = make_shared<ProcessComputationNode>(*(*(*itN)), index);
      nnode->setProperty("event", NodeEvent::speciationEvent);
      associateNode(nnode, ptree->getNodeGraphid(*(*itN)));
      setNodeIndex(nnode, index);
    }

    auto itE = ptree->allEdgesIterator();
    for (itE->start(); !itE->end(); itE->next())
    {
      auto index = ptree->getEdgeIndex(*(*itE));
      auto model = process->getModelForNode(index);
      size_t nmodel = process->getModelNumberForNode(index);
      auto nedge = make_shared<ProcessComputationEdge>(model, nmodel, index);
      associateEdge(nedge, ptree->getEdgeGraphid(*(*itE)));
      setEdgeIndex(nedge, index);
    }

    return;
  }

  // Map of the mrca of the MixedTransitionModel split in several paths
  map<shared_ptr<MixedTransitionModelInterface>, unsigned int> mMrca;

  auto vMod = scenario->getModels();
  map<shared_ptr<MixedTransitionModelInterface>, vector<shared_ptr<PhyloNode>>> mnodes;

  auto vNodes = ptree->getAllNodes();

  auto root = ptree->getRoot();

  // first the nodes that carry the models
  for (const auto& node:vNodes)
  {
    if (node == root)
      continue;

    const auto medge = process_->getModelForNode(ptree->getNodeIndex(node));
    shared_ptr<MixedTransitionModelInterface> mok(0);
    for (const auto& mod:vMod)
    {
      if (mod.get() == medge.get())
      {
        mok = mod;
        break;
      }
    }
    if (mok == 0)
      continue;

    if (mnodes.find(mok) == mnodes.end())
      mnodes[mok] = vector<shared_ptr<PhyloNode>>();
    mnodes[mok].push_back(node);
  }

  // then the mrca

  for (const auto& mnode:mnodes)
  {
    auto nrca = ptree->MRCA(mnode.second);
    mMrca[mnode.first] = ptree->getNodeIndex(nrca);
  }

  // Then construction of the tree

  auto nroot = make_shared<ProcessComputationNode>(*root, ptree->getRootIndex());
  createNode(nroot);
  addNodeIndex(nroot);
  buildFollowingScenario_(nroot, *scenario, mMrca);

  rootAt(nroot);
}


void ProcessComputationTree::buildFollowingScenario_(
    shared_ptr<ProcessComputationNode> father,
    const ModelScenario& scenario,
    map<shared_ptr<MixedTransitionModelInterface>,
    unsigned int>& mMrca)
{
  auto spInd = father->getSpeciesIndex();

  size_t nbpath = scenario.getNumberOfModelPaths();

  auto ptree = process_->getParametrizablePhyloTree();

  shared_ptr<MixedTransitionModelInterface> mrca(0); // One MixedModel which is
  // "split" at father node

  for (const auto& mod:mMrca)
  {
    if (mod.second == spInd)
    {
      mrca = mod.first;
      mMrca.erase(mod.first);
      break;
    }
  }


  if (mrca == 0) // it is a speciation node
  {
    // set the event
    father->setProperty("event", NodeEvent::speciationEvent);

    // and then look at the branches
    auto vbrInd = ptree->getBranches(spInd);

    for (const auto& brInd:vbrInd)
    {
      auto sonInd = ptree->getSon(brInd);
      auto son = ptree->getNode(sonInd);

      auto modSon = process_->getModelForNode(sonInd);
      size_t nmodel = process_->getModelNumberForNode(sonInd);

      auto sonnode = make_shared<ProcessComputationNode>(*son, sonInd);
      createNode(sonnode);
      addNodeIndex(sonnode);

      // Check if model is a mixture
      auto mixMod = dynamic_pointer_cast<const MixedTransitionModelInterface>(modSon);

      if (!mixMod) // No mixture, then a branch with a full model
      {
        auto nedge = make_shared<ProcessComputationEdge>(modSon, nmodel, sonInd);
        link(father, sonnode, nedge);
        addEdgeIndex(nedge);
        buildFollowingScenario_(sonnode, scenario, mMrca);
        continue;
      }

      // If it is a mixture model, check its decomposition in the modelpaths

      map<Vuint, vector<shared_ptr<ModelPath>>> vMP;

      auto v0 = Vuint(); // ie model not seen in the model path
      for (size_t i = 0; i < nbpath; i++)
      {
        auto smp = make_shared<ModelPath>(*scenario.getModelPath(i));
        if (!smp->hasModel(mixMod))
        {
          if (vMP.find(Vuint()) == vMP.end())
            vMP[v0] = vector<shared_ptr<ModelPath>>(1, smp);
          else
            vMP[v0].push_back(smp);
          continue;
        }

        const Vuint np(smp->getPathNode(mixMod));
        if (vMP.find(np) == vMP.end())
          vMP[np] = vector<shared_ptr<ModelPath>>(1, smp);
        else
          vMP[np].push_back(smp);
      }

      if (vMP.size() > 1) // mixture on edges, build a mixture node
      {
        sonnode->setProperty("event", NodeEvent::mixtureEvent);

        auto vedge = make_shared<ProcessComputationEdge>(nullptr, 0, sonInd, true);
        link(father, sonnode, vedge);
        addEdgeIndex(vedge);

        // and then below sonnode
        for (auto vmp:vMP)
        {
          // edge for proportion
          auto ssonnode = make_shared<ProcessComputationNode>(*son, sonInd);
          createNode(ssonnode);
          addNodeIndex(ssonnode);
          ssonnode->setProperty("event", NodeEvent::speciationEvent);
          auto ssedge = make_shared<ProcessComputationEdge>(mixMod, nmodel, sonInd, true, vmp.first);
          link(sonnode, ssonnode, ssedge);
          addEdgeIndex(ssedge);

          // edge for transition
          auto s3onnode = make_shared<ProcessComputationNode>(*son, sonInd);
          createNode(s3onnode);
          addNodeIndex(s3onnode);
          auto s3edge = make_shared<ProcessComputationEdge>(mixMod, nmodel, sonInd, false, vmp.first);
          link(ssonnode, s3onnode, s3edge);
          addEdgeIndex(s3edge);

          if (vmp.second.size() > 1)
          {
            ModelScenario scen(vmp.second);
            buildFollowingScenario_(s3onnode, scen, mMrca);
          }
          else
            buildFollowingPath_(s3onnode, *vmp.second[0]);
        }
      }
      else // all modelpaths agree on this path, no mixture here
           // (probably done before)
      {
        auto beg = vMP.begin();
        auto nedge = make_shared<ProcessComputationEdge>(mixMod, nmodel, sonInd, false, beg->first);
        link(father, sonnode, nedge);
        addEdgeIndex(nedge);

        if (beg->second.size() == 1) // only one path
          buildFollowingPath_(sonnode, *beg->second[0]);
        else
          buildFollowingScenario_(sonnode, scenario, mMrca);
      }
    }
  }
  else // it is a mixture node
  {
    // build the node
    father->setProperty("event", NodeEvent::mixtureEvent);

    // Check decomposition in the modelpaths
    map<Vuint, vector<shared_ptr<ModelPath>>> vMP;

    auto v0 = Vuint(); // ie model not seen in the model path
    for (size_t i = 0; i < nbpath; i++)
    {
      auto smp = make_shared<ModelPath>(*scenario.getModelPath(i));
      if (!smp->hasModel(mrca))
      {
        if (vMP.find(Vuint()) == vMP.end())
          vMP[v0] = vector<shared_ptr<ModelPath>>(1, smp);
        else
          vMP[v0].push_back(smp);
        continue;
      }

      const Vuint np(smp->getPathNode(mrca));
      if (vMP.find(np) == vMP.end())
        vMP[np] = vector<shared_ptr<ModelPath>>(1, smp);
      else
        vMP[np].push_back(smp);
    }

    size_t nmodel = 0; // number model of the model split at this mrca node
    bool modok(false);

    auto modNbs = process_->getModelNumbers();
    for (auto nb: modNbs)
    {
      if (process_->getModel(nb) == mrca)
      {
        nmodel = nb;
        modok = true;
      }
    }

    if (!modok)
      throw Exception("ProcessComputationTree::ProcessComputationTree : unknown model  for process : " + mrca->getName());


    auto node = ptree->getNode(spInd);

    for (auto vmp:vMP)
    {
      auto sonnode = make_shared<ProcessComputationNode>(*node, spInd);
      createNode(sonnode);
      addNodeIndex(sonnode);

      auto sedge = make_shared<ProcessComputationEdge>(mrca, nmodel, spInd, true, vmp.first);
      link(father, sonnode, sedge);
      addEdgeIndex(sedge);

      if (vmp.second.size() > 1)
      {
        ModelScenario scen(vmp.second);
        buildFollowingScenario_(sonnode, scen, mMrca);
      }
      else
        buildFollowingPath_(sonnode, *vmp.second[0]);
    }
  }
}

void ProcessComputationTree::buildFollowingPath_(shared_ptr<ProcessComputationNode> father, const ModelPath& path)
{
  father->setProperty("event", NodeEvent::speciationEvent);

  auto nodeIndex = father->getSpeciesIndex();

  auto ptree = process_->getParametrizablePhyloTree();

  // and look at the branches

  auto vbrInd = ptree->getBranches(nodeIndex);

  for (const auto brInd:vbrInd)
  {
    // create the son
    auto sonInd = ptree->getSon(brInd);
    auto son = ptree->getNode(sonInd);

    auto sonnode = make_shared<ProcessComputationNode>(*son, sonInd);
    createNode(sonnode);
    addNodeIndex(sonnode);


    // then the edge
    auto modSon = process_->getModelForNode(sonInd);
    auto nmodel = process_->getModelNumberForNode(sonInd);

    shared_ptr<ProcessComputationEdge> nedge;
    auto mixmod = dynamic_pointer_cast<const MixedTransitionModelInterface>(modSon);

    if (mixmod && path.hasModel(mixmod))
      nedge = make_shared<ProcessComputationEdge>(modSon, nmodel, sonInd, false, (const Vuint&)path.getPathNode(mixmod));
    else
      nedge = make_shared<ProcessComputationEdge>(modSon, nmodel, sonInd);

    link(father, sonnode, nedge);
    addEdgeIndex(nedge);

    buildFollowingPath_(sonnode, path);
  }
}
