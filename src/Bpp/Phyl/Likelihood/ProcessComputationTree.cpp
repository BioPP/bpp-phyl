//
// File: ProcessComputationTree.cpp
// Authors:
//   Laurent GuÃÂ©guen
// Created: mardi 9 juillet 2013, ÃÂ  15h 37
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#include <numeric>

#include "ProcessComputationTree.h"

using namespace bpp;
using namespace std;

ProcessComputationTree::ProcessComputationTree(const SubstitutionProcess& process) :
  BaseTree(process.getModelScenario() ? 0 : process.getParametrizablePhyloTree()->getGraph()),
  process_(process)
{
  std::shared_ptr<const ParametrizablePhyloTree> ptree = process.getParametrizablePhyloTree();
  // if no model scenario, copy the basic tree
  auto scenario = process.getModelScenario();

  if (!scenario)
  {
    auto itN = ptree->allNodesIterator();
    for (itN->start(); !itN->end(); itN->next())
    {
      auto index = ptree->getNodeIndex(*(*itN));
      auto nnode = std::make_shared<ProcessComputationNode>(*(*(*itN)), index);
      nnode->setProperty("event", NodeEvent::speciationEvent);
      associateNode(nnode, ptree->getNodeGraphid(*(*itN)));
      setNodeIndex(nnode, index);
    }

    auto itE = ptree->allEdgesIterator();
    for (itE->start(); !itE->end(); itE->next())
    {
      auto index = ptree->getEdgeIndex(*(*itE));
      auto model = process.getModelForNode(index);
      size_t nmodel = process.getModelNumberForNode(index);
      auto nedge = std::make_shared<ProcessComputationEdge>(model.get(), nmodel, index);
      associateEdge(nedge, ptree->getEdgeGraphid(*(*itE)));
      setEdgeIndex(nedge, index);
    }

    return;
  }

  // Map of the mrca of the MixedTransitionModel split in several paths
  std::map<std::shared_ptr<MixedTransitionModel>, uint> mMrca;

  auto vMod = scenario->getModels();
  std::map<std::shared_ptr<MixedTransitionModel>, std::vector<std::shared_ptr<PhyloNode> > > mnodes;

  auto vNodes = ptree->getAllNodes();

  auto root = ptree->getRoot();

  // first the nodes that carry the models
  for (const auto& node:vNodes)
  {
    if (node == root)
      continue;

    const auto medge = process_.getModelForNode(ptree->getNodeIndex(node));
    std::shared_ptr<MixedTransitionModel> mok(0);
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
      mnodes[mok] = std::vector<std::shared_ptr<PhyloNode> >();
    mnodes[mok].push_back(node);
  }

  // then the mrca

  for (const auto& mnode:mnodes)
  {
    auto nrca = ptree->MRCA(mnode.second);
    mMrca[mnode.first] = ptree->getNodeIndex(nrca);
  }

  // Then construcion of the tree

  auto nroot = std::make_shared<ProcessComputationNode>(*root, ptree->getRootIndex());
  createNode(nroot);
  addNodeIndex(nroot);
  _build_following_scenario(nroot, *scenario, mMrca);

  rootAt(nroot);
}


void ProcessComputationTree::_build_following_scenario(shared_ptr<ProcessComputationNode> father, const ModelScenario& scenario,  std::map<std::shared_ptr<MixedTransitionModel>, uint>& mMrca)
{
  auto spInd = father->getSpeciesIndex();

  size_t nbpath = scenario.getNumberOfModelPaths();

  auto ptree = process_.getParametrizablePhyloTree();

  std::shared_ptr<MixedTransitionModel> mrca(0); // One MixedModel which is
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

      auto modSon = process_.getModelForNode(sonInd);
      size_t nmodel = process_.getModelNumberForNode(sonInd);

      auto sonnode = std::make_shared<ProcessComputationNode>(*son, sonInd);
      createNode(sonnode);
      addNodeIndex(sonnode);

      // Check if model is a mixture
      auto mixMod = dynamic_pointer_cast<const MixedTransitionModel>(modSon);

      if (!mixMod) // No mixture, then a branch with a full model
      {
        auto nedge = std::make_shared<ProcessComputationEdge>(modSon.get(), nmodel, sonInd);
        link(father, sonnode, nedge);
        addEdgeIndex(nedge);
        _build_following_scenario(sonnode, scenario, mMrca);
        continue;
      }

      // If it is a mixture model, check its decomposition in the modelpaths

      map<Vuint, vector<shared_ptr<ModelPath> > > vMP;

      auto v0 = Vuint(); // ie model not seen in the model path
      for (size_t i = 0; i < nbpath; i++)
      {
        auto smp = std::make_shared<ModelPath>(*scenario.getModelPath(i));
        if (!smp->hasModel(mixMod))
        {
          if (vMP.find(Vuint()) == vMP.end())
            vMP[v0] = vector<shared_ptr<ModelPath> >(1, smp);
          else
            vMP[v0].push_back(smp);
          continue;
        }

        const Vuint np(smp->getPathNode(mixMod));
        if (vMP.find(np) == vMP.end())
          vMP[np] = vector<shared_ptr<ModelPath> >(1, smp);
        else
          vMP[np].push_back(smp);
      }

      if (vMP.size() > 1) // mixture on edges, build a mixture node
      {
        sonnode->setProperty("event", NodeEvent::mixtureEvent);

        auto vedge = std::make_shared<ProcessComputationEdge>(nullptr, 0, sonInd, true);
        link(father, sonnode, vedge);
        addEdgeIndex(vedge);

        // and then below sonnode
        for (auto vmp:vMP)
        {
          // edge for proportion
          auto ssonnode = std::make_shared<ProcessComputationNode>(*son, sonInd);
          createNode(ssonnode);
          addNodeIndex(ssonnode);
          ssonnode->setProperty("event", NodeEvent::speciationEvent);
          auto ssedge = std::make_shared<ProcessComputationEdge>(mixMod.get(), nmodel, sonInd, true, vmp.first);
          link(sonnode, ssonnode, ssedge);
          addEdgeIndex(ssedge);

          // edge for transition
          auto s3onnode = std::make_shared<ProcessComputationNode>(*son, sonInd);
          createNode(s3onnode);
          addNodeIndex(s3onnode);
          auto s3edge = std::make_shared<ProcessComputationEdge>(mixMod.get(), nmodel, sonInd, false, vmp.first);
          link(ssonnode, s3onnode, s3edge);
          addEdgeIndex(s3edge);

          if (vmp.second.size() > 1)
          {
            ModelScenario scen(vmp.second);
            _build_following_scenario(s3onnode, scen, mMrca);
          }
          else
            _build_following_path(s3onnode, *vmp.second[0]);
        }
      }
      else // all modelpaths agree on this path, no mixture here
           // (probably done before)
      {
        auto beg = vMP.begin();
        auto nedge = std::make_shared<ProcessComputationEdge>(mixMod.get(), nmodel, sonInd, false, beg->first);
        link(father, sonnode, nedge);
        addEdgeIndex(nedge);

        if (beg->second.size() == 1) // only one path
          _build_following_path(sonnode, *beg->second[0]);
        else
          _build_following_scenario(sonnode, scenario, mMrca);
      }
    }
  }
  else // it is a mixture node
  {
    // build the node
    father->setProperty("event", NodeEvent::mixtureEvent);

    // Check decomposition in the modelpaths
    map<Vuint, vector<shared_ptr<ModelPath> > > vMP;

    auto v0 = Vuint(); // ie model not seen in the model path
    for (size_t i = 0; i < nbpath; i++)
    {
      auto smp = std::make_shared<ModelPath>(*scenario.getModelPath(i));
      if (!smp->hasModel(mrca))
      {
        if (vMP.find(Vuint()) == vMP.end())
          vMP[v0] = vector<shared_ptr<ModelPath> >(1, smp);
        else
          vMP[v0].push_back(smp);
        continue;
      }

      const Vuint np(smp->getPathNode(mrca));
      if (vMP.find(np) == vMP.end())
        vMP[np] = vector<shared_ptr<ModelPath> >(1, smp);
      else
        vMP[np].push_back(smp);
    }

    size_t nmodel; // number model of the model split at this mrca node
    bool modok(false);

    auto modNbs = process_.getModelNumbers();
    for (auto nb: modNbs)
    {
      if (process_.getModel(nb) == mrca)
      {
        nmodel = nb;
        modok = true;
      }
    }

    if (!modok)
      throw Exception("ProcessComputationTree::ProcessComputationTree : unknown model  for process : " + mrca.get()->getName());


    auto node = ptree->getNode(spInd);

    for (auto vmp:vMP)
    {
      auto sonnode = std::make_shared<ProcessComputationNode>(*node, spInd);
      createNode(sonnode);
      addNodeIndex(sonnode);

      auto sedge = std::make_shared<ProcessComputationEdge>(mrca.get(), nmodel, spInd, true, vmp.first);
      link(father, sonnode, sedge);
      addEdgeIndex(sedge);

      if (vmp.second.size() > 1)
      {
        ModelScenario scen(vmp.second);
        _build_following_scenario(sonnode, scen, mMrca);
      }
      else
        _build_following_path(sonnode, *vmp.second[0]);
    }
  }
}

void ProcessComputationTree::_build_following_path(shared_ptr<ProcessComputationNode> father, const ModelPath& path)
{
  father->setProperty("event", NodeEvent::speciationEvent);

  auto nodeIndex = father->getSpeciesIndex();

  auto ptree = process_.getParametrizablePhyloTree();

  // and look at the branches

  auto vbrInd = ptree->getBranches(nodeIndex);

  for (const auto brInd:vbrInd)
  {
    // create the son
    auto sonInd = ptree->getSon(brInd);
    auto son = ptree->getNode(sonInd);

    auto sonnode = std::make_shared<ProcessComputationNode>(*son, sonInd);
    createNode(sonnode);
    addNodeIndex(sonnode);


    // then the edge
    auto modSon = process_.getModelForNode(sonInd);
    auto nmodel = process_.getModelNumberForNode(sonInd);

    shared_ptr<ProcessComputationEdge> nedge;
    auto mixmod = dynamic_pointer_cast<const MixedTransitionModel>(modSon);

    if (mixmod && path.hasModel(mixmod))
      nedge = std::make_shared<ProcessComputationEdge>(modSon.get(), nmodel, sonInd, false, (const Vuint&)path.getPathNode(mixmod));
    else
      nedge = std::make_shared<ProcessComputationEdge>(modSon.get(), nmodel, sonInd);

    link(father, sonnode, nedge);
    addEdgeIndex(nedge);

    _build_following_path(sonnode, path);
  }
}
