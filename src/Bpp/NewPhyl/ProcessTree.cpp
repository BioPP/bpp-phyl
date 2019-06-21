//
// File: ProcessTree.cpp
// Created by: Laurent Guéguen
// Created on: mardi 11 juin 2019, à 09h 39
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use, 
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info". 

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability. 

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or 
  data to be ensured and,  more generally, to use and operate it in the 
  same conditions as regards security. 

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>

#include <Bpp/Graph/AssociationTreeGraphImplObserver.h>

#include "ProcessTree.h"
#include "Model.h"
#include "Parameter.h"
#include <Bpp/Phyl/Model/MixedTransitionModel.h>

//From the stl:
#include <string>

using namespace bpp;
using namespace bpp::dataflow;
using namespace std;

ProcessTree::ProcessTree(Context& context, const ParametrizablePhyloTree& tree, const BrLenMap& vrefmap) :
  AssociationTreeGlobalGraphObserver<ProcessNode,ProcessEdge>(tree.getGraph()), context_(context)
{
  vector<uint> vNodesId=tree.getGraph()->getAllNodes();
  
  for (auto& index:vNodesId)
  {
    auto pn=make_shared<ProcessNode>(*tree.getNodeFromGraphid(index));
    pn->setProperty("event",NodeEvent::speciationEvent);
    associateNode(pn, index);
    setNodeIndex(pn, tree.getNodeIndex(tree.getNodeFromGraphid(index)));
  }
  
  // Ids of the branches in the graph may be different from the
  // ids in the phyloTree
  
  vector<uint> vEdgesId=tree.getGraph()->getAllEdges();
  
  for (auto& index:vEdgesId)
  {
    // retrieve PhyloBranch id
    const auto pb=tree.getEdgeFromGraphid(index);
    uint ids=tree.getEdgeIndex(pb);
    
    if (vrefmap.find(ids)!=vrefmap.end())
    {
      auto brref=make_shared<ProcessEdge>(vrefmap.at(ids));
            associateEdge(move(brref), index);
            setEdgeIndex(getEdgeFromGraphid(index),ids);
    }
    else
      throw Exception("ProcessTree::ProcessTree missing reference for branch " + TextTools::toString(ids));
  }
}


void ProcessTree::buildUnderNode_(const ParametrizablePhyloTree& tree, const BrLenMap& vrefmap, ModelMap& modelmap, shared_ptr<PhyloNode> node, shared_ptr<ProcessNode> newFather, shared_ptr<ProcessEdge> newEdge)
{
  // build a new node as a speciation node, and sets all edges below it
  auto newNode=make_shared<ProcessNode>(*node);
  addNodeIndex(newNode);
  createNode(newNode);
  
  if (!newFather)
  {
    if (isRooted())
      rootAt(newNode);
    // check mixture node at root
    auto id = tree.getNodeIndex(node);    
    auto vId=tree.getOutgoingNeighbors(id);
    ModelAssign* pmodelass(0);
    for (auto ids:vId)
    {
      if (modelmap.find(ids)==modelmap.end())
        throw Exception("ProcessTree::buildUnderEdgeFromNode_ missing model on branch " + TextTools::toString(ids));
      if (!pmodelass)
        pmodelass=&modelmap[ids];
      else
        if (pmodelass->model_!=modelmap[ids].model_)
        {
          pmodelass=0;
          break;
        }
    }
    if (pmodelass)
    {
      auto model=pmodelass->model_;
      const auto assign=pmodelass->modelNum_;
      if (assign.size()>1) // mixture
      {
        newNode->setProperty("event",NodeEvent::mixtureEvent);
        auto probas=ProbabilitiesFromMixedModel::create(context_, {model});
        newNode->setProba(probas);
        
        for (const auto ass:assign)
        {
          // replace all numbers of submodel of model submodel in
          // inheriting branches
          ModelMap modelmap2(modelmap);
          for (auto& id2:modelmap2)
            if (id2.second.model_==model)
              id2.second.modelNum_={ass};

          auto nMod=NumericConstant<size_t>::create(context_, ass);
          auto newMixEdge=make_shared<ProcessEdge>(std::shared_ptr<ConfiguredParameter>(0), model, nMod);
          buildUnderNode_(tree, vrefmap, modelmap2, node, newNode, newMixEdge);
        }
        return;
      }
    }
  }
  
  // if not mixture node
  newNode->setProperty("event",NodeEvent::speciationEvent);
  if (newFather) // no root node
  {
    link(newFather, newNode, newEdge);
    if (newEdge)
      addEdgeIndex(newEdge);
  }
  
  // then the edges below
  auto vEdges=tree.getOutgoingEdges(node);
  for (const auto& oldEdge:vEdges)
    buildUnderEdgeFromNode_(tree, vrefmap, modelmap, oldEdge, newNode);
}

void ProcessTree::buildUnderEdgeFromNode_(const ParametrizablePhyloTree& tree, const BrLenMap& vrefmap, ModelMap& modelmap, shared_ptr<PhyloBranchParam> oldEdge, shared_ptr<ProcessNode> newNode)
{
  uint id=tree.getEdgeIndex(oldEdge);
  
  if (modelmap.find(id)==modelmap.end())
    throw Exception("ProcessTree::buildUnderEdgeFromNode_ missing model on branch " + TextTools::toString(id));
  
  auto model=modelmap[id].model_;
  const auto assign=modelmap[id].modelNum_;
  auto newEdge=make_shared<ProcessEdge>();
  
  if (assign.size()<=1)
  {
    // speciation node, a length is required
    newEdge->setModel(model);
    newEdge->setBrLen(vrefmap.at(id));
    if (assign.size()==1)
    {
      auto nMod=NumericConstant<size_t>::create(context_, assign[0]);
      newEdge->setNMod(nMod);
    }
    
    buildUnderNode_(tree, vrefmap, modelmap, tree.getSon(oldEdge), newNode, newEdge);
  }
  else
  {
    // build a mixture node at distance 0, with edges carrying the submodels
    auto mixNode=make_shared<ProcessNode>(*newNode);
    mixNode->setProperty("event",NodeEvent::mixtureEvent);
    createNode(mixNode);
    auto probas=ProbabilitiesFromMixedModel::create(context_, {model});
    mixNode->setProba(probas);
    addNodeIndex(mixNode);
    link(newNode, mixNode, newEdge);
    addEdgeIndex(newEdge);
    
    for (const auto ass:assign)
    {
      // replace all numbers of submodel of model submodel in
      // inheriting branches
      ModelMap modelmap2(modelmap);
      for (auto& id2:modelmap2)
        if (id2.second.model_==model)
          id2.second.modelNum_={ass};

      auto nMod=NumericConstant<size_t>::create(context_, ass);
      auto newMixEdge=make_shared<ProcessEdge>(vrefmap.at(id), model, nMod);
      
      buildUnderNode_(tree, vrefmap, modelmap2, tree.getSon(oldEdge), mixNode, newMixEdge);
    }
  }
}

