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
#include "Parametrizable.h"
#include <Bpp/Phyl/Model/MixedTransitionModel.h>

//From the stl:
#include <string>

using namespace bpp;
using namespace bpp::dataflow;
using namespace std;

ProcessTree::ProcessTree(Context& context,
                         const SubstitutionProcess& process,
                         ParameterList& parList,
                         const BrLenMap& vrefmap) :
  AssociationTreeGlobalGraphObserver<ProcessNode,ProcessEdge>(0), context_(context)
{
  ProcessComputationTree tree(process);
  setGraph(tree.getGraph());
  
  auto vNodes=tree.getAllNodes();
  
  for (const auto& node:vNodes)
  {
    auto pn=make_shared<ProcessNode>(*node);
    associateNode(pn, tree.getNodeGraphid(node));
    setNodeIndex(pn, tree.getNodeIndex(node));
  }
  
  // Build the ConfiguredModels from the TransitionModels

  auto vnMod=process.getModelNumbers();
        
  std::map<size_t, std::shared_ptr<ConfiguredModel>> modelmap;

  for (auto nMod:vnMod)
    modelmap[nMod] = ConfiguredParametrizable::createConfigured<TransitionModel, ConfiguredModel>(context_, *process.getModel(nMod), parList, "_"+ TextTools::toString(nMod));


  // Assign References on all branches
  
  auto vEdges=tree.getAllEdges();

  for (const auto& edge:vEdges)
  {
    auto id=tree.getEdgeGraphid(edge);
    uint spIndex=edge->getSpeciesIndex(); // index of the edge that
                                        // represents the index in the
                                        // ParametrizablePhyloTree
    auto pmodel = modelmap[process.getModelNumberForNode(spIndex)];
    auto vNb=edge->subModelNumbers();

    if (vNb.size()>1)
      throw Exception("ProcessTree::ProcessTree : only simple submodels are used, not combinations.");
    
    std::shared_ptr<ProcessEdge> brref;

    if (!edge->useProb()) // model transition is used
    {
      if (vrefmap.find(spIndex)==vrefmap.end()) // ie there is no branch length
        throw Exception("ProcessTree::ProcessTree missing branch length for branch " + TextTools::toString(spIndex));
      
      if (vNb.size()==0) // ie full model 
        brref=make_shared<ProcessEdge>(spIndex, vrefmap.at(spIndex), pmodel);
      else
      {
        auto nMod=NumericConstant<size_t>::create(context_, vNb[0]); // Only 1st submodel is used        
        brref=make_shared<ProcessEdge>(spIndex, vrefmap.at(spIndex), pmodel, nMod);
      }
    }
    else // a branch issued from a mixture.
    {
      if (vNb.size()==0)
        brref=make_shared<ProcessEdge>(spIndex, ConstantOne<double>::create (context_, Dimension<double>()), pmodel);
      else                         
      {
        auto nMod=NumericConstant<size_t>::create(context_, vNb[0]); // Only 1st submodel is used
        auto brprob=ProbabilityFromMixedModel::create(context_, {pmodel}, vNb[0]);
        brref=make_shared<ProcessEdge>(spIndex, vrefmap.at(spIndex), pmodel, nMod);
      }
    }

    associateEdge(brref, id);
    setEdgeIndex(brref,tree.getEdgeIndex(edge));
  }
}


/*void ProcessTree::buildUnderNode_(const ParametrizablePhyloTree& tree, const BrLenMap& vrefmap, ModelMap& modelmap, shared_ptr<PhyloNode> node, shared_ptr<ProcessNode> newFather, shared_ptr<ProcessEdge> newEdge)
{
  // build a new node as a speciation node, and sets all edges below it
  auto newNode=make_shared<ProcessNode>(*node);
  addNodeIndex(newNode);
  createNode(newNode);
  
  // Specific case of mixture at root
  if (!newFather)
  {
    newNode->setProba(ConstantOne<double>::create(context_, Dimension<double>()));
    if (isRooted())
      rootAt(newNode);

    auto id = tree.getNodeIndex(node);    
    auto vId = tree.getOutgoingNeighbors(id);
    ModelAssign* pmodelass(0);
    
    // root node is a mixture only if the model is the same in the
    // outgoing branches.
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

    // test is mixture
    if (pmodelass)
    {
      auto model=pmodelass->model_;
      const auto assign=pmodelass->modelNum_;
      if (assign.size()>1) // mixture
      {
        newNode->setProperty("event",NodeEvent::mixtureEvent);
        
        for (const auto ass:assign)
        {
          // replace all numbers of model submodels in inheriting
          // branches
          ModelMap modelmap2(modelmap);
          for (auto& id2:modelmap2)
            if (id2.second.model_==model)
              id2.second.modelNum_={ass};

          auto nMod=NumericConstant<size_t>::create(context_, ass);

          auto newMixEdge=make_shared<ProcessEdge>(std::shared_ptr<ConfiguredParameter>(0), model, nMod);
          
          auto proba=ProbabilityFromMixedModel::create(context_, {model}, ass);
          newMixEdge->setProba(proba);

          // build subtree for each submodel 
          buildUnderNode_(tree, vrefmap, modelmap2, node, newNode, newMixEdge);
        }
        return;
      }
    }
  }
  
  // Otherwise not a mixture node (mixture is handled in
  // buildUnderEdgeFromNode_)
  
  newNode->setProperty("event", NodeEvent::speciationEvent);
  if (newFather) // no root node
  {
    link(newFather, newNode, newEdge);
    if (newEdge)
    {
      addEdgeIndex(newEdge);
      newNode->setProba(newEdge->getProba());
    }
    else 
      newNode->setProba(newFather->getProba());
  }
  
  // then the edges below
  auto vEdges=tree.getOutgoingEdges(node);
  for (const auto& oldEdge:vEdges)
    buildUnderEdgeFromNode_(tree, vrefmap, modelmap, oldEdge, newNode, newEdge);
}

void ProcessTree::buildUnderEdgeFromNode_(const ParametrizablePhyloTree& tree, const BrLenMap& vrefmap, ModelMap& modelmap, shared_ptr<PhyloBranchParam> oldEdge, shared_ptr<ProcessNode> newNode, shared_ptr<ProcessEdge> aboveEdge)
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

    // set proba equal to the father node proba
    newEdge->setProba(newNode->getProba());
    
    buildUnderNode_(tree, vrefmap, modelmap, tree.getSon(oldEdge), newNode, newEdge);
  }
  else
  {
    // build a mixture node at distance 0, with edges carrying the submodels
    auto mixNode=make_shared<ProcessNode>(*newNode, tree.getNodeIndex(tree.getFatherOfEdge(oldEdge)));
    mixNode->setProperty("event",NodeEvent::mixtureEvent);
    mixNode->setProba(newEdge->getProba());
    
    createNode(mixNode);
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

      ValueRef<double> proba=ProbabilityFromMixedModel::create(context_, {model}, ass);
      // if a mixture inside a mixture
      if (aboveEdge && aboveEdge->getProba())
      {
        NodeRefVec deps={proba,aboveEdge->getProba()};
        
        proba=CWiseMul<double,std::tuple<double,double>>::create(context_, std::move(deps), Dimension<double>());
      }
      newMixEdge->setProba(proba);
      
      buildUnderNode_(tree, vrefmap, modelmap2, tree.getSon(oldEdge), mixNode, newMixEdge);
    }
  }
}
*/
