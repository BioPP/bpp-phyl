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
using namespace std;

ProcessTree::ProcessTree(Context& context,
                         const SubstitutionProcess& process,
                         const ProcessComputationTree& tree,
                         ParameterList& parList,
                         const BrLenMap& vrefmap) :
  AssociationTreeGlobalGraphObserver<ProcessNode,ProcessEdge>(tree.getGraph()), context_(context)
{
  auto vNodes=tree.getAllNodes();
  
  for (const auto& node:vNodes)
  {
    auto pn=make_shared<ProcessNode>(*node);
    associateNode(pn, tree.getNodeGraphid(node));
    setNodeIndex(pn, tree.getNodeIndex(node));
  }

  // Build the ConfiguredModels from the BranchModels

  auto vnMod=process.getModelNumbers();
        
  std::map<const BranchModel*, std::shared_ptr<ConfiguredModel>> modelmap;

  for (auto nMod:vnMod)
  {
    auto mod=process.getModel(nMod);
    modelmap[mod] = ConfiguredParametrizable::createConfigured<BranchModel, ConfiguredModel>(context_, *mod, parList, (nMod==1?"":"_"+ TextTools::toString(nMod))); // suffix "_1" will be added if necessary 
  }

  // Assign References on all branches

  auto vEdges=tree.getAllEdges();

  const auto partree = process.getParametrizablePhyloTree();

  for (const auto& edge:vEdges)
  {
    std::shared_ptr<ProcessEdge> brref;

    auto id=tree.getEdgeGraphid(edge);
    uint spIndex=edge->getSpeciesIndex(); // index of the matching
                                          // edge in the
                                        // ParametrizablePhyloTree

    auto model=edge->getModel();
    if (!model) // ie empty branch
    {
      brref=make_shared<ProcessEdge>(spIndex, nullptr);
      associateEdge(brref, id);
      setEdgeIndex(brref,tree.getEdgeIndex(edge));
      continue;
    }
    
    auto modelit=modelmap.find(model);
    if (modelit==modelmap.end())
      throw Exception("ProcessTree::ProcessTree : Model unknown " + model->getName() + "  for node " + TextTools::toString(spIndex));
    
    std::shared_ptr<ConfiguredModel> pmodel= modelit->second;
    
    auto vNb=edge->subModelNumbers();

    if (vNb.size()>1)
      throw Exception("ProcessTree::ProcessTree : only simple submodels are used, not combinations. Ask developpers");
    
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
        brref=make_shared<ProcessEdge>(spIndex, ConstantOne<double>::create (context_, Dimension<double>()));
      else                         
      {
        auto nMod=NumericConstant<size_t>::create(context_, (size_t)vNb[0]); // Only 1st submodel is used
        auto brprob=ProbabilityFromMixedModel::create(context_, {pmodel}, (size_t)vNb[0]);
        brref=make_shared<ProcessEdge>(spIndex, brprob);
      }
    }

    associateEdge(brref, id);
    setEdgeIndex(brref,tree.getEdgeIndex(edge));
  }

}

