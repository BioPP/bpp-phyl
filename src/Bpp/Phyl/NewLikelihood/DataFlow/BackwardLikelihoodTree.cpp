//
// File: BackwardLikelihoodTree.cpp
// Created by: Laurent Guéguen
// Created on: vendredi 21 juin 2019, à 18h 41
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

#include "BackwardLikelihoodTree.h"
#include "Parametrizable.h"

using namespace bpp;
using namespace std;

BackwardLikelihoodAboveRef BackwardLikelihoodTree::makeBackwardLikelihoodAtEdge (PhyloTree::EdgeIndex edgeIndex)
{
  if (!forwardTree_)
    throw Exception("BackwardLikelihoodTree::makeBackwardLikelihoodAtEdge: forwardTree_ is missing.");

  auto edgeForward = forwardTree_->getEdge(edgeIndex);
  auto fatherIndex = forwardTree_->getFatherOfEdge(edgeIndex);
  auto fatherForward = forwardTree_->getNode(fatherIndex);

  // get/check if node with backward likelihood exists
  auto backNode = hasNode(fatherIndex)
    ? getNode(fatherIndex)
    : makeBackwardLikelihoodAtNode(fatherIndex);

  auto fatherNode = forwardTree_->getProcessTree()->getNode(fatherIndex);

  BackwardLikelihoodAboveRef backwardEdge;
  
  if (fatherNode->isSpeciation())
  {
    // get forward likelihoods of brothers
    auto edgeIds = forwardTree_->getOutgoingEdges(fatherIndex);

    NodeRefVec deps;
    for (auto eId : edgeIds)
    {
      if (eId!=edgeIndex)
        deps.push_back(forwardTree_->getEdge(eId));
    }
    deps.push_back(backNode);
    
    backwardEdge = SpeciationBackward::create (context_, std::move (deps), likelihoodMatrixDim_);  
  }
  else if (fatherNode->isMixture()) // Transmit array with no modification
  {
    backwardEdge = backNode;
  }
  else
    throw Exception("BackwardLikelihoodTree::makeBackwardLikelihoodAtEdge : event not recognized for node " + TextTools::toString(fatherNode->getSpeciesIndex()));

  // put object in the tree
  if (!hasEdge(backwardEdge))
  {
    associateEdge(backwardEdge, forwardTree_->getEdgeGraphid(edgeForward));
    setEdgeIndex(backwardEdge, edgeIndex);
  }

  return backwardEdge;
}


ConditionalLikelihoodRef BackwardLikelihoodTree::makeBackwardLikelihoodAtNode (PhyloTree::NodeIndex nodeIndex)
{
  if (!forwardTree_)
    throw Exception("BackwardLikelihoodTree::makeBackwardLikelihoodAtNode: forwardTree_ is missing.");

  auto forwardNode = forwardTree_->getNode(nodeIndex);
  // if root
  if (nodeIndex==forwardTree_->getRootIndex())
    return setRootFrequencies(rFreqs_);
    
  // else get incoming edges
  const auto edgesIndexes = forwardTree_->getIncomingEdges(nodeIndex);

  // get upper dependencies
  NodeRefVec deps;
  ConditionalLikelihoodRef backwardNode(0);

  for (const auto& edgeIndex:edgesIndexes)
  {
    auto backEdge = hasEdge(edgeIndex)
      ? getEdge(edgeIndex)
      : makeBackwardLikelihoodAtEdge(edgeIndex);

    auto iedge = forwardTree_->getEdge(edgeIndex);
    
    const auto processEdge = forwardTree_->getProcessTree()->getEdge(edgeIndex);
    
    const auto brlen= processEdge->getBrLen();
    const auto model= processEdge->getModel();
    const auto nMod = processEdge->getNMod();
    const auto brprob = processEdge->getProba();

    ConditionalLikelihoodRef backLikeEdge;
    
    if (brlen) // Branch with transition through a model
    {      
      // useless, the transitionMatrix already exists, but is not
      // available directly through a tree. This new object will be
      // deleted since it already exists in the context
      //
      // ToDo : find another way to get it

      auto zero=NumericConstant<size_t>::create(context_, 0);

      auto transitionMatrix =
        ConfiguredParametrizable::createMatrix<ConfiguredModel, TransitionMatrixFromModel> (context_, {model, brlen, zero, nMod}, transitionMatrixDimension (nbState_));

      // Uses the transposed transition matrix to compute the bottom
      // of the edge

      backLikeEdge=BackwardTransition::create (context_, {transitionMatrix, backEdge}, likelihoodMatrixDim_);
    }
    else if (brprob)
      backLikeEdge = BackwardProportion::create(context_, {brprob, backEdge}, likelihoodMatrixDim_);
    else
      throw Exception("BackwardLikelihoodTree::makeBackwardLikelihoodAtNode : missing information on edge " + processEdge->getSpeciesIndex());
    
    if (edgesIndexes.size()==1)
    {
      backwardNode = backLikeEdge;
      break;
    }
    else
      deps.push_back(backLikeEdge);
  }

  // Then compute sum if several incoming nodes
  if (deps.size()>1)
    backwardNode = MixtureBackward::create(context_, std::move(deps), likelihoodMatrixDim_);

  
  if (!hasNode(backwardNode))
  {
    associateNode(backwardNode, forwardTree_->getNodeGraphid(forwardNode));
    setNodeIndex(backwardNode, nodeIndex);
  }

  return backwardNode;
}


