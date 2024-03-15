// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
      if (eId != edgeIndex)
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
  if (nodeIndex == forwardTree_->getRootIndex())
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

    const auto brprob = processEdge->getProba();

    ConditionalLikelihoodRef backLikeEdge;

    if (processEdge->getBrLen()) // Branch with transition through a model
    {
      auto transitionMatrix = processEdge->getTransitionMatrix();

      // Uses the transposed transition matrix to compute the bottom
      // of the edge

      backLikeEdge = BackwardTransition::create (context_, {transitionMatrix, backEdge}, likelihoodMatrixDim_);
    }
    else if (brprob)
      backLikeEdge = BackwardProportion::create(context_, {brprob, backEdge}, likelihoodMatrixDim_);
    else
      throw Exception("BackwardLikelihoodTree::makeBackwardLikelihoodAtNode : missing information on edge " + TextTools::toString(processEdge->getSpeciesIndex()));

    if (edgesIndexes.size() == 1)
    {
      backwardNode = backLikeEdge;
      break;
    }
    else
      deps.push_back(backLikeEdge);
  }

  // Then compute sum if several incoming nodes
  if (deps.size() > 1)
    backwardNode = MixtureBackward::create(context_, std::move(deps), likelihoodMatrixDim_);


  if (!hasNode(backwardNode))
  {
    associateNode(backwardNode, forwardTree_->getNodeGraphid(forwardNode));
    setNodeIndex(backwardNode, nodeIndex);
  }

  return backwardNode;
}
