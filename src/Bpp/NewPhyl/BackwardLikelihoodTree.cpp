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
using namespace bpp::dataflow;
using namespace std;


BackwardLikelihoodAboveRef BackwardLikelihoodTree::makeBackwardAboveLikelihoodEdge (PhyloTree::EdgeIndex edgeIndex)
{
  cerr << "makeBackwardAboveLikelihoodEdge " << edgeIndex << endl;
  if (!forwardTree_)
    throw Exception("BackwardLikelihoodTree::makeBackwardAboveLikelihoodEdge: forwardTree_ is missing.");
  
  auto fatherIndex = forwardTree_->getFatherOfEdge(edgeIndex);
  auto fatherForward = forwardTree_->getNode(fatherIndex);

  // get/check if edge with backward likelihood exists
  ValueRef<Eigen::MatrixXd> backNode = hasNode(fatherIndex)
    ? getNode(fatherIndex)
    : makeConditionalAboveLikelihoodNode(fatherIndex);

  // get forward likelihoods of brothers
  auto edgeIds = forwardTree_->getOutgoingEdges(fatherIndex);

  NodeRefVec deps;
  for (auto eId : edgeIds)
  {
    if (eId!=edgeIndex)
      deps.push_back(forwardTree_->getEdge(eId));
  }
  deps.push_back(backNode);

  auto backwardEdge = ConditionalLikelihoodFromUpper::create (context_, std::move (deps),
                                                              likelihoodMatrixDim_);  

  associateEdge(backwardEdge, forwardTree_->getEdgeGraphid(forwardTree_->getEdge(edgeIndex)));
  setEdgeIndex(backwardEdge, edgeIndex);
  writeGraphToDot("backwardEdge_"+TextTools::toString(edgeIndex)+".dot",{backwardEdge.get()});
  cerr << "makebackwardabovelikelihoodedge " << edgeIndex << endl;
  return backwardEdge;
}


ConditionalLikelihoodRef BackwardLikelihoodTree::makeConditionalAboveLikelihoodNode (PhyloTree::NodeIndex nodeIndex)
{
  cerr << "makeConditionalAboveLikelihoodNode "  << nodeIndex << endl;
  if (!forwardTree_)
    throw Exception("BackwardLikelihoodTree::makeBackwardAboveLikelihoodEdge: forwardTree_ is missing.");

  auto forwardNode = forwardTree_->getNode(nodeIndex);
  cerr << nodeIndex << "=" << forwardNode << endl;
  // if root
  if (nodeIndex==forwardTree_->getRootIndex())
    return setRootFrequencies(rFreqs_);
    
  // else get fathers
  const auto fathersIndexes = forwardTree_->getIncomingNeighbors(nodeIndex);

  // get upper dependencies
  NodeRefVec deps;
  NodeRefVec probaDeps;
  ConditionalLikelihoodRef backwardNode(0);

  for (const auto& fatherIndex:fathersIndexes)
  {
    auto forwardFather = forwardTree_->getNode(fatherIndex);
    cerr << "node: " << forwardTree_->nodeToString(forwardFather) << endl;
    auto forwardEdgeToFather = forwardTree_->getEdgeLinking(forwardFather, forwardNode);

    cerr << "edge: " << forwardTree_->edgeToString(forwardEdgeToFather) << endl;
    // if an edge exists on forward tree (ie a model)
    if (forwardEdgeToFather)
    {
      const auto processEdge = forwardTree_->getProcessEdge(forwardEdgeToFather); 
      const auto brlen = processEdge->getBrLen();
      const auto model = processEdge->getModel();
      probaDeps.push_back(processEdge->getProba());
      
      // useless, the transitionMatrix already exists, but is not
      // available directly through a tree
      // ToDo : find another way to get it
      auto transitionMatrix =
        ConfiguredParametrizable::createMatrix<ConfiguredModel, TransitionMatrixFromModel> (context_, {model, brlen}, transitionMatrixDimension (nbState_));

      
      // build the backward top edge to this father if it does not exist
      const auto edgeToFatherIndex = forwardTree_->getEdgeIndex(forwardEdgeToFather);
      BackwardLikelihoodAboveRef backEdge = hasEdge(edgeToFatherIndex)
        ? getEdge(edgeToFatherIndex)
        : makeBackwardAboveLikelihoodEdge(edgeToFatherIndex);

      // and uses the transposed transition matrix to compute the
      // bottom of the edge
      auto backLikeEdge=BackwardLikelihoodFromConditional::create (
        context_, {transitionMatrix, backEdge}, likelihoodMatrixDim_);

      if (fathersIndexes.size()==1)
      {
        backwardNode = backLikeEdge;
        break;
      }
      else
        deps.push_back(backLikeEdge);
    }
    else
    {
      if (deps.size()!=0)
        throw Exception("BackwardLikelihoodFromConditional::makeConditionalAboveLikelihoodNode mixing edges and non-edges.");
      backwardNode = hasNode(fatherIndex)
        ? getNode(fatherIndex)
        : makeConditionalAboveLikelihoodNode(fatherIndex);
    }
  }
  // Then mean if several incoming nodes
  if (deps.size()>1)
  {
    for (auto& probaDep:probaDeps)
      deps.push_back(probaDep);
    backwardNode = MixtureFromUpperBackward::create(context_, std::move(deps), likelihoodMatrixDim_);
  }

  if (!hasNode(backwardNode))
  {
    associateNode(backwardNode, forwardTree_->getNodeGraphid(forwardTree_->getNode(nodeIndex)));
    setNodeIndex(backwardNode, nodeIndex);
  }

  writeGraphToDot("backwardNode_"+TextTools::toString(nodeIndex)+".dot",{backwardNode.get()});
  cerr << "makeconditionalabovelikelihoodnode " << nodeIndex << endl;
  return backwardNode;
}


