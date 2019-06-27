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


dataflow::ValueRef<Eigen::MatrixXd> BackwardLikelihoodTree::makeBackwardAboveLikelihoodEdge (PhyloTree::EdgeIndex index) {
  
  auto fatherIndex = processTree_->getFatherOfEdge(index);

  // get/check if edge with backward likelihood exists
  dataflow::ValueRef<Eigen::MatrixXd> backNode;
  
  if (hasNode(fatherIndex))
    backNode = getNode(fatherIndex);
  else
    backNode = makeConditionalAboveLikelihoodNode(fatherIndex);
  
  // get forward likelihoods of brothers
  if (!forwardTree_)
    throw Exception("BackwardLikelihoodTree::makeBackwardAboveLikelihoodEdge: forwardTree_ is missing.");
  
  auto edgeIds = forwardTree_->getOutgoingEdges(fatherIndex);
  dataflow::NodeRefVec deps;
  
  for (auto eId : edgeIds)
  {
    if (eId!=index)
      deps.push_back(forwardTree_->getEdge(eId));
  }
  deps.push_back(backNode);
  
  auto r= dataflow::ConditionalLikelihoodFromBrothersBackward::create (context_, std::move (deps),
                                                                       likelihoodMatrixDim_);
  
  associateEdge(r,index);
  setEdgeIndex(r, index);
  return r;
}

dataflow::ValueRef<Eigen::MatrixXd> BackwardLikelihoodTree::makeConditionalAboveLikelihoodNode (PhyloTree::NodeIndex index) {
  //!!!! must be initialized before
  if (index==processTree_->getRootIndex())
  {
    setRootFrequencies(rFreqs_);
    return getNode(index);
  }
  
  // get Edge with model
  const auto edgeToFather = processTree_->getEdgeToFather(index);
//  const auto edgeToFatherIndex = processTree_->getEdgeIndex(edgeToFather);
  
  // get/check if edge with backward likelihood exists
  dataflow::ValueRef<Eigen::MatrixXd> backEdge;
  
  // if (hasEdge(edgeToFatherIndex))
  //   backEdge = getEdgeToFather(index);
  // else
  //   backEdge=makeBackwardAboveLikelihoodEdge(index);
  
  const auto brlen = edgeToFather->getBrLen();
  const auto model= edgeToFather->getModel();
  // useless, the transitionMatrix already exists, but is not available directly through a tree
  auto transitionMatrix =
    dataflow::ConfiguredParametrizable::createMatrix<dataflow::ConfiguredModel, dataflow::TransitionMatrixFromModel> (context_, {model, brlen}, transitionMatrixDimension (nbState_));
  
  auto r= dataflow::BackwardLikelihoodFromConditional::create (
    context_, {transitionMatrix, backEdge}, likelihoodMatrixDim_);
  
  associateNode(r, index);
  setNodeIndex(r, index);
  return r;
}


