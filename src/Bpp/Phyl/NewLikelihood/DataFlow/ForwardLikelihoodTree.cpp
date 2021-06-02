//
// File: ForwardLikelihoodTree.cpp
// Created by: Laurent Guéguen
// Created on: mardi 11 juin 2019, à 10h 09
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

#include "ForwardLikelihoodTree.h"

#include "Model.h"
#include "Parametrizable.h"

#include "Sequence_DF.h"

using namespace bpp;
using namespace std;

ConditionalLikelihoodForwardRef ForwardLikelihoodTree::makeInitialConditionalLikelihood (const string & sequenceName, const AlignedValuesContainer & sites)
{
  size_t nbSites=sites.getNumberOfSites();
  const auto sequenceIndex = sites.getSequencePosition (sequenceName);
  Eigen::MatrixXd initCondLik (nbState_, nbSites);
  for (size_t site = 0; site < nbSites; ++site) {
    for (auto state = 0; state < nbState_; ++state) {
      initCondLik (Eigen::Index (state), Eigen::Index (site)) =
        sites (site, sequenceIndex, statemap_.getAlphabetStateAsInt(size_t(state)));
    }
  }

  return Sequence_DF::create (context_, std::move(initCondLik), sequenceName);
}

ForwardLikelihoodBelowRef ForwardLikelihoodTree::makeForwardLikelihoodAtEdge (shared_ptr<ProcessEdge> processEdge, const AlignedValuesContainer & sites)
{
  const auto brlen= processEdge->getBrLen();
  const auto model= processEdge->getModel();
  const auto nMod = processEdge->getNMod();
  const auto brprob = processEdge->getProba();
  
  auto childConditionalLikelihood = makeForwardLikelihoodAtNode (processTree_->getSon(processEdge), sites);

  ForwardLikelihoodBelowRef forwardEdge;

  auto zero=NumericConstant<size_t>::create(context_, size_t(0));

  if (brlen) // Branch with transition through a model
  {
    if (dynamic_cast<const TransitionModel*>(model->getTargetValue()))
    {
      auto transitionMatrix = ConfiguredParametrizable::createMatrix<ConfiguredModel, TransitionMatrixFromModel> (context_, {model, brlen, zero, nMod}, transitionMatrixDimension (size_t(nbState_)));
      
      processEdge->setTransitionMatrix(transitionMatrix);
      forwardEdge = ForwardTransition::create (
        context_, {transitionMatrix, childConditionalLikelihood}, likelihoodMatrixDim_);
    }
    else{
      auto transitionFunction = TransitionFunctionFromModel::create(context_, {model, brlen, zero}, transitionFunctionDimension(nbState_));
      
      forwardEdge = ForwardTransitionFunction::create(context_ , {childConditionalLikelihood, transitionFunction}, likelihoodMatrixDim_);
    }
  }
  else if (brprob)
  {
    forwardEdge = ForwardProportion::create(
      context_, {brprob, childConditionalLikelihood}, likelihoodMatrixDim_);
  }
  else // junction branch above a mixture node
  {
    forwardEdge = childConditionalLikelihood;
  }

  
  if (!hasEdgeIndex(forwardEdge)) // ie this edge does not exist before
  {
    // false (root) top-node of the edge (only way to build an edge in
    // the DAG). Correct top-node will be set afterwards.

    link(getRoot(),childConditionalLikelihood,forwardEdge);
    setEdgeIndex(forwardEdge, processTree_->getEdgeIndex(processEdge)); // gets the index of the corresponding branch in the processTree_

    auto spIndex = processEdge->getSpeciesIndex();
    
    if (brprob || brlen)
    {
      if (mapEdgesIndexes_.find(spIndex)==mapEdgesIndexes_.end())
        mapEdgesIndexes_[spIndex]=DAGindexes();
      mapEdgesIndexes_[spIndex].push_back(getEdgeIndex(forwardEdge));
    }
  }

#ifdef DEBUG
  std::cerr << "E " << getEdgeIndex(forwardEdge) << " : forwardEdge " << forwardEdge << endl;
  std::cerr << "   -> N " << processTree_->getNodeIndex(processTree_->getSon(processEdge)) << std::endl;
#endif
  
  return forwardEdge;
}

ConditionalLikelihoodForwardRef ForwardLikelihoodTree::makeForwardLikelihoodAtNode (shared_ptr<ProcessNode> processNode, const AlignedValuesContainer & sites)
{
  const auto childBranches = processTree_->getBranches (processNode);

  auto spIndex=processNode->getSpeciesIndex();

  ConditionalLikelihoodForwardRef forwardNode;

  if (childBranches.empty ())
  {
    forwardNode = makeInitialConditionalLikelihood (processNode->getName (), sites);
    if (!hasNodeIndex(forwardNode)) 
    {
      createNode(forwardNode);
      setNodeIndex(forwardNode, processTree_->getNodeIndex(processNode));
      if (mapNodesIndexes_.find(spIndex)==mapNodesIndexes_.end())
        mapNodesIndexes_[spIndex]=DAGindexes();
      mapNodesIndexes_[spIndex].push_back(getNodeIndex(forwardNode));
    }
  }
  else {
    // depE are edges used to link ForwardLikelihoodTree edges to this
    // node
    std::vector<ForwardLikelihoodBelowRef> depE(childBranches.size());
    NodeRefVec deps(childBranches.size());
    
    for (size_t i = 0; i < childBranches.size (); ++i) {
      depE[i]=makeForwardLikelihoodAtEdge (childBranches[i], sites);
      deps[i]=depE[i];
    }

    if (processNode->isSpeciation())
      forwardNode = SpeciationForward::create(context_, std::move(deps),
                                                          likelihoodMatrixDim_);
    else if (processNode->isMixture())
      forwardNode = MixtureForward::create(context_, std::move(deps),
                                                       likelihoodMatrixDim_);
    else
      throw Exception("ForwardLikelihoodTree::makeConditionalLikelihoodAtNode : event not recognized for node " + TextTools::toString(processNode->getSpeciesIndex()));

    // Fix the DAG 
    if (!hasNodeIndex(forwardNode)) 
    {
      createNode(forwardNode);
      setNodeIndex(forwardNode,processTree_->getNodeIndex(processNode));

      // Put the node in speciesIndex map if it is not son of a
      // mixture node (which would hold the species index)

      bool fathmixture(false);

      if (processTree_->hasFather(processNode))
      {
        auto fatherNode = processTree_->getFatherOfNode (processNode);
        fathmixture = fatherNode->isMixture();
      }
      
      if (!fathmixture)
      {
        if (mapNodesIndexes_.find(spIndex)==mapNodesIndexes_.end())
          mapNodesIndexes_[spIndex]=DAGindexes();
        mapNodesIndexes_[spIndex].push_back(getNodeIndex(forwardNode));
      }

      for (size_t i = 0; i < depE.size (); ++i)
      {
        auto fs=getNodes(depE[i]);
        // fix the top node of the edge (because it was set to bidon
        unlink(fs.first,fs.second);
        link(forwardNode, fs.second, depE[i]);
      }
    }
    
  }

#ifdef DEBUG
  std::cerr << "N " << getNodeIndex(forwardNode) << " : forwardNode " << forwardNode << endl;
  for (size_t i = 0; i < childBranches.size (); ++i) {
    std::cerr << "  -> E " << processTree_->getEdgeIndex(childBranches[i]) << std::endl;
  }

#endif
  
  return(forwardNode);
}

/************************************************************/
/************************************************************
 * For ProbabilityDAG
 */

ProbabilityDAG::ProbabilityDAG(std::shared_ptr<ForwardLikelihoodTree> forwardTree) :
  DAProb(forwardTree->getGraph()), context_(forwardTree->getContext())
{
  auto rootProb = ConstantOne<double>::create(context_, Dimension<double>());

  associateNode(rootProb, forwardTree->getNodeGraphid(forwardTree->getRoot()));
  setNodeIndex(rootProb, forwardTree->getRootIndex());
  
  auto allIndex = forwardTree->getAllNodesIndexes();
  std::vector<const Node_DF*> vN;

  for (auto id: allIndex)
  {
    // Start from external nodes (which may be not leaves)
    if (forwardTree->getOutgoingEdges(id).size()==0) 
    {
      auto n=makeProbaAtNode_(id, forwardTree).get();
      n->getTargetValue();
      vN.push_back(n);
    }
  }

  // using bpp::DotOptions;
  // writeGraphToDot("proba.dot", vN, DotOptions::DetailedNodeInfo);

}


ProbaRef ProbabilityDAG::makeProbaAtEdge_ (PhyloTree::EdgeIndex edgeIndex, std::shared_ptr<ForwardLikelihoodTree> forwardTree)
{
  auto fatherIndex = forwardTree->getFatherOfEdge(edgeIndex);
  auto edgeForward = forwardTree->getEdge(edgeIndex);

  // get/check if node with backward likelihood exists
  auto probaNode = hasNode(fatherIndex)
    ? getNode(fatherIndex)
    : makeProbaAtNode_(fatherIndex, forwardTree);

  const auto processEdge = forwardTree->getProcessTree()->getEdge(edgeIndex);
  
  const auto brprob = processEdge->getProba();

  ProbaRef probaEdge;

  if (brprob)
  {
    probaEdge = ProbaMul::create(context_, {brprob, probaNode}, Dimension<double>());
  }
  else
  {
    probaEdge = Identity<double>::create(context_, {probaNode}, Dimension<double>(), edgeIndex);
  }

  // put object in the tree
  if (!hasEdge(probaEdge))
  {
    associateEdge(probaEdge, forwardTree->getEdgeGraphid(edgeForward));
    setEdgeIndex(probaEdge, edgeIndex);
  }

  return probaEdge;
}


ProbaRef ProbabilityDAG::makeProbaAtNode_ (PhyloTree::NodeIndex nodeIndex, std::shared_ptr<ForwardLikelihoodTree> forwardTree)
{
  // else get incoming edges
  const auto edgesIndexes = forwardTree->getIncomingEdges(nodeIndex);

  // get upper dependencies
  NodeRefVec deps;

  ProbaRef probaNode;
  auto forwardNode = forwardTree->getNode(nodeIndex);
  
  for (const auto& edgeIndex:edgesIndexes)
  {
    auto backEdge = hasEdge(edgeIndex)
      ? getEdge(edgeIndex)
      : makeProbaAtEdge_(edgeIndex, forwardTree);

    if (edgesIndexes.size()==1)
    {
      probaNode = Identity<double>::create(context_, {backEdge}, Dimension<double>(), nodeIndex);
      break;
    }
    else
      deps.push_back(backEdge);
  }

  // Then compute sum if several incoming nodes
  if (deps.size()>1)
    probaNode = ProbaSum::create(context_, std::move(deps), Dimension<double>());

  if (!hasNode(probaNode))
  {
    associateNode(probaNode, forwardTree->getNodeGraphid(forwardNode));
    setNodeIndex(probaNode, nodeIndex);
  }

  return probaNode;
}


