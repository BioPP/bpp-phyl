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

using namespace bpp;
using namespace dataflow;
using namespace std;

ConditionalLikelihoodForwardRef ForwardLikelihoodTree::makeInitialConditionalLikelihood (const string & sequenceName, const AlignedValuesContainer & sites)
{
  size_t nbSites=sites.getNumberOfSites();
  const auto sequenceIndex = sites.getSequencePosition (sequenceName);
  Eigen::MatrixXd initCondLik (nbState_, nbSites);
  for (size_t site = 0; site < nbSites; ++site) {
    for (size_t state = 0; state < nbState_; ++state) {
      initCondLik (Eigen::Index (state), Eigen::Index (site)) =
        sites.getStateValueAt (site, sequenceIndex, statemap_.getAlphabetStateAsInt(state));
    }
  }
  return NumericConstant<Eigen::MatrixXd>::create (context_, move (initCondLik));
}

ForwardLikelihoodBelowRef ForwardLikelihoodTree::makeForwardLikelihoodEdge (shared_ptr<ProcessEdge> processEdge, const AlignedValuesContainer & sites)
{
  cerr << "MAKEFORWARDLIKELIHOODEDGE " << endl;
  const auto brlen= processEdge->getBrLen();
  const auto model= processEdge->getModel();
  const auto nMod = processEdge->getNMod();
  
  auto childConditionalLikelihood = makeConditionalLikelihoodNode (processTree_->getSon(processEdge), sites);
  auto transitionMatrix = ConfiguredParametrizable::createMatrix<ConfiguredModel, TransitionMatrixFromModel> (context_, {model, brlen, nMod}, transitionMatrixDimension (nbState_));
  auto forwardEdge = ForwardLikelihoodFromConditional::create (
    context_, {transitionMatrix, childConditionalLikelihood}, likelihoodMatrixDim_);
  
  if (!hasEdgeIndex(forwardEdge)) // ie this edge does not exist before
  {
    // false (root) top-node of the edge (only way to build an edge in
    // the DAG). Correct top-node will be set afterwards.
    
    link(getRoot(),childConditionalLikelihood,forwardEdge);
    addEdgeIndex(forwardEdge);
    mapEdge_[forwardEdge]=processEdge;
  }
  cerr << "makeforwardlikelihoodedge " << endl;
  return forwardEdge;
}

ConditionalLikelihoodForwardRef ForwardLikelihoodTree::makeConditionalLikelihoodNode (shared_ptr<ProcessNode> processNode, const AlignedValuesContainer & sites)
{
  cerr << "makeConditionalLikelihoodNode " << endl;
  const auto childBranches = processTree_->getBranches (processNode);
  const auto nbChildren = childBranches.size();
  
  ConditionalLikelihoodForwardRef forwardNode;
  
  if (childBranches.empty ())
  {
    forwardNode = makeInitialConditionalLikelihood (processNode->getName (), sites);
    if (!hasNodeIndex(forwardNode)) 
    {
      createNode(forwardNode);
      addNodeIndex(forwardNode);
      mapNode_[forwardNode]=processNode;
    }    
  }
  else {
    auto prop=dynamic_cast<NodeEvent*>(processNode->getProperty("event"));
    if (!prop) 
      throw Exception("ForwardLikelihoodTree::makeConditionalLikelihoodNode : Node has no event associated: Node id " + TextTools::toString(processTree_->getNodeIndex(processNode)));

    if (prop->isSpeciation())
    {
      // depE is used to link ForwardLikelihoodTree edges to the node
      std::vector<ForwardLikelihoodBelowRef> depE(nbChildren);
  
      NodeRefVec deps(nbChildren);
      for (size_t i = 0; i < childBranches.size (); ++i) {
        depE[i] = makeForwardLikelihoodEdge (childBranches[i], sites);
        deps[i] = depE[i];
      }
    
      forwardNode = SpeciationFromChildrenForward::create (context_, std::move(deps),
                                                           likelihoodMatrixDim_);
      if (!hasNodeIndex(forwardNode)) 
      {
        createNode(forwardNode);
        addNodeIndex(forwardNode);
        cerr << "nodeIndex " << getNodeIndex(forwardNode) << "=" << forwardNode << endl;
        mapNode_[forwardNode]=processNode;
    
        for (size_t i = 0; i < childBranches.size (); ++i)
        {
          auto fs=getNodes(depE[i]);
          unlink(fs.first,fs.second);
          link(forwardNode, fs.second, depE[i]);
          cerr << getEdgeLinking(forwardNode, fs.second) << ":" << getNodeIndex(fs.second) << "=" << fs.second << endl;
        }
      }
    }
    else if (prop->isMixture())
    {            
      // depE is used to link ForwardLikelihoodTree edges to the node
      std::vector<ConditionalLikelihoodForwardRef> depN(nbChildren);
      NodeRefVec deps(2*nbChildren);  // for arrays and probas
      for (size_t i = 0; i < nbChildren; ++i) { 
        depN[i] = makeConditionalLikelihoodNode (processTree_->getSon(childBranches[i]), sites);
        deps[i] = depN[i];
      }

      for (size_t i = 0; i < nbChildren; ++i)
      {
        if (!childBranches[i]->getProba())
          throw Exception("ForwardLikelihoodTree::makeConditionalLikelihoodNode: Missing proba for branch " + TextTools::toString(i) + " under mixture node " + TextTools::toString(processTree_->getNodeIndex(processNode)) + ".");
        deps[i+nbChildren] = childBranches[i]->getProba();
      }
      
      forwardNode = MixtureFromChildrenForward::create (context_, std::move(deps),
                                                        likelihoodMatrixDim_);
      if (!hasNodeIndex(forwardNode)) 
      {
        createNode(forwardNode);
        addNodeIndex(forwardNode);
        mapNode_[forwardNode]=processNode;
    
        for (size_t i = 0; i < childBranches.size (); ++i)
          link(forwardNode, depN[i]);  // link with no edge because no time
      }
    }
    else
      throw Exception("ForwardLikelihoodTree::makeConditionalLikelihoodNode : event not recognized.");
  }
  
  cerr << "makeconditionallikelihoodnode " << endl;
  outputToDot("forwardNode_"+TextTools::toString(getNodeIndex(forwardNode))+".dot","forwardTree");
  return(forwardNode);
}

void ForwardLikelihoodTree::setSpeciesMapIndexes_()
{
  auto nodeIter = allNodesIterator();
  while (!nodeIter->end())
  {
    auto dagId = getNodeIndex(**nodeIter);
    auto speciesId = processTree_->getNodeIndex(mapNode_[**nodeIter]);
    
    if (mapIndexes_.find(speciesId)!=mapIndexes_.end())
      mapIndexes_[speciesId].push_back(dagId);
    else
      mapIndexes_[speciesId]=DAGindexes(1,dagId);

    nodeIter->next();
  }
}
 
