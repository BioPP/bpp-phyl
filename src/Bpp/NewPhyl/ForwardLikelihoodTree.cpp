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

  auto v = Sequence_DF::create (context_, std::move(initCondLik), sequenceName);
  return v;
}

ForwardLikelihoodBelowRef ForwardLikelihoodTree::makeForwardLikelihoodAtEdge (shared_ptr<ProcessEdge> processEdge, const AlignedValuesContainer & sites)
{
  const auto brlen= processEdge->getBrLen();
  const auto model= processEdge->getModel();
  const auto nMod = processEdge->getNMod();
  const auto brprob = processEdge->getProba();
  
  auto childConditionalLikelihood = makeForwardLikelihoodAtNode (processTree_->getSon(processEdge), sites);

  ForwardLikelihoodBelowRef forwardEdge;
  auto zero=NumericConstant<size_t>::create(context_, 0);
  
  if (brlen) // Branch with transition through a model
  {
    if (dynamic_cast<const TransitionModel*>(model->getTargetValue()))
    {
      auto transitionMatrix = ConfiguredParametrizable::createMatrix<ConfiguredModel, TransitionMatrixFromModel> (context_, {model, brlen, zero, nMod}, transitionMatrixDimension (nbState_));
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
  }

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
      if (mapIndexes_.find(spIndex)==mapIndexes_.end())
        mapIndexes_[spIndex]=DAGindexes();
      mapIndexes_[spIndex].push_back(getNodeIndex(forwardNode));
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
        auto propfath=dynamic_cast<NodeEvent*>(fatherNode->getProperty("event"));
        fathmixture = propfath && propfath->isMixture();
      }
      
      if (!fathmixture)
      {
        if (mapIndexes_.find(spIndex)==mapIndexes_.end())
          mapIndexes_[spIndex]=DAGindexes();
        mapIndexes_[spIndex].push_back(getNodeIndex(forwardNode));
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

  return(forwardNode);
}

// void ForwardLikelihoodTree::setSpeciesMapIndexes_()
// {
//   auto nodeIter = allNodesIterator();
//   while (!nodeIter->end())
//   {
//     auto dagId = getNodeIndex(**nodeIter);
//     auto speciesId = processTree_->getNodeIndex(mapNode_[**nodeIter]);
    
//     if (mapIndexes_.find(speciesId)!=mapIndexes_.end())
//       mapIndexes_[speciesId].push_back(dagId);
//     else
//       mapIndexes_[speciesId]=DAGindexes(1,dagId);

//     nodeIter->next();
//   }
// }
 
