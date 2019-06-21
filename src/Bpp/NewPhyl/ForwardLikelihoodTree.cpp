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
  auto r=NumericConstant<Eigen::MatrixXd>::create (context_, move (initCondLik));
  if (!hasNodeIndex(r))
  {
    createNode(r);
    addNodeIndex(r);
  }
  
  return r;
}

ForwardLikelihoodBelowRef ForwardLikelihoodTree::makeForwardLikelihoodNode (shared_ptr<ProcessEdge> edge, const AlignedValuesContainer & sites)
{
  const auto brlen = edge->getBrLen();
  const auto model= edge->getModel();
  const auto nMod = edge->getNMod();
  
  auto childConditionalLikelihood = makeConditionalLikelihoodNode (processTree_->getSon(edge), sites);
  auto transitionMatrix = ConfiguredParametrizable::createMatrix<ConfiguredModel, TransitionMatrixFromModel> (context_, {model, brlen, nMod}, transitionMatrixDimension (nbState_));
  auto r= ForwardLikelihoodFromConditional::create (
    context_, {transitionMatrix, childConditionalLikelihood}, likelihoodMatrixDim_);
  
  if (!hasNodeIndex(r))
  {
    createNode(r);
    addNodeIndex(r);
    link(r,childConditionalLikelihood);
  }
  return r;
}

ConditionalLikelihoodForwardRef ForwardLikelihoodTree::makeConditionalLikelihoodNode (shared_ptr<ProcessNode> node, const AlignedValuesContainer & sites)
{
  const auto childBranches = processTree_->getBranches (node);
  ConditionalLikelihoodForwardRef r;
  NodeRefVec deps(childBranches.size ());
  std::vector<ValueRef<Eigen::MatrixXd>> depM (childBranches.size ());
  
  if (childBranches.empty ()) 
    r=makeInitialConditionalLikelihood (node->getName (), sites);
  else {
    auto prop=dynamic_cast<NodeEvent*>(node->getProperty("event"));
    if (!prop) 
      throw Exception("ForwardLikelihoodTree::makeConditionalLikelihoodNode : Node has no event associated: Node id " + TextTools::toString(processTree_->getNodeIndex(node)));

    if (prop->isSpeciation())
    {
      for (size_t i = 0; i < childBranches.size (); ++i) { 
        depM[i] = makeForwardLikelihoodNode (childBranches[i], sites);
        deps[i] = depM[i];
      }
    
      r= SpeciationFromChildrenForward::create (context_, std::move(deps),
                                                likelihoodMatrixDim_);
    }
    else if (prop->isMixture())
    {            
      for (size_t i = 0; i < childBranches.size (); ++i) { 
        depM[i] = makeConditionalLikelihoodNode (processTree_->getSon(childBranches[i]), sites);
        deps[i] = depM[i];
      }

      auto probas=node->getProba();
      if (!probas)
        throw Exception("ForwardLikelihoodTree::makeConditionalLikelihoodNode : mixture Node has no mixture probabilities associated: Node id " + TextTools::toString(processTree_->getNodeIndex(node)));
      deps.push_back(probas);

      r= MixtureFromChildrenForward::create (context_, std::move(deps),
                                             likelihoodMatrixDim_);
    }
  }
  
  if (!hasNodeIndex(r)) 
  {
    createNode(r);
    addNodeIndex(r);
    for (size_t i = 0; i < childBranches.size (); ++i)
      link(r, depM[i]);      
  }
  return(r);
}

