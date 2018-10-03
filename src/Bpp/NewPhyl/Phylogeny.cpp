//
// File: Phylogeny.cpp
// Authors:
// Created: 2017-06-06
// Last modified: 2017-06-06
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/NewPhyl/Likelihood.h>
#include <Bpp/NewPhyl/LinearAlgebra.h> // allow conversion to nodeRef FIXME
#include <Bpp/NewPhyl/Model.h>
#include <Bpp/NewPhyl/Phylogeny.h>
#include <Bpp/NewPhyl/PhylogenyTypes.h>
#include <Bpp/NewPhyl/Utils.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <utility>

namespace bpp {

std::size_t SequenceNodesInilialisedFromNames::getNbSites () const {
	return sequences_.getNumberOfSites ();
}

// DF graph construction functions

LikelihoodDataDimension likelihoodDataDimension (const SequenceNodeAccess & sequences,
                                                 const ModelNodeAccess & models) {
	return {sequences.getNbSites (), models.getNbStates ()};
}

DF::ValueRef<MatrixDouble> makeConditionalLikelihoodNode (const TreeTopologyInterface & tree,
                                                          const SequenceNodeAccess & sequenceNodes,
                                                          const BranchLengthNodeAccess & brlenNodes,
                                                          const ModelNodeAccess & modelNodes,
                                                          TopologyNodeIndex node) {
	auto dim = likelihoodDataDimension (sequenceNodes, modelNodes);
	auto childBranches = tree.childBranches (node);
	if (childBranches.empty ()) {
		// Leaf : lik data from sequences
		return DF::makeNode<DF::ConditionalLikelihoodFromSequence> (
		    {sequenceNodes.getSequenceNode (node)}, dim);
	} else {
		// Combine forward likelihoods
		return DF::makeNode<DF::ConditionalLikelihoodFromChildrens> (
		    mapToVector (childBranches,
		                 [&](TopologyBranchIndex branch) -> DF::NodeRef {
			                 return makeForwardLikelihoodNode (tree, sequenceNodes, brlenNodes,
			                                                   modelNodes, branch);
		                 }),
		    dim);
	}
}

DF::ValueRef<MatrixDouble> makeForwardLikelihoodNode (const TreeTopologyInterface & tree,
                                                      const SequenceNodeAccess & sequenceNodes,
                                                      const BranchLengthNodeAccess & brlenNodes,
                                                      const ModelNodeAccess & modelNodes,
                                                      TopologyBranchIndex branch) {
	auto dim = likelihoodDataDimension (sequenceNodes, modelNodes);
	auto transitionMatrix = DF::makeNode<DF::TransitionMatrixFromModel> (
	    {modelNodes.getModelNode (branch), brlenNodes.getBranchLengthNode (branch)},
	    TransitionMatrixDimension (modelNodes.getNbStates ()));
	return DF::makeNode<DF::ForwardLikelihoodFromChild> (
	    {std::move (transitionMatrix),
	     makeConditionalLikelihoodNode (tree, sequenceNodes, brlenNodes, modelNodes,
	                                    tree.childNode (branch))},
	    dim);
}
} // namespace bpp
