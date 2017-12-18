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
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <utility>

namespace bpp {

// Index comparisons
bool operator== (TopologyNodeIndex lhs, TopologyNodeIndex rhs) {
	return lhs.value == rhs.value;
}
bool operator< (TopologyNodeIndex lhs, TopologyNodeIndex rhs) {
	return lhs.value < rhs.value;
}
bool operator== (TopologyBranchIndex lhs, TopologyBranchIndex rhs) {
	return lhs.value == rhs.value;
}
bool operator< (TopologyBranchIndex lhs, TopologyBranchIndex rhs) {
	return lhs.value < rhs.value;
}

// Access class implementations
SameModelForAllBranches::SameModelForAllBranches (DF::ValueRef<const SubstitutionModel *> model)
    : model_ (std::move (model)) {}
DF::ValueRef<const SubstitutionModel *>
SameModelForAllBranches::getModelNode (TopologyBranchIndex) const {
	return model_;
}
SizeType SameModelForAllBranches::getNbStates () const {
	return static_cast<SizeType> (model_->getValue ()->getNumberOfStates ());
}

BranchLengthsInitializedFromValues::BranchLengthsInitializedFromValues (
    const BranchLengthValueAccess & values)
    : values_ (values) {}
DF::ValueRef<double>
BranchLengthsInitializedFromValues::getBranchLengthNode (TopologyBranchIndex id) const {
	return getBranchLengthMutableNode (id);
}
DF::MutableRef<double>
BranchLengthsInitializedFromValues::getBranchLengthMutableNode (TopologyBranchIndex id) const {
	auto it = mutableNodes_.find (id);
	if (it != mutableNodes_.end ()) {
		return it->second;
	} else {
		auto mutableNode = DF::makeNode<DF::Mutable<double>> (values_.getBranchLengthValue (id));
		mutableNodes_.emplace (id, mutableNode);
		return mutableNode;
	}
}

SequenceNodesInilialisedFromNames::SequenceNodesInilialisedFromNames (
    const SequenceNameValueAccess & names, const SiteContainer & sequences)
    : names_ (names), sequences_ (sequences) {}
DF::ValueRef<const Sequence *>
SequenceNodesInilialisedFromNames::getSequenceNode (TopologyNodeIndex id) const {
	auto it = sequenceNodes_.find (id);
	if (it != sequenceNodes_.end ()) {
		return it->second;
	} else {
		auto sequence = DF::makeNode<DF::Constant<const Sequence *>> (
		    &sequences_.getSequence (names_.getSequenceName (id)));
		sequenceNodes_.emplace (id, sequence);
		return sequence;
	}
}
SizeType SequenceNodesInilialisedFromNames::getNbSites () const {
	return static_cast<SizeType> (sequences_.getNumberOfSites ());
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

DF::ValueRef<double> makeLogLikelihoodNode (const TreeTopologyInterface & tree,
                                            const SequenceNodeAccess & sequenceNodes,
                                            const BranchLengthNodeAccess & brlenNodes,
                                            const ModelNodeAccess & modelNodes) {
	auto dim = likelihoodDataDimension (sequenceNodes, modelNodes);
	auto root = tree.rootNode ();
	auto rootBranchModel = modelNodes.getModelNode (tree.fatherBranch (root));
	auto rootEquilibriumFrequencies =
	    DF::makeNode<DF::EquilibriumFrequenciesFromModel> ({rootBranchModel}, dim.nbStates ());
	auto likelihood = DF::makeNode<DF::Likelihood> (
	    {makeConditionalLikelihoodNode (tree, sequenceNodes, brlenNodes, modelNodes, root),
	     std::move (rootEquilibriumFrequencies)},
	    VectorDimension (dim.nbSites ()));
	auto logLik = DF::makeNode<DF::TotalLogLikelihood> ({std::move (likelihood)});
	return DF::makeNode<DF::NegDouble> ({std::move (logLik)});
}
} // namespace bpp
