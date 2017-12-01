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
// Access classes

bool operator== (TreeTopologyView::NodeIndex lhs, TreeTopologyView::NodeIndex rhs) {
	return lhs.value == rhs.value;
}
bool operator< (TreeTopologyView::NodeIndex lhs, TreeTopologyView::NodeIndex rhs) {
	return lhs.value < rhs.value;
}
bool operator== (TreeTopologyView::BranchIndex lhs, TreeTopologyView::BranchIndex rhs) {
	return lhs.value == rhs.value;
}
bool operator< (TreeTopologyView::BranchIndex lhs, TreeTopologyView::BranchIndex rhs) {
	return lhs.value < rhs.value;
}

SameModelForAllBranches::SameModelForAllBranches (DF::ValueRef<const SubstitutionModel *> model)
    : model_ (std::move (model)) {}
DF::ValueRef<const SubstitutionModel *>
SameModelForAllBranches::getModelNode (TreeTopologyView::BranchIndex) const {
	return model_;
}
SizeType SameModelForAllBranches::getNbStates () const {
	return static_cast<SizeType> (model_->getValue ()->getNumberOfStates ());
}

BranchLengthParametersInitializedFromValues::BranchLengthParametersInitializedFromValues (
    const BranchLengthValueAccess & values)
    : values_ (values) {}
DF::ValueRef<double> BranchLengthParametersInitializedFromValues::getBranchLengthNode (
    TreeTopologyView::BranchIndex id) const {
	return getBranchLengthParameter (id);
}
DF::ParameterRef<double> BranchLengthParametersInitializedFromValues::getBranchLengthParameter (
    TreeTopologyView::BranchIndex id) const {
	auto it = parameterNodes_.find (id);
	if (it != parameterNodes_.end ()) {
		return it->second;
	} else {
		auto parameter = DF::makeNode<DF::Parameter<double>> (values_.getBranchLengthValue (id));
		parameterNodes_.emplace (id, parameter);
		return parameter;
	}
}

SequenceNodesInilialisedFromNames::SequenceNodesInilialisedFromNames (
    const SequenceNameValueAccess & names, const SiteContainer & sequences)
    : names_ (names), sequences_ (sequences) {}
DF::ValueRef<const Sequence *>
SequenceNodesInilialisedFromNames::getSequenceNode (TreeTopologyView::NodeIndex id) const {
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

DF::ValueRef<MatrixDouble> makeConditionalLikelihoodNode (const TreeTopologyView & tree,
                                                          const SequenceNodeAccess & sequenceNodes,
                                                          const BranchLengthNodeAccess & brlenNodes,
                                                          const ModelNodeAccess & modelNodes,
                                                          TreeTopologyView::NodeIndex node) {
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
		                 [&](TreeTopologyView::BranchIndex branch) -> DF::NodeRef {
			                 return makeForwardLikelihoodNode (tree, sequenceNodes, brlenNodes,
			                                                   modelNodes, branch);
		                 }),
		    dim);
	}
}

DF::ValueRef<MatrixDouble> makeForwardLikelihoodNode (const TreeTopologyView & tree,
                                                      const SequenceNodeAccess & sequenceNodes,
                                                      const BranchLengthNodeAccess & brlenNodes,
                                                      const ModelNodeAccess & modelNodes,
                                                      TreeTopologyView::BranchIndex branch) {
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

DF::ValueRef<double> makeLogLikelihoodNode (const TreeTopologyView & tree,
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

namespace Phyl {
	namespace {
		DF::ValueRef<MatrixDouble> makeConditionalLikelihoodNode (const LikelihoodParameters & params,
		                                                          Topology::Node node);
		DF::ValueRef<MatrixDouble> makeForwardLikelihoodNode (const LikelihoodParameters & params,
		                                                      Topology::Branch branch);

		LikelihoodDataDimension dimensions (const LikelihoodParameters & params) {
			return {params.leafData.nbSites, params.process.nbStates};
		}

		DF::ValueRef<MatrixDouble> makeConditionalLikelihoodNode (const LikelihoodParameters & params,
		                                                          Topology::Node node) {
			auto dim = dimensions (params);
			if (node.nbChildBranches () == 0) {
				auto & sequence = params.leafData.sequences->access (node).value ();
				return DF::makeNode<DF::ConditionalLikelihoodFromSequence> ({sequence}, dim);
			} else {
				DF::NodeRefVec deps;
				node.foreachChildBranch ([&deps, &params](Topology::Branch && branch) {
					deps.emplace_back (makeForwardLikelihoodNode (params, std::move (branch)));
				});
				return DF::makeNode<DF::ConditionalLikelihoodFromChildrens> (std::move (deps), dim);
			}
		}

		DF::ValueRef<MatrixDouble> makeForwardLikelihoodNode (const LikelihoodParameters & params,
		                                                      Topology::Branch branch) {
			auto dim = dimensions (params);
			auto conditionalLikelihood = makeConditionalLikelihoodNode (params, branch.childNode ());

			auto & branchModel = params.process.modelByBranch->access (branch).value ();
			auto & brlen = params.process.branchLengths->access (branch).value ();
			auto modelTransitionMatrix = DF::makeNode<DF::TransitionMatrixFromModel> (
			    {branchModel, brlen}, TransitionMatrixDimension (dim.nbStates ()));

			return DF::makeNode<DF::ForwardLikelihoodFromChild> (
			    {std::move (modelTransitionMatrix), std::move (conditionalLikelihood)}, dim);
		}
	} // namespace

	DF::ValueRef<double> makeLogLikelihoodNode (const LikelihoodParameters & params) {
		auto dim = dimensions (params);
		auto & rootBranchModel =
		    params.process.modelByBranch->access (params.process.tree->rootNode ().fatherBranch ())
		        .value ();
		auto rootEquilibriumFrequencies =
		    DF::makeNode<DF::EquilibriumFrequenciesFromModel> ({rootBranchModel}, dim.nbStates ());
		auto likelihood = DF::makeNode<DF::Likelihood> (
		    {makeConditionalLikelihoodNode (params, params.process.tree->rootNode ()),
		     std::move (rootEquilibriumFrequencies)},
		    dim.nbSites ());
		auto logLik = DF::makeNode<DF::TotalLogLikelihood> ({std::move (likelihood)});
		return DF::makeNode<DF::NegDouble> ({std::move (logLik)});
	}

} // namespace Phyl
} // namespace bpp
