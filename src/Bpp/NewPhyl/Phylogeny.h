//
// File: Phylogeny.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-12
// Last modified: 2017-05-12
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

#pragma once
#ifndef BPP_NEWPHYL_PHYLOGENY_H
#define BPP_NEWPHYL_PHYLOGENY_H

#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/FrozenPtr.h>
#include <Bpp/NewPhyl/Likelihood.h>
#include <Bpp/NewPhyl/Model.h>
#include <Bpp/NewPhyl/NodeSpecification.h>
#include <Bpp/NewPhyl/TopologyMap.h>

namespace bpp {

class Sequence;
class SubstitutionModel;

namespace Phyl {

	struct Process {
		const FrozenPtr<Topology::Tree> tree;
		const FrozenPtr<Topology::BranchValueMap<DF::Parameter<double>>> branchLengths;
		const FrozenPtr<Topology::BranchValueMap<DF::Value<const SubstitutionModel *>>> modelByBranch;
		const std::size_t nbStates;
	};

	struct LikelihoodParameters {
		const Process process;
		const FrozenPtr<Topology::NodeValueMap<DF::Parameter<const Sequence *>>> leafData;
		const std::size_t nbSites;
	};

	// SPECS

	struct ConditionalLikelihoodSpec {
		const LikelihoodParameters likParams;
		const Topology::Node node;

		bool computed_from_data () const { return node.nbChildBranches () == 0; }
		DF::NodeSpecificationVec computeDependencies () const;
		DF::Node buildNode (DF::NodeVec deps) const {
			if (computed_from_data ())
				return DF::Node::create<ComputeConditionalLikelihoodFromDataNode> (
				    std::move (deps), likParams.nbSites, LikelihoodVector (likParams.process.nbStates));
			else
				return DF::Node::create<ComputeConditionalLikelihoodFromChildrensNode> (
				    std::move (deps), likParams.nbSites, LikelihoodVector (likParams.process.nbStates));
		}
		std::type_index nodeType () const {
			return computed_from_data () ? typeid (ComputeConditionalLikelihoodFromDataNode)
			                             : typeid (ComputeConditionalLikelihoodFromChildrensNode);
		}
		static std::string description () { return prettyTypeName<ConditionalLikelihoodSpec> (); }
	};

	struct ForwardLikelihoodSpec {
		const LikelihoodParameters likParams;
		const Topology::Branch branch;

		DF::NodeSpecificationVec computeDependencies () const {
			return DF::makeNodeSpecVec (
			    ConditionalLikelihoodSpec{likParams, branch.childNode ()},
			    ModelTransitionMatrixSpec (likParams.process.modelByBranch->access (branch).value (),
			                               likParams.process.branchLengths->access (branch).value (),
			                               likParams.process.nbStates));
		}
		DF::Node buildNode (DF::NodeVec deps) const {
			return DF::Node::create<ComputeForwardLikelihoodNode> (
			    std::move (deps), likParams.nbSites, LikelihoodVector (likParams.process.nbStates));
		}
		static std::type_index nodeType () { return typeid (ComputeForwardLikelihoodNode); }
		static std::string description () { return prettyTypeName<ForwardLikelihoodSpec> (); }
	};

	inline DF::NodeSpecificationVec ConditionalLikelihoodSpec::computeDependencies () const {
		if (computed_from_data ()) {
			return DF::makeNodeSpecVec (
			    DF::NodeSpecReturnParameter{likParams.leafData->access (node).value ()});
		} else {
			DF::NodeSpecificationVec depSpecs;
			node.foreachChildBranch ([this, &depSpecs](Topology::Branch && branch) {
				depSpecs.emplace_back (ForwardLikelihoodSpec{likParams, branch});
			});
			return depSpecs;
		}
	}

	struct LogLikelihoodSpec : DF::NodeSpecAlwaysGenerate<ComputeLogLikelihoodNode> {
		const LikelihoodParameters likParams;

		LogLikelihoodSpec (const LikelihoodParameters & params) : likParams (params) {}

		DF::NodeSpecificationVec computeDependencies () const {
			return DF::makeNodeSpecVec (
			    ConditionalLikelihoodSpec{likParams, likParams.process.tree->rootNode ()},
			    ModelEquilibriumFrequenciesSpec (
			        likParams.process.modelByBranch
			            ->access (likParams.process.tree->rootNode ().fatherBranch ())
			            .value (),
			        likParams.process.nbStates));
		}
	};
}
}

#endif // BPP_NEWPHYL_PHYLOGENY_H
