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

#ifndef BPP_NEWPHYL_PHYLOGENY_H
#define BPP_NEWPHYL_PHYLOGENY_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h> // Matrix types
#include <Bpp/NewPhyl/DataFlowTemplates.h>
#include <Bpp/NewPhyl/FrozenPtr.h>
#include <Bpp/NewPhyl/Signed.h>
#include <Bpp/NewPhyl/Topology.h>
#include <Bpp/NewPhyl/TopologyMap.h>
#include <string>

namespace bpp {

// Forward declarations
class Sequence;
class SubstitutionModel;

namespace Phyl {

	// Phylogeny encoding structs
	// TODO remove frozen ptr, just bundle values
	struct Process {
		const FrozenPtr<Topology::Tree> tree;
		const FrozenPtr<Topology::BranchValueMap<DF::ParameterRef<double>>> branchLengths;
		const FrozenPtr<Topology::BranchValueMap<DF::ValueRef<const SubstitutionModel *>>>
		    modelByBranch;
		const SizeType nbStates;
	};

	struct SequenceMap {
		const FrozenPtr<Topology::NodeValueMap<DF::ParameterRef<const Sequence *>>> sequences;
		const SizeType nbSites;
	};

	struct LikelihoodParameters {
		const Process process;
		const SequenceMap leafData;
	};

	// Build functions

	DF::ValueRef<double> makeLogLikelihoodNode (const LikelihoodParameters & params);
	DF::ValueRef<DF::MatrixDouble> makeConditionalLikelihoodNode (const LikelihoodParameters & params,
	                                                              Topology::Node node);
	DF::ValueRef<DF::MatrixDouble> makeForwardLikelihoodNode (const LikelihoodParameters & params,
	                                                          Topology::Branch branch);
} // namespace Phyl
} // namespace bpp

#endif // BPP_NEWPHYL_PHYLOGENY_H
