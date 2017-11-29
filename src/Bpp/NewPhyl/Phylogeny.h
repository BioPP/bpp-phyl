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
#include <Bpp/NewPhyl/DataFlowTemplates.h>
#include <Bpp/NewPhyl/Signed.h>
#include <Bpp/NewPhyl/Vector.h>

#include <Bpp/NewPhyl/FrozenPtr.h>
#include <Bpp/NewPhyl/Topology.h>
#include <Bpp/NewPhyl/TopologyMap.h>
#include <string>

namespace bpp {

// Forward declarations
class Sequence;
class SubstitutionModel;

/* Virtual views/access classes.
 *
 * These classes are interface which describe capabilities used by DF graph construction methods.
 */
class TreeTopologyView {
	// Allows to move around a tree topology through node/branch indexes.
public:
	virtual ~TreeTopologyView () = default;
	virtual IndexType rootNode () const = 0;
	virtual bool validBranchIndex (IndexType branchId) const = 0;
	virtual IndexType branchFatherNode (IndexType branchId) const = 0;
	virtual IndexType branchChildNode (IndexType branchId) const = 0;
	virtual bool validNodeIndex (IndexType nodeId) const = 0;
	virtual IndexType nodeFatherBranch (IndexType nodeId) const = 0;
	virtual Vector<IndexType> nodeChildBranches (IndexType nodeId) const = 0;
};
class BranchLengthNodeAccess {
	// Can access Branch length value node by branch id
public:
	virtual ~BranchLengthNodeAccess () = default;
	virtual const DF::ValueRef<double> & getBranchLengthNode (IndexType branchId) const = 0;
};
class ModelNodeAccess {
	// Can access a Model node by branch id, has a defined number of states.
public:
	virtual ~ModelNodeAccess () = default;
	virtual const DF::ValueRef<const SubstitutionModel *> &
	getModelNode (IndexType branchId) const = 0;
	virtual SizeType getNbStates () const = 0; // tree constant
};
class SequenceNodeAccess {
	// Can access Sequence nodes at leaves, has a defined number of sites.
public:
	virtual ~SequenceNodeAccess () = default;
	virtual const DF::ValueRef<const Sequence *> & getSequenceNode (IndexType nodeId) const = 0;
	virtual SizeType getNbSites () const = 0;
};

/* Common useful access classes.
 */

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
} // namespace Phyl
} // namespace bpp

#endif // BPP_NEWPHYL_PHYLOGENY_H
