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
#include <Bpp/NewPhyl/LinearAlgebraFwd.h>
#include <Bpp/NewPhyl/Signed.h>
#include <Bpp/NewPhyl/Vector.h>
#include <cstdint> // intptr_t
#include <map>
#include <string>

namespace bpp {

// Forward declarations
class Sequence;
class SubstitutionModel;
class SiteContainer;

/** Data structure independent index for nodes in a topology.
 * Implemented as a intptr_t to be able to hold both integers and pointer-like values.
 * Requires an explicit conversion to avoid unintended implicit conversions.
 */
struct TopologyNodeIndex {
	std::intptr_t value;
	explicit TopologyNodeIndex (std::intptr_t v) : value (v) {}
};

/// Data structure independent index for branches. Equivalent to TopologyNodeIndex.
struct TopologyBranchIndex {
	std::intptr_t value;
	explicit TopologyBranchIndex (std::intptr_t v) : value (v) {}
};

// Comparison operators
bool operator== (TopologyNodeIndex lhs, TopologyNodeIndex rhs);
bool operator< (TopologyNodeIndex lhs, TopologyNodeIndex rhs);
bool operator== (TopologyBranchIndex lhs, TopologyBranchIndex rhs);
bool operator< (TopologyBranchIndex lhs, TopologyBranchIndex rhs);

/** Virtual interface to a structure representing a read-only tree topology.
 * Describes how to navigate from root to leaves, and back up.
 * Navigation is encoded with TopologyNodeIndex and TopologyBranchIndex types.
 * Indexes have no "invalid" state: the tree topology must be valid and non empty.
 */
class TreeTopologyInterface {
public:
	virtual ~TreeTopologyInterface () = default;
	virtual TopologyNodeIndex rootNode () const = 0;
	virtual TopologyNodeIndex fatherNode (TopologyBranchIndex id) const = 0;
	virtual TopologyNodeIndex childNode (TopologyBranchIndex id) const = 0;
	virtual TopologyBranchIndex fatherBranch (TopologyNodeIndex id) const = 0;
	virtual Vector<TopologyBranchIndex> childBranches (TopologyNodeIndex id) const = 0;
};

/// Interface: Can access Branch length fixed value by branch index.
class BranchLengthValueAccess {
public:
	virtual ~BranchLengthValueAccess () = default;
	virtual double getBranchLengthValue (TopologyBranchIndex id) const = 0;
};

/// Interface: Can access Branch length as a Value<double> DF node by branch index.
class BranchLengthNodeAccess {
public:
	virtual ~BranchLengthNodeAccess () = default;
	virtual DF::ValueRef<double> getBranchLengthNode (TopologyBranchIndex id) const = 0;
};

/// Interface: Can access a Model DF node by branch id, has a defined number of states.
class ModelNodeAccess {
public:
	virtual ~ModelNodeAccess () = default;
	virtual DF::ValueRef<const SubstitutionModel *> getModelNode (TopologyBranchIndex id) const = 0;
	virtual SizeType getNbStates () const = 0; // tree constant
};

/// Interface: Can access Sequence names at leaves.
class SequenceNameValueAccess {
public:
	virtual ~SequenceNameValueAccess () = default;
	virtual std::string getSequenceName (TopologyNodeIndex id) const = 0;
};

/// Interface: Can access Sequence DF node at leaves, has a defined number of sites.
class SequenceNodeAccess {
public:
	virtual ~SequenceNodeAccess () = default;
	virtual DF::ValueRef<const Sequence *> getSequenceNode (TopologyNodeIndex id) const = 0;
	virtual SizeType getNbSites () const = 0;
};

/** ModelNodeAccess implementation using one model node for all branches.
 */
class SameModelForAllBranches : public ModelNodeAccess {
public:
	SameModelForAllBranches (DF::ValueRef<const SubstitutionModel *> model);
	DF::ValueRef<const SubstitutionModel *> getModelNode (TopologyBranchIndex) const override;
	SizeType getNbStates () const override;

private:
	DF::ValueRef<const SubstitutionModel *> model_;
};

/** BranchLengthNodeAccess impl: creates one DF::Parameter for each branch from values.
 * Associate a DF::Parameter<double> to each branch.
 * Parameters are initialised by values from a BranchLengthValueAccess object.
 * Parameters are created lazily (when accessed).
 */
class BranchLengthParametersInitializedFromValues : public BranchLengthNodeAccess {
public:
	BranchLengthParametersInitializedFromValues (const BranchLengthValueAccess & values);
	DF::ValueRef<double> getBranchLengthNode (TopologyBranchIndex id) const override;
	DF::ParameterRef<double> getBranchLengthParameter (TopologyBranchIndex id) const;

private:
	const BranchLengthValueAccess & values_;
	mutable std::map<TopologyBranchIndex, DF::ParameterRef<double>> parameterNodes_;
};

/** SequenceNodeAccess impl: creates one sequence node for each leave from names.
 * Associate a DF::Constant<const Sequence*> to each leaf.
 * Sequence are selected from a SiteContainer by names from a SequenceNameValueAccess.
 * DF::Constant nodes are created lazily (when accessed).
 */
class SequenceNodesInilialisedFromNames : public SequenceNodeAccess {
public:
	SequenceNodesInilialisedFromNames (const SequenceNameValueAccess & names,
	                                   const SiteContainer & sequences);
	DF::ValueRef<const Sequence *> getSequenceNode (TopologyNodeIndex id) const override;
	SizeType getNbSites () const override;

private:
	const SequenceNameValueAccess & names_;
	const SiteContainer & sequences_;
	mutable std::map<TopologyNodeIndex, DF::ValueRef<const Sequence *>> sequenceNodes_;
};

/* Phylogeny DF graph construction functions.
 * TODO doc
 */
DF::ValueRef<MatrixDouble> makeConditionalLikelihoodNode (const TreeTopologyInterface & tree,
                                                          const SequenceNodeAccess & sequenceNodes,
                                                          const BranchLengthNodeAccess & brlenNodes,
                                                          const ModelNodeAccess & modelNodes,
                                                          TopologyNodeIndex node);
DF::ValueRef<MatrixDouble> makeForwardLikelihoodNode (const TreeTopologyInterface & tree,
                                                      const SequenceNodeAccess & sequenceNodes,
                                                      const BranchLengthNodeAccess & brlenNodes,
                                                      const ModelNodeAccess & modelNodes,
                                                      TopologyBranchIndex branch);
DF::ValueRef<double> makeLogLikelihoodNode (const TreeTopologyInterface & tree,
                                            const SequenceNodeAccess & sequenceNodes,
                                            const BranchLengthNodeAccess & brlenNodes,
                                            const ModelNodeAccess & modelNodes);
} // namespace bpp

#endif // BPP_NEWPHYL_PHYLOGENY_H
