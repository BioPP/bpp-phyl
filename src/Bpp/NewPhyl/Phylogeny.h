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

/* Virtual views/access classes.
 *
 * These classes are interface which describe capabilities used by DF graph construction methods.
 */
class TreeTopologyView {
	/* Allows to move around a tree topology through node/branch indexes.
	 *
	 * These indexes abstract away what the underlying data structure uses.
	 * They are defined as a std::intptr_t.
	 * This means they can store any integer or pointer, which should be general enough.
	 * Indexes are represented by an *Index struct, which prevents implicit conversions (error prone).
	 */

public:
	struct NodeIndex {
		std::intptr_t value;
		explicit NodeIndex (std::intptr_t v) : value (v) {}
	};
	struct BranchIndex {
		std::intptr_t value;
		explicit BranchIndex (std::intptr_t v) : value (v) {}
	};

	virtual ~TreeTopologyView () = default;
	virtual NodeIndex rootNode () const = 0;
	virtual NodeIndex fatherNode (BranchIndex id) const = 0;
	virtual NodeIndex childNode (BranchIndex id) const = 0;
	virtual BranchIndex fatherBranch (NodeIndex id) const = 0;
	virtual Vector<BranchIndex> childBranches (NodeIndex id) const = 0;
};
bool operator== (TreeTopologyView::NodeIndex lhs, TreeTopologyView::NodeIndex rhs);
bool operator< (TreeTopologyView::NodeIndex lhs, TreeTopologyView::NodeIndex rhs);
bool operator== (TreeTopologyView::BranchIndex lhs, TreeTopologyView::BranchIndex rhs);
bool operator< (TreeTopologyView::BranchIndex lhs, TreeTopologyView::BranchIndex rhs);

class BranchLengthValueAccess {
	// Can access Branch length fixed value by branch id
public:
	virtual ~BranchLengthValueAccess () = default;
	virtual double getBranchLengthValue (TreeTopologyView::BranchIndex id) const = 0;
};

class BranchLengthNodeAccess {
	// Can access Branch length as a Value<double> DF node by branch id
public:
	virtual ~BranchLengthNodeAccess () = default;
	virtual DF::ValueRef<double> getBranchLengthNode (TreeTopologyView::BranchIndex id) const = 0;
};

class ModelNodeAccess {
	// Can access a Model DF node by branch id, has a defined number of states.
public:
	virtual ~ModelNodeAccess () = default;
	virtual DF::ValueRef<const SubstitutionModel *>
	getModelNode (TreeTopologyView::BranchIndex id) const = 0;
	virtual SizeType getNbStates () const = 0; // tree constant
};

class SequenceNameValueAccess {
	// Can access Sequence names at leaves.
public:
	virtual ~SequenceNameValueAccess () = default;
	virtual std::string getSequenceName (TreeTopologyView::NodeIndex id) const = 0;
};

class SequenceNodeAccess {
	// Can access Sequence DF node at leaves, has a defined number of sites.
public:
	virtual ~SequenceNodeAccess () = default;
	virtual DF::ValueRef<const Sequence *> getSequenceNode (TreeTopologyView::NodeIndex id) const = 0;
	virtual SizeType getNbSites () const = 0;
};

/* Common useful access classes.
 */
class SameModelForAllBranches : public ModelNodeAccess {
	// Use a single model for all branches of a tree.
public:
	SameModelForAllBranches (DF::ValueRef<const SubstitutionModel *> model);
	DF::ValueRef<const SubstitutionModel *>
	    getModelNode (TreeTopologyView::BranchIndex) const override;
	SizeType getNbStates () const override;

private:
	DF::ValueRef<const SubstitutionModel *> model_;
};

class BranchLengthParametersInitializedFromValues : public BranchLengthNodeAccess {
	/* Associate a DF::Parameter<double> to each branch.
	 * Parameters are initialised by values from a BranchLengthValueAccess object.
	 * Parameters are created lazily (when accessed).
	 */
public:
	BranchLengthParametersInitializedFromValues (const BranchLengthValueAccess & values);
	DF::ValueRef<double> getBranchLengthNode (TreeTopologyView::BranchIndex id) const override;
	DF::ParameterRef<double> getBranchLengthParameter (TreeTopologyView::BranchIndex id) const;

private:
	const BranchLengthValueAccess & values_;
	mutable std::map<TreeTopologyView::BranchIndex, DF::ParameterRef<double>> parameterNodes_;
};

class SequenceNodesInilialisedFromNames : public SequenceNodeAccess {
	/* Associate a DF::Constant<const Sequence*> to each leaf.
	 * Sequence are selected from a SiteContainer by names from a SequenceNameValueAccess.
	 * DF::Constant nodes are created lazily (when accessed).
	 */
public:
	SequenceNodesInilialisedFromNames (const SequenceNameValueAccess & names,
	                                   const SiteContainer & sequences);
	DF::ValueRef<const Sequence *> getSequenceNode (TreeTopologyView::NodeIndex id) const override;
	SizeType getNbSites () const override;

private:
	const SequenceNameValueAccess & names_;
	const SiteContainer & sequences_;
	mutable std::map<TreeTopologyView::NodeIndex, DF::ValueRef<const Sequence *>> sequenceNodes_;
};

/* Phylogeny DF graph construction functions.
 */
DF::ValueRef<MatrixDouble> makeConditionalLikelihoodNode (const TreeTopologyView & tree,
                                                          const SequenceNodeAccess & sequenceNodes,
                                                          const BranchLengthNodeAccess & brlenNodes,
                                                          const ModelNodeAccess & modelNodes,
                                                          TreeTopologyView::NodeIndex node);
DF::ValueRef<MatrixDouble> makeForwardLikelihoodNode (const TreeTopologyView & tree,
                                                      const SequenceNodeAccess & sequenceNodes,
                                                      const BranchLengthNodeAccess & brlenNodes,
                                                      const ModelNodeAccess & modelNodes,
                                                      TreeTopologyView::BranchIndex branch);
DF::ValueRef<double> makeLogLikelihoodNode (const TreeTopologyView & tree,
                                            const SequenceNodeAccess & sequenceNodes,
                                            const BranchLengthNodeAccess & brlenNodes,
                                            const ModelNodeAccess & modelNodes);
} // namespace bpp

#endif // BPP_NEWPHYL_PHYLOGENY_H
