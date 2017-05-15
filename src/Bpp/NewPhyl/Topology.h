//
// File: Topology.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-26
// Last modified: 2017-04-26
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
#ifndef BPP_NEWPHYL_TOPOLOGY_H
#define BPP_NEWPHYL_TOPOLOGY_H

#include <limits>
#include <stdexcept> // std::runtime_error
#include <string>
#include <vector>

namespace bpp {
namespace Topology {
	using IndexType = std::size_t;
	constexpr IndexType invalid{std::numeric_limits<IndexType>::max ()};

	class NodeRef;
	class BranchRef;

	/* Tree class, stores structure.
	 * TODO improve...
	 */
	class Tree {
	public:
		struct Node {
			IndexType fatherId_{invalid};
			std::vector<IndexType> childrenIds_{};
			std::string nodeName_{};

			Node () = default;
		};

		std::size_t nbNodes () const noexcept { return nodes_.size (); }
		const Node & node (IndexType id) const noexcept { return nodes_[id]; }
		const IndexType & rootId () const noexcept { return rootId_; }
		IndexType & rootId () noexcept { return rootId_; }

		IndexType createNode () {
			auto id = IndexType (nbNodes ());
			nodes_.emplace_back ();
			return id;
		}
		void createEdge (IndexType fatherId, IndexType childId) {
			auto & father = nodes_[fatherId];
			auto & child = nodes_[childId];
			if (child.fatherId_ != invalid)
				throw std::runtime_error ("child already has a father");
			child.fatherId_ = fatherId;
			father.childrenIds_.emplace_back (childId);
		}
		void setNodeName (IndexType id, std::string name) noexcept {
			nodes_[id].nodeName_ = std::move (name);
		}
		IndexType createNode (std::vector<IndexType> childrens, std::string name = std::string ()) {
			auto id = createNode ();
			setNodeName (id, std::move (name));
			for (auto i : childrens)
				createEdge (id, i);
			return id;
		}

		NodeRef nodeRef (IndexType nodeId) const noexcept;

	private:
		IndexType rootId_{invalid};
		std::vector<Node> nodes_{};
	};

	/* NodeRef and BranchRef.
	 * Kind of iterators on the tree, allow to inspect structure.
	 */
	class NodeRef {
	public:
		NodeRef (const Tree & tree, IndexType nodeId) noexcept : tree_ (tree), nodeId_ (nodeId) {}

		IndexType nodeId () const noexcept { return nodeId_; }
		IndexType nbChildBranches () const noexcept {
			return tree_.node (nodeId_).childrenIds_.size ();
		}
		const std::string & name () const { return tree_.node (nodeId_).nodeName_; }

		// Navigate
		BranchRef fatherBranch () const;
		BranchRef childBranch (IndexType id) const;

	private:
		const Tree & tree_;
		IndexType nodeId_;
	};

	inline NodeRef Tree::nodeRef (IndexType nodeId) const noexcept { return NodeRef (*this, nodeId); }

	class BranchRef {
	public:
		BranchRef (const Tree & tree, IndexType childNodeId) noexcept
		    : tree_ (tree), childNodeId_ (childNodeId) {}

		IndexType fatherNodeId () const noexcept { return tree_.node (childNodeId_).fatherId_; }
		IndexType childNodeId () const noexcept { return childNodeId_; }

		// Navigate
		NodeRef fatherNode () const {
			auto id = fatherNodeId ();
			if (id == invalid)
				throw std::runtime_error ("branch has no father node");
			return NodeRef (tree_, id);
		}
		NodeRef childNode () const { return NodeRef (tree_, childNodeId_); }

	private:
		const Tree & tree_;
		IndexType childNodeId_;
	};

	inline BranchRef NodeRef::fatherBranch () const { return BranchRef (tree_, nodeId_); }
	inline BranchRef NodeRef::childBranch (IndexType id) const {
		return BranchRef (tree_, tree_.node (nodeId_).childrenIds_[id]);
	}
}
}

#endif // BPP_NEWPHYL_TOPOLOGY_H
