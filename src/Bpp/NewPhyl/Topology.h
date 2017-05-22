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
#include <memory>
#include <stdexcept> // std::runtime_error
#include <string>
#include <vector>

namespace bpp {
namespace Topology {
	using IndexType = std::size_t; // For Node and Branch
	constexpr IndexType invalid{std::numeric_limits<IndexType>::max ()};

	class Node;
	class Branch;

	/* Tree class, stores structure.
	 * TODO improve, make this a virtual class with impls for subtrees, etc
   * TODO split branchs from nodes (for dag like structures), ids too
	 * TODO add "observers" (other name needed !), a value one, and an indexing one (names).
	 */
	class Tree {
	public:
		struct Node {
			IndexType fatherId_{invalid};
			std::vector<IndexType> childrenIds_{};
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
		IndexType createEdge (IndexType fatherId, IndexType childId) {
			auto & father = nodes_[fatherId];
			auto & child = nodes_[childId];
			if (child.fatherId_ != invalid)
				throw std::runtime_error ("child already has a father");
			child.fatherId_ = fatherId;
			father.childrenIds_.emplace_back (childId);
			return childId; // Edge id is id of child node
		}
		IndexType createNode (std::vector<IndexType> childrens) {
			auto id = createNode ();
			for (auto i : childrens)
				createEdge (id, i);
			return id;
		}

		static std::unique_ptr<Tree> create () { return std::unique_ptr<Tree>{new Tree}; }
		static std::shared_ptr<const Tree> finalize (std::unique_ptr<Tree> && tree) {
			return {std::move (tree)};
		}

	protected:
		Tree () = default;

	private:
		IndexType rootId_{invalid};
		std::vector<Node> nodes_{};
	};

	/* NodeRef and BranchRef.
	 * Kind of iterators on the tree, allow to inspect structure.
	 * TODO add && navigation ops (reuse tree shared_ptr mostly)
	 */
	class Node {
	public:
		Node (std::shared_ptr<const Tree> tree, IndexType nodeId) noexcept
		    : tree_ (tree), nodeId_ (nodeId) {}

		IndexType nodeId () const noexcept { return nodeId_; }
		std::size_t nbChildBranches () const noexcept {
			return tree_->node (nodeId_).childrenIds_.size ();
		}

		// Navigate
		Branch fatherBranch () const;
		Branch childBranch (std::size_t index) const;
		template <typename Callable> void foreachChildBranch (Callable callable) const;

	private:
		std::shared_ptr<const Tree> tree_;
		IndexType nodeId_;
	};

	class Branch {
	public:
		Branch (std::shared_ptr<const Tree> tree, IndexType childNodeId) noexcept
		    : tree_ (tree), childNodeId_ (childNodeId) {}

		IndexType fatherNodeId () const noexcept { return tree_->node (childNodeId_).fatherId_; }
		IndexType childNodeId () const noexcept { return childNodeId_; }

		// Navigate
		Node fatherNode () const;
		Node childNode () const;

	private:
		std::shared_ptr<const Tree> tree_;
		IndexType childNodeId_;
	};

	// Node navigation
	inline Branch Node::fatherBranch () const { return Branch (tree_, nodeId_); }
	inline Branch Node::childBranch (std::size_t index) const {
		return Branch (tree_, tree_->node (nodeId_).childrenIds_.at (index));
	}
	template <typename Callable> void Node::foreachChildBranch (Callable callable) const {
		for (auto childId : tree_->node (nodeId_).childrenIds_)
			callable (Branch (tree_, childId));
	}
	// Branch navigation
	inline Node Branch::fatherNode () const {
		auto id = fatherNodeId ();
		if (id == invalid)
			throw std::runtime_error ("branch has no father node");
		return Node (tree_, id);
	}
	inline Node Branch::childNode () const { return Node (tree_, childNodeId_); }
}
}

#endif // BPP_NEWPHYL_TOPOLOGY_H
