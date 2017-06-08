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

#include <Bpp/NewPhyl/FrozenPtr.h>
#include <Bpp/NewPhyl/Range.h>
#include <limits>
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
	 */
	class Tree : public std::enable_shared_from_this<Tree> {
	public:
		struct NodeData {
			IndexType fatherId_{invalid};
			std::vector<IndexType> childrenIds_{};
			NodeData () = default;
		};

		Tree () = default;

		// Base info
		IndexType nbNodes () const { return nodes_.size (); }
		IndexType nbBranches () const { return nbNodes (); }
		IndexType rootNodeId () const { return rootNodeId_; }

		// Branch info
		IndexType branchFatherNode (IndexType branchId) const {
			return nodes_.at (branchChildNode (branchId)).fatherId_;
		}
		IndexType branchChildNode (IndexType branchId) const { return branchId; }

		// Node info
		IndexType nodeFatherBranch (IndexType nodeId) const { return nodeId; }
		std::size_t nodeNbChildrenBranches (IndexType nodeId) const {
			return nodes_.at (nodeId).childrenIds_.size ();
		}
		IndexType nodeChildBranch (IndexType nodeId, std::size_t childBranchIndexInNode) const {
			return nodes_.at (nodeId).childrenIds_.at (childBranchIndexInNode);
		}

		// Setup
		void setRootNodeId (IndexType nodeId) {
			if (rootNodeId_ != invalid)
				throw std::runtime_error ("root node has already been set");
			rootNodeId_ = nodeId;
		}
		IndexType createNode () {
			auto id = nbNodes ();
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

		// "Iterators" only use if is a shared_ptr...
		Node node (IndexType id) const;
		Node rootNode () const;
		Branch branch (IndexType id) const;

	private:
		IndexType rootNodeId_{invalid};
		std::vector<NodeData> nodes_{};
	};

	/* NodeRef and BranchRef.
	 * Kind of iterators on the tree, allow to inspect structure.
	 */
	class Node {
	public:
		Node (FrozenSharedPtr<Tree> tree, IndexType nodeId) noexcept
		    : tree_ (std::move (tree)), nodeId_ (nodeId) {}

		IndexType nodeId () const noexcept { return nodeId_; }
		IndexType fatherBranchId () const { return tree_->nodeFatherBranch (nodeId ()); }
		std::size_t nbChildBranches () const { return tree_->nodeNbChildrenBranches (nodeId ()); }
		IndexType childBranchId (std::size_t branchIndex) const {
			return tree_->nodeChildBranch (nodeId (), branchIndex);
		}
		const FrozenSharedPtr<Tree> & tree () const noexcept { return tree_; }

		// Navigate
		Branch fatherBranch () const &;
		Branch fatherBranch () &&;
		Branch childBranch (std::size_t index) const &;
		Branch childBranch (std::size_t index) &&;
		template <typename Callable> void foreachChildBranch (Callable callable) const;

	private:
		Branch buildBranch (IndexType branchId) const &;
		Branch buildBranch (IndexType branchId) &&;

		FrozenSharedPtr<Tree> tree_;
		IndexType nodeId_;
	};

	class Branch {
	public:
		Branch (FrozenSharedPtr<Tree> tree, IndexType branchId) noexcept
		    : tree_ (std::move (tree)), branchId_ (branchId) {}

		IndexType branchId () const noexcept { return branchId_; }
		IndexType fatherNodeId () const { return tree_->branchFatherNode (branchId ()); }
		IndexType childNodeId () const { return tree_->branchChildNode (branchId ()); }
		const FrozenSharedPtr<Tree> & tree () const noexcept { return tree_; }

		// Navigate
		Node fatherNode () const &;
		Node fatherNode () &&;
		Node childNode () const &;
		Node childNode () &&;

	private:
		Node buildNode (IndexType nodeId) const &;
		Node buildNode (IndexType nodeId) &&;

		FrozenSharedPtr<Tree> tree_;
		IndexType branchId_;
	};

	// Tree "iterators"
	inline Node Tree::node (IndexType id) const {
		return Node{FrozenSharedPtr<Tree>::shared_from_this (*this), id};
	}
	inline Node Tree::rootNode () const {
		return Node{FrozenSharedPtr<Tree>::shared_from_this (*this), rootNodeId ()};
	}
	inline Branch Tree::branch (IndexType id) const {
		return Branch{FrozenSharedPtr<Tree>::shared_from_this (*this), id};
	}

	// Node navigation
	inline Branch Node::fatherBranch () const & { return buildBranch (fatherBranchId ()); }
	inline Branch Node::fatherBranch () && {
		return std::move (*this).buildBranch (fatherBranchId ());
	}
	inline Branch Node::childBranch (std::size_t index) const & {
		return buildBranch (childBranchId (index));
	}
	inline Branch Node::childBranch (std::size_t index) && {
		return std::move (*this).buildBranch (childBranchId (index));
	}
	template <typename Callable> void Node::foreachChildBranch (Callable callable) const {
		for (auto branchIndex : bpp::range (nbChildBranches ()))
			callable (childBranch (branchIndex));
	}
	inline Branch Node::buildBranch (IndexType branchId) const & { return Branch (tree_, branchId); }
	inline Branch Node::buildBranch (IndexType branchId) && {
		return Branch (std::move (tree_), branchId);
	}

	// Branch navigation
	inline Node Branch::fatherNode () const & { return buildNode (fatherNodeId ()); }
	inline Node Branch::fatherNode () && { return std::move (*this).buildNode (fatherNodeId ()); }
	inline Node Branch::childNode () const & { return buildNode (childNodeId ()); }
	inline Node Branch::childNode () && { return std::move (*this).buildNode (childNodeId ()); }
	inline Node Branch::buildNode (IndexType nodeId) const & { return Node (tree_, nodeId); }
	inline Node Branch::buildNode (IndexType nodeId) && { return Node (std::move (tree_), nodeId); }
}
}

#endif // BPP_NEWPHYL_TOPOLOGY_H
