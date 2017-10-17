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
#include <Bpp/NewPhyl/Signed.h>
#include <Bpp/NewPhyl/Vector.h>
#include <limits>
#include <stdexcept> // std::runtime_error
#include <string>

namespace bpp {
namespace Topology {
	// Index references a branch or a node in a tree
	using Index = IndexType;
	constexpr Index invalid{-1}; // TODO move to a struct to constrain ?

	class Node;
	class Branch;

	/* Tree class, stores structure.
	 * TODO improve, make this a virtual class with impls for subtrees, etc
	 */
	class Tree : public std::enable_shared_from_this<Tree> {
	public:
		struct NodeData {
			Index fatherId_{invalid};
			Vector<Index> childrenIds_{};
			NodeData () = default;
		};

		Tree () = default;

		// Base info
		SizeType nbNodes () const { return nodes_.size (); }
		SizeType nbBranches () const { return nbNodes (); }
		Index rootNodeId () const { return rootNodeId_; }

		// Branch info
		Index branchFatherNode (Index branchId) const {
			return nodes_.at (branchChildNode (branchId)).fatherId_;
		}
		Index branchChildNode (Index branchId) const { return branchId; }

		// Node info
		Index nodeFatherBranch (Index nodeId) const { return nodeId; }
		SizeType nodeNbChildrenBranches (Index nodeId) const {
			return nodes_.at (nodeId).childrenIds_.size ();
		}
		Index nodeChildBranch (Index nodeId, IndexType childBranchIndexInNode) const {
			return nodes_.at (nodeId).childrenIds_.at (childBranchIndexInNode);
		}

		// Setup
		void setRootNodeId (Index nodeId) {
			if (rootNodeId_ != invalid)
				throw std::runtime_error ("root node has already been set");
			rootNodeId_ = nodeId;
		}
		Index createNode () {
			auto id = nbNodes ();
			nodes_.emplace_back ();
			return id;
		}
		Index createEdge (Index fatherId, Index childId) {
			auto & father = nodes_[fatherId];
			auto & child = nodes_[childId];
			if (child.fatherId_ != invalid)
				throw std::runtime_error ("child already has a father");
			child.fatherId_ = fatherId;
			father.childrenIds_.emplace_back (childId);
			return childId; // Edge id is id of child node
		}
		Index createNode (Vector<Index> childrens) {
			auto id = createNode ();
			for (auto i : childrens)
				createEdge (id, i);
			return id;
		}

		// "Iterators" only use if is a shared_ptr...
		Node node (Index id) const;
		Node rootNode () const;
		Branch branch (Index id) const;

	private:
		Index rootNodeId_{invalid};
		Vector<NodeData> nodes_{};
	};

	/* NodeRef and BranchRef.
	 * Kind of iterators on the tree, allow to inspect structure.
	 */
	class Node {
	public:
		Node (FrozenPtr<Tree> tree, Index nodeId) noexcept
		    : tree_ (std::move (tree)), nodeId_ (nodeId) {}

		Index nodeId () const noexcept { return nodeId_; }
		Index fatherBranchId () const { return tree_->nodeFatherBranch (nodeId ()); }
		SizeType nbChildBranches () const { return tree_->nodeNbChildrenBranches (nodeId ()); }
		Index childBranchId (IndexType branchIndex) const {
			return tree_->nodeChildBranch (nodeId (), branchIndex);
		}
		const FrozenPtr<Tree> & tree () const noexcept { return tree_; }

		// Navigate
		Branch fatherBranch () const &;
		Branch fatherBranch () &&;
		Branch childBranch (IndexType index) const &;
		Branch childBranch (IndexType index) &&;
		template <typename Callable> void foreachChildBranch (Callable callable) const;

	private:
		Branch buildBranch (Index branchId) const &;
		Branch buildBranch (Index branchId) &&;

		FrozenPtr<Tree> tree_;
		Index nodeId_;
	};

	class Branch {
	public:
		Branch (FrozenPtr<Tree> tree, Index branchId) noexcept
		    : tree_ (std::move (tree)), branchId_ (branchId) {}

		Index branchId () const noexcept { return branchId_; }
		Index fatherNodeId () const { return tree_->branchFatherNode (branchId ()); }
		Index childNodeId () const { return tree_->branchChildNode (branchId ()); }
		const FrozenPtr<Tree> & tree () const noexcept { return tree_; }

		// Navigate
		Node fatherNode () const &;
		Node fatherNode () &&;
		Node childNode () const &;
		Node childNode () &&;

	private:
		Node buildNode (Index nodeId) const &;
		Node buildNode (Index nodeId) &&;

		FrozenPtr<Tree> tree_;
		Index branchId_;
	};

	// Tree "iterators"
	inline Node Tree::node (Index id) const {
		return Node{FrozenPtr<Tree>::shared_from_this (*this), id};
	}
	inline Node Tree::rootNode () const {
		return Node{FrozenPtr<Tree>::shared_from_this (*this), rootNodeId ()};
	}
	inline Branch Tree::branch (Index id) const {
		return Branch{FrozenPtr<Tree>::shared_from_this (*this), id};
	}

	// Node navigation
	inline Branch Node::fatherBranch () const & { return buildBranch (fatherBranchId ()); }
	inline Branch Node::fatherBranch () && {
		return std::move (*this).buildBranch (fatherBranchId ());
	}
	inline Branch Node::childBranch (IndexType index) const & {
		return buildBranch (childBranchId (index));
	}
	inline Branch Node::childBranch (IndexType index) && {
		return std::move (*this).buildBranch (childBranchId (index));
	}
	template <typename Callable> void Node::foreachChildBranch (Callable callable) const {
		for (auto branchIndex : bpp::range (nbChildBranches ()))
			callable (childBranch (branchIndex));
	}
	inline Branch Node::buildBranch (Index branchId) const & { return Branch (tree_, branchId); }
	inline Branch Node::buildBranch (Index branchId) && {
		return Branch (std::move (tree_), branchId);
	}

	// Branch navigation
	inline Node Branch::fatherNode () const & { return buildNode (fatherNodeId ()); }
	inline Node Branch::fatherNode () && { return std::move (*this).buildNode (fatherNodeId ()); }
	inline Node Branch::childNode () const & { return buildNode (childNodeId ()); }
	inline Node Branch::childNode () && { return std::move (*this).buildNode (childNodeId ()); }
	inline Node Branch::buildNode (Index nodeId) const & { return Node (tree_, nodeId); }
	inline Node Branch::buildNode (Index nodeId) && { return Node (std::move (tree_), nodeId); }
} // namespace Topology
} // namespace bpp

#endif // BPP_NEWPHYL_TOPOLOGY_H
