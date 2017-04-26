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

#include <stdexcept>
#include <string>
#include <vector>

namespace bpp {
namespace Topology {
	using IndexType = int;
	constexpr IndexType invalid{-1};

	class Element;
	class NodeRef;
	class BranchRef;

	class Tree {
	public:
		struct Node {
			IndexType fatherId_{invalid};
			std::vector<IndexType> childrenIds_{};

			Node () = default;
		};

		std::size_t nbNodes () const { return nodes_.size (); }
		const Node & node (IndexType id) const { return nodes_[i]; }
		const IndexType & rootId () const { return rootId_; }
		IndexType & rootId () { return rootId_; }

		IndexType createNode () {
			auto id = IndexType (nbNodes ());
			nodes_.emplace (std::move (n));
			return id;
		}
		void createEdge (IndexType fatherId, IndexType childId) {
			auto & father = nodes_[fatherId];
			auto & child = nodes_[childId];
			if (child.fatherId_ != invalid)
				throw std::runtime_error ("child already has a father");
			child.fatherId_ = fatherId;
			father.childrenIds_.emplace (childId);
		}
		IndexType createNode (std::vector<IndexType> childrens) {
			auto id = createNode ();
			for (auto i : childrens)
				createEdge (id, i);
			return id;
		}

		Element nodeRef (IndexType nodeId) const;
		Element upwardBranchRef (IndexType nodeId) const;

	private:
		IndexType rootId_{invalid};
		std::vector<Node> nodes_{};
	};

	class NodeRef {
	public:
		NodeRef (const Tree & tree, IndexType nodeId) : tree_ (tree), nodeId_ (nodeId) {}

	private:
		const Tree & tree_;
		IndexType nodeId_;
	};

	class BranchRef {
	public:
		BranchRef (const Tree & tree, IndexType childNodeId)
		    : tree_ (tree), childNodeId_ (childNodeId) {}

	private:
		const Tree & tree_;
		IndexType childNodeId_;
	};

	class Element {
	public:
		enum Type { Node, Branch };
		// Able to create elements from other parts of tree, probably a polymorphic ref
		// Value type
	private:
		Type type_;
		union {
			NodeRef nodeRef;
			BranchRef branchRef;
		} d_;
	};
}
}

#endif // BPP_NEWPHYL_TOPOLOGY_H
