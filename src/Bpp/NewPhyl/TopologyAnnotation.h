//
// File: TopologyAnnotation.h
// Authors:
// Created: 2017-05-25
// Last modified: 2017-05-25
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
#ifndef BPP_NEWPHYL_TOPOLOGYANNOTATION_H
#define BPP_NEWPHYL_TOPOLOGYANNOTATION_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Optional.h>
#include <Bpp/NewPhyl/Topology.h>
#include <cassert>
#include <utility>
#include <vector>

namespace bpp {
class PhyloTree; // Forward declaration

namespace Topology {
/* Map: associate existing DF nodes to tree elements (by id).
 * Can be used for both Node or Branch.
 */
#if 0
	class HomogeneousMap {
	public:
		HomogeneousMap (DF::Node node) : node_ (std::move (node)) {}

		DF::Node getNode (IndexType) const { return node_; }

	private:
		DF::Node node_;
	};

	class HeterogeneousMap {
	public:
		HeterogeneousMap (IndexType size) : nodeById_ (size) {}
		void setNode (IndexType id, DF::Node node) { nodeById_.at (id) = std::move (node); }

		DF::Node getNode (IndexType id) const { return nodeById_.at (id); }

	private:
		std::vector<DF::Node> nodeById_;
	};

	/* (Node|Branch)Map<T> associates T values to a tree's elements.
	 */
	class NodeMap {
	private:
		struct Interface {
			virtual ~Interface () = default;
			virtual DF::Node getNode (IndexType id) const = 0;
		};
		template <typename T> struct Map final : public Interface {
			T map_;
			Map (const T & map) : map_ (map) {}
			Map (T && map) : map_ (std::move (map)) {}
			DF::Node getNode (IndexType id) const override { return map_.getNode (id); }
		};
    std::shared_ptr<const Interface> map_;
	};
#endif

	template <typename T> class NodeMap {
	public:
		explicit NodeMap (std::shared_ptr<const Tree> tree)
		    : tree_ (std::move (tree)), data_ (tree_->nbNodes ()) {}

		bool hasParameter (const Node & node) const noexcept { return access (node).has_value (); }
		DF::Parameter<T> & getParameter (const Node & node) { return access (node).value (); }
		const DF::Parameter<T> & getParameter (const Node & node) const {
			return access (node).value ();
		}

		template <typename... Args>
		DF::Parameter<T> & createParameter (const Node & node, Args &&... args) {
			return access (node).emplace (DF::Parameter<T>::create (std::forward<Args> (args)...));
		}

	private:
		Optional<DF::Parameter<T>> & access (const Node & node) noexcept {
			assert (node.tree () == tree_);
			return data_[node.nodeId ()];
		}
		const Optional<DF::Parameter<T>> & access (const Node & node) const noexcept {
			assert (node.tree () == tree_);
			return data_[node.nodeId ()];
		}

		std::shared_ptr<const Tree> tree_;
		std::vector<Optional<DF::Parameter<T>>> data_;
	};

	/* (Node|Branch)Index<T> associates T values with bijection to a tree's elements.
	 */

	// Retrieve info from PhyloTree
	struct ConvertedPhyloTreeData {
		std::shared_ptr<const Tree> topology;
		// TODO add data name index, brlens map (move to other file)
	};
	ConvertedPhyloTreeData convertPhyloTree (const bpp::PhyloTree & phyloTree);
}
}

#endif // BPP_NEWPHYL_TOPOLOGYANNOTATION_H
