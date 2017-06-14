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
#include <Bpp/NewPhyl/FrozenPtr.h>
#include <Bpp/NewPhyl/Optional.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/NewPhyl/Topology.h>
#include <cassert>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace bpp {
// Forward declarations
class PhyloTree;

namespace Topology {
	/* ValueMap: associate existing DF nodes to tree elements (by id).
	 * Can be used for both Node or Branch.
	 */
	template <typename T> class ValueMapBase {
	public:
		explicit ValueMapBase (std::size_t size) : data_ (size) {}

		Optional<T> & value (IndexType id) noexcept {
			assert (id < data_.size ());
			return data_[id];
		}
		const Optional<T> & value (IndexType id) const noexcept {
			assert (id < data_.size ());
			return data_[id];
		}

		IndexType size () const noexcept { return IndexType (data_.size ()); }

	private:
		std::vector<Optional<T>> data_;
	};

	template <typename T> class NodeValueMap : public ValueMapBase<T> {
	public:
		explicit NodeValueMap (const FrozenSharedPtr<Tree> & tree)
		    : ValueMapBase<T> (tree->nbNodes ()) {}
		explicit NodeValueMap (ValueMapBase<T> && map) : ValueMapBase<T> (std::move (map)) {}
		using ValueMapBase<T>::value;
		Optional<T> & value (const Node & node) noexcept { return this->value (node.nodeId ()); }
		const Optional<T> & value (const Node & node) const noexcept {
			return this->value (node.nodeId ());
		}
	};

	template <typename T> class BranchValueMap : public ValueMapBase<T> {
	public:
		explicit BranchValueMap (const FrozenSharedPtr<Tree> & tree)
		    : ValueMapBase<T> (tree->nbBranches ()) {}
		explicit BranchValueMap (ValueMapBase<T> && map) : ValueMapBase<T> (std::move (map)) {}
		using ValueMapBase<T>::value;
		Optional<T> & value (const Branch & branch) noexcept {
			return this->value (branch.branchId ());
		}
		const Optional<T> & value (const Branch & branch) const noexcept {
			return this->value (branch.branchId ());
		}
	};

	/* Utility functions to create a parameter map from a value map, with parameters initialized to
	 * the given values.
	 */
	template <typename T>
	ValueMapBase<DF::Parameter<T>> make_parameter_map_from_value_map (ValueMapBase<T> & valueMap) {
		ValueMapBase<DF::Parameter<T>> map (valueMap.size ());
		for (auto i : map.size ())
			map.value (i) =
			    valueMap.value (i).map ([](const T & t) { return DF::Parameter<T>::create (t); });
	}
	template <typename T>
	NodeValueMap<DF::Parameter<T>>
	make_node_parameter_map_from_value_map (NodeValueMap<T> & valueMap) {
		return NodeValueMap<DF::Parameter<T>>{make_parameter_map_from_value_map (valueMap)};
	}
	template <typename T>
	BranchValueMap<DF::Parameter<T>>
	make_branch_parameter_map_from_value_map (BranchValueMap<T> & valueMap) {
		return BranchValueMap<DF::Parameter<T>>{make_parameter_map_from_value_map (valueMap)};
	}

	/* (Node|Branch)Index<T> associates T values with bijection to a tree's elements.
	 * TODO
	 */

	// Retrieve info from PhyloTree
	struct ConvertedPhyloTreeData {
		FrozenSharedPtr<Tree> topology;
		FreezableUniquePtr<BranchValueMap<double>> branchLengths;
		FreezableUniquePtr<NodeValueMap<std::string>> nodeNames;
	};
	ConvertedPhyloTreeData convertPhyloTree (const bpp::PhyloTree & phyloTree);
}
}

#endif // BPP_NEWPHYL_TOPOLOGYANNOTATION_H
