//
// File: TopologyMap.h
// Authors:
// Created: 2017-05-25 00:00:00
// Last modified: 2017-06-15
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

#ifndef BPP_NEWPHYL_TOPOLOGYMAP_H
#define BPP_NEWPHYL_TOPOLOGYMAP_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/FrozenPtr.h>
#include <Bpp/NewPhyl/Optional.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/NewPhyl/Topology.h>
#include <cassert>
#include <memory>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <vector>

namespace bpp {

namespace Topology {
	// Error function
	void failureIndexMapIndexAlreadySet (const std::type_info & mapType, IndexType id);
	void failureIndexMapValueAlreadySet (const std::type_info & mapType, std::string value);

	/* ValueMapBase: associate values to tree elements (by id).
	 * Base class (as both nodes and branches use IndexType indexes).
	 * Created with a fixed size, so types are not required to be movable / copyable.
	 */
	template <typename T> class ValueMapBase {
	public:
		explicit ValueMapBase (IndexType size) : data_ (size) {}

		Optional<T> & access (IndexType id) noexcept {
			assert (id < data_.size ());
			return data_[id];
		}
		const Optional<T> & access (IndexType id) const noexcept {
			assert (id < data_.size ());
			return data_[id];
		}

		using value_iterator = typename std::vector<Optional<T>>::iterator;
		using value_const_iterator = typename std::vector<Optional<T>>::const_iterator;
		Range::Range<value_iterator> value_range () noexcept { return range (data_); }
		Range::Range<value_const_iterator> value_range () const noexcept { return range (data_); }

		IndexType size () const noexcept { return IndexType (data_.size ()); }

	private:
		std::vector<Optional<T>> data_;
	};

	/* IndexMapBase: associate (both direction) values to tree elements (by id).
	 * Indexes AND values can only be set once (required to have a bijection).
	 * All set values must be different, and hashable (same system as std::unordered_map).
	 * This class is useful to add names to tree elements.
	 */
	template <typename T, typename Hash = std::hash<T>> class IndexMapBase {
	public:
		explicit IndexMapBase (IndexType size) : valueMap_ (size), indexMap_ (size) {}

		template <typename... Args> void set (IndexType id, Args &&... args) {
			auto & opt = valueMap_.access (id);
			if (opt)
				failureIndexMapIndexAlreadySet (typeid (IndexMapBase), id);
			opt.emplace (std::forward<Args> (args)...);
			auto result = indexMap_.emplace (opt.value (), id);
			if (!result.second)
				failureIndexMapValueAlreadySet (typeid (IndexMapBase), debug_to_string (opt.value ()));
		}

		const Optional<const T> & access (IndexType id) const noexcept { return valueMap_.access (id); }

		Optional<IndexType> index (const T & value) const noexcept {
			auto it = indexMap_.find (value);
			if (it != indexMap_.end ())
				return it->second;
			else
				return {};
		}

		using value_iterator = typename ValueMapBase<const T>::value_const_iterator;
		using index_iterator = typename std::unordered_map<T, IndexType, Hash>::const_iterator;
		Range::Range<value_iterator> value_range () const noexcept { return indexMap_.value_range (); }
		Range::Range<index_iterator> index_range () const noexcept { return range (indexMap_); }

		IndexType size () const noexcept { return valueMap_.size (); }

	private:
		ValueMapBase<const T> valueMap_;
		std::unordered_map<T, IndexType, Hash> indexMap_;
	};

	/* Node / Branch specific wrappers for ValueMapBase.
	 * Add overloads to build from a Tree, and access by Node / Branch iterators.
	 */
	template <typename T> class NodeValueMap : public ValueMapBase<T> {
	public:
		explicit NodeValueMap (const FrozenPtr<Tree> & tree) : ValueMapBase<T> (tree->nbNodes ()) {}
		explicit NodeValueMap (ValueMapBase<T> && map) : ValueMapBase<T> (std::move (map)) {}
		using ValueMapBase<T>::access;
		Optional<T> & access (const Node & node) noexcept { return access (node.nodeId ()); }
		const Optional<T> & access (const Node & node) const noexcept {
			return access (node.nodeId ());
		}
	};

	template <typename T> class BranchValueMap : public ValueMapBase<T> {
	public:
		explicit BranchValueMap (const FrozenPtr<Tree> & tree)
		    : ValueMapBase<T> (tree->nbBranches ()) {}
		explicit BranchValueMap (ValueMapBase<T> && map) : ValueMapBase<T> (std::move (map)) {}
		using ValueMapBase<T>::access;
		Optional<T> & access (const Branch & branch) noexcept { return access (branch.branchId ()); }
		const Optional<T> & access (const Branch & branch) const noexcept {
			return access (branch.branchId ());
		}
	};

	/* Node / Branch specific wrappers for IndexMapBase.
	 */
	template <typename T, typename Hash = std::hash<T>>
	class NodeIndexMap : public IndexMapBase<T, Hash> {
	public:
		explicit NodeIndexMap (FrozenPtr<Tree> tree)
		    : IndexMapBase<T, Hash> (tree->nbNodes ()), tree_ (std::move (tree)) {}
		explicit NodeIndexMap (FrozenPtr<Tree> tree, IndexMapBase<T, Hash> && map)
		    : IndexMapBase<T, Hash> (std::move (map)), tree_ (std::move (tree)) {
			assert (this->size () == tree_->nbNodes ());
		}

		using IndexMapBase<T, Hash>::set;
		using IndexMapBase<T, Hash>::access;
		using IndexMapBase<T, Hash>::index;
		template <typename... Args> void set (const Node & node, Args &&... args) {
			set (node.nodeId (), std::forward<Args> (args)...);
		}
		const Optional<const T> & access (const Node & node) const noexcept {
			return access (node.nodeId ());
		}
		Optional<Node> node (const T & value) const noexcept {
			return index (value).map ([this](IndexType id) { return Node{tree_, id}; });
		}
		const FrozenPtr<Tree> & tree () const noexcept { return tree_; }

	private:
		FrozenPtr<Tree> tree_;
	};

	template <typename T, typename Hash = std::hash<T>>
	class BranchIndexMap : public IndexMapBase<T, Hash> {
	public:
		explicit BranchIndexMap (FrozenPtr<Tree> tree)
		    : IndexMapBase<T, Hash> (tree->nbBranches ()), tree_ (std::move (tree)) {}
		explicit BranchIndexMap (FrozenPtr<Tree> tree, IndexMapBase<T, Hash> && map)
		    : IndexMapBase<T, Hash> (std::move (map)), tree_ (std::move (tree)) {
			assert (this->size () == tree_->nbBranches ());
		}
		using IndexMapBase<T, Hash>::set;
		using IndexMapBase<T, Hash>::access;
		using IndexMapBase<T, Hash>::index;
		template <typename... Args> void set (const Branch & branch, Args &&... args) {
			set (branch.branchId (), std::forward<Args> (args)...);
		}
		const Optional<const T> & access (const Branch & branch) const noexcept {
			return access (branch.branchId ());
		}
		Optional<Branch> branch (const T & value) const noexcept {
			return index (value).map ([this](IndexType id) { return Branch{tree_, id}; });
		}
		const FrozenPtr<Tree> & tree () const noexcept { return tree_; }

	private:
		FrozenPtr<Tree> tree_;
	};

	/* Create a ValueMap of DF::Parameters initialized by values found in another ValueMap.
	 * Parameters are only created for existing values.
	 */
	template <typename T>
	ValueMapBase<DF::Parameter<T>>
	make_parameter_map_from_value_map (const ValueMapBase<T> & valueMap) {
		ValueMapBase<DF::Parameter<T>> map (valueMap.size ());
		for (auto i : index_range (map))
			map.access (i) =
			    valueMap.access (i).map ([](const T & t) { return DF::Parameter<T>::create (t); });
		return map;
	}
	template <typename T>
	NodeValueMap<DF::Parameter<T>>
	make_node_parameter_map_from_value_map (const NodeValueMap<T> & valueMap) {
		return NodeValueMap<DF::Parameter<T>>{make_parameter_map_from_value_map (valueMap)};
	}
	template <typename T>
	BranchValueMap<DF::Parameter<T>>
	make_branch_parameter_map_from_value_map (const BranchValueMap<T> & valueMap) {
		return BranchValueMap<DF::Parameter<T>>{make_parameter_map_from_value_map (valueMap)};
	}

	/* Create a uniform map from a value, filling all possible value slots.
	 */
	template <typename T> ValueMapBase<T> make_uniform_value_map (IndexType size, const T & t) {
		ValueMapBase<T> map (size);
		for (auto & v : map.value_range ())
			v = t;
		return map;
	}
	template <typename T>
	NodeValueMap<T> make_uniform_node_value_map (const FrozenPtr<Tree> & tree, const T & t) {
		return NodeValueMap<T>{make_uniform_value_map (tree->nbNodes (), t)};
	}
	template <typename T>
	BranchValueMap<T> make_uniform_branch_value_map (const FrozenPtr<Tree> & tree, const T & t) {
		return BranchValueMap<T>{make_uniform_value_map (tree->nbBranches (), t)};
	}
}
}
#endif // BPP_NEWPHYL_TOPOLOGYMAP_H
