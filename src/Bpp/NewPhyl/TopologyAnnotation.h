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
#include <memory>
#include <utility>
#include <vector>

namespace bpp {
// Forward declarations
class PhyloTree;

namespace Topology {
	// TODO reference base tree to check it is the same ?

	/* Map: associate existing DF nodes to tree elements (by id).
	 * Can be used for both Node or Branch.
	 */
	template <typename T> class BaseMap {
	public:
		explicit BaseMap (std::size_t size) : data_ (size) {}

		Optional<T> & access (IndexType id) noexcept {
			assert (id < data_.size ());
			return data_[id];
		}
		const Optional<T> & access (IndexType id) const noexcept {
			assert (id < data_.size ());
			return data_[id];
		}

	private:
		std::vector<Optional<T>> data_;
	};
	template <typename T> class NodeMap : public BaseMap<T> {
	public:
		explicit NodeMap (const std::shared_ptr<const Tree> & tree) : BaseMap<T> (tree->nbNodes ()) {}
		Optional<T> & operator[] (const Node & node) noexcept { return this->access (node.nodeId ()); }
		const Optional<T> & operator[] (const Node & node) const noexcept {
			return this->access (node.nodeId ());
		}
	};
	template <typename T> class BranchMap : public BaseMap<T> {
	public:
		explicit BranchMap (const std::shared_ptr<const Tree> & tree)
		    : BaseMap<T> (tree->nbBranches ()) {}
		Optional<T> & operator[] (const Branch & branch) noexcept {
			return this->access (branch.branchId ());
		}
		const Optional<T> & operator[] (const Branch & branch) const noexcept {
			return this->access (branch.branchId ());
		}
	};

	/* (Node|Branch)Index<T> associates T values with bijection to a tree's elements.
	 */

	// Retrieve info from PhyloTree
	struct ConvertedPhyloTreeData {
		std::shared_ptr<const Tree> topology;
		//BranchMap<double> branchLengths;
		//NodeMap<std::string> leafNames;
		// TODO add data name index, brlens map (move to other file)
	};
	ConvertedPhyloTreeData convertPhyloTree (const bpp::PhyloTree & phyloTree);
}
}

#endif // BPP_NEWPHYL_TOPOLOGYANNOTATION_H
