//
// File: Registry.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-25
// Last modified: 2017-04-25
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
#ifndef BPP_NEWPHYL_REGISTRY_H
#define BPP_NEWPHYL_REGISTRY_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Topology.h>
#include <typeindex>
#include <unordered_map>
#include <utility>

namespace bpp {
namespace DF {

	// TODO add dataset for multi model and similar
	// DataSet should know some parameters like nb site / alphabet

	class RegistryKey {
	public:
		RegistryKey (const Topology::Element & treeElement, const std::type_index & operationType)
		    : treeElement_ (treeElement), operationType_ (operationType) {}

		bool operator== (const RegistryKey & key) const noexcept {
			return treeElement_ == key.treeElement_ && operationType_ == key.operationType_;
		}
		std::size_t hashCode () const noexcept {
			auto a = treeElement_.hashCode ();
			auto b = operationType_.hash_code ();
			return a ^ (b << 1);
		}
		struct Hash {
			std::size_t operator() (const RegistryKey & key) const noexcept { return key.hashCode (); }
		};

		const Topology::Element & element () const { return treeElement_; }
		const std::type_index & operation () const { return operationType_; }

	private:
		Topology::Element treeElement_;
		std::type_index operationType_;
	};

	class Registry {
	public:
		template <typename F> Node node (const RegistryKey & key, F && createIfNotFound) {
			auto it = dataflowNodes_.find (key);
			if (it != dataflowNodes_.end ())
				return it->second;
			auto n = std::forward<F> (createIfNotFound) ();
			dataflowNodes_.emplace (key, n);
			return std::move (n);
		}
		template <typename T, typename F>
		Node node (const Topology::Element & element, F && createIfNotFound) {
			return node (RegistryKey (element, typeid (T)), std::forward<F> (createIfNotFound));
		}

		const std::unordered_map<RegistryKey, Node, RegistryKey::Hash> & rawAccess () const {
			return dataflowNodes_;
		}

	private:
		std::unordered_map<RegistryKey, Node, RegistryKey::Hash> dataflowNodes_;
	};
}
}

#endif // BPP_NEWPHYL_REGISTRY_H
