//
// File: NodeSpecification.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-15
// Last modified: 2017-05-15
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
#ifndef BPP_NEWPHYL_NODESPECIFICATION_H
#define BPP_NEWPHYL_NODESPECIFICATION_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <memory>
#include <typeindex>
#include <unordered_map>
#include <vector>

namespace bpp {
namespace DF {

	class Registry {
	public:
		class Key {
			// Key is type of node, and dependencies
		public:
			Key (std::type_index nodeType, const std::vector<Node> & dependencies)
			    : nodeType_ (nodeType), dependencies_ (dependencies) {}

			bool operator== (const Key & other) const noexcept {
				return nodeType_ == other.nodeType_ && dependencies_ == other.dependencies_;
			}
			std::size_t hashCode () const noexcept {
				auto nodeTypeHash = std::hash<std::type_index>{}(nodeType_);
				std::size_t vecHash = dependencies_.size ();
				for (auto & n : dependencies_)
					vecHash ^= n.hashCode () + 0x9e3779b9 + (vecHash << 6) + (vecHash >> 2);
				return vecHash ^ (nodeTypeHash << 1);
			}

		private:
			std::type_index nodeType_;
			std::vector<Node> dependencies_;
		};

		// May be null
		Node get (const Key & key) const {
			auto it = nodes_.find (key);
			if (it != nodes_.end ())
				return it->second;
			else
				return Node ();
		}

		void set (Key key, Node node) {
			auto result = nodes_.emplace (std::move (key), std::move (node));
			if (!result.second)
				throw std::runtime_error ("Node already set for key");
		}

	private:
		struct Hash {
			std::size_t operator() (const Key & k) const noexcept { return k.hashCode (); }
		};
		std::unordered_map<Key, Node, Hash> nodes_;
	};

	class NodeSpecification {
		// Abstracted spec value
		// TODO doc
	public:
		// FIXME add nonself
		// TODO use perfect fwd
		template <typename T>
		explicit NodeSpecification (const T & spec) : specification_ (new Specification<T> (spec)) {}

		std::vector<NodeSpecification> computeDependencies () const {
			return specification_->computeDependencies ();
		}
		Node buildNode (std::vector<Node> dependencies) const {
			return specification_->buildNode (std::move (dependencies));
		}
		std::type_index nodeType () const { return specification_->nodeType (); }

		// Build DF graph recursively without merging
		Node instantiate () const {
			std::vector<Node> deps;
			for (auto & depSpec : computeDependencies ())
				deps.emplace_back (depSpec.instantiate ());
			return buildNode (std::move (deps));
		}

		// Build DF graph while merging using the given registry
		Node instantiateWithReuse (Registry & registry) const {
			Node n{};
			auto depSpecs = computeDependencies ();
			if (depSpecs.empty ()) {
				// Parameters are not stored in the registry.
				// buildNode should just create a new reference.
				n = buildNode ({});
			} else {
				// Instantiate dependencies
				std::vector<Node> deps;
				for (auto & depSpec : depSpecs)
					deps.emplace_back (depSpec.instantiateWithReuse (registry));
				// Check the registry
				Registry::Key key{nodeType (), deps};
				n = registry.get (key);
				if (!n.hasNode ()) {
					// Build it if not found
					n = buildNode (std::move (deps));
					registry.set (std::move (key), n);
				}
			}
			return n;
		}

	private:
		struct Interface {
			virtual ~Interface () = default;
			virtual std::vector<NodeSpecification> computeDependencies () const = 0;
			virtual Node buildNode (std::vector<Node> dependencies) const = 0;
			virtual std::type_index nodeType () const = 0;
		};
		template <typename T> struct Specification final : public Interface {
			T spec_;
			Specification (const T & spec) : spec_ (spec) {}
			std::vector<NodeSpecification> computeDependencies () const {
				return spec_.computeDependencies ();
			}
			Node buildNode (std::vector<Node> dependencies) const {
				return spec_.buildNode (std::move (dependencies));
			}
			std::type_index nodeType () const { return spec_.nodeType (); }
		};
		std::unique_ptr<Interface> specification_;
	};

	using NodeSpecificationVec = std::vector<NodeSpecification>;
}
}

#endif // BPP_NEWPHYL_NODESPECIFICATION_H
