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
#include <vector>

namespace bpp {
namespace DF {

	class NodeSpecification {
		// Abstracted spec value
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

    // Build DF graph recursively without merging
		Node instantiate () const {
			std::vector<Node> deps;
			for (auto & depSpec : computeDependencies ())
				deps.emplace_back (depSpec.instantiate ());
			return buildNode (std::move (deps));
		}

		// TODO add registry version
		// Registry should now be a map from (type, deps) -> node
		// Should only depend on the DF graph and not the spec

	private:
		struct Interface {
			virtual ~Interface () = default;
			virtual std::vector<NodeSpecification> computeDependencies () const = 0;
			virtual Node buildNode (std::vector<Node> dependencies) const = 0;
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
		};
		std::unique_ptr<Interface> specification_;
	};

	using NodeSpecificationVec = std::vector<NodeSpecification>;
}
}

#endif // BPP_NEWPHYL_NODESPECIFICATION_H
