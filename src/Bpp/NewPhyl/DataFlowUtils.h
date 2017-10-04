//
// File: DataFlowUtils.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-10-04
// Last modified: 2017-10-04
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

#ifndef BPP_NEWPHYL_DATAFLOWUTILS_H
#define BPP_NEWPHYL_DATAFLOWUTILS_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Range.h> // For checkDependenciesAreAllValueNode
#include <Bpp/NewPhyl/Signed.h>
#include <cassert>
#include <typeinfo>

namespace bpp {
namespace DF {
	// Error functions / utils
	[[noreturn]] void failureDependencyTypeMismatch (const std::type_info & inNodeType,
	                                                 IndexType depIndex,
	                                                 const std::type_info & expectedType,
	                                                 const Node & givenNode);
	void checkDependencyNumber (const std::type_info & computeNodeType, SizeType expectedSize,
	                            SizeType givenSize);

	/* Node value access.
	 */
	template <typename T> bool isValueNode (const Node & n) noexcept {
		return dynamic_cast<const Value<T> *> (&n) != nullptr;
	}
	template <typename T> const T & accessValueUnsafe (const Node & n) noexcept {
		assert (isValueNode<T> (n));
		return static_cast<const Value<T> &> (n).value ();
	}

	/* Dependency type checking.
	 * In a separate namespace to not clutter the DF namespace with small names.
	 *
	 * Usage is to call check(dependencyVector, typeid(NodeType), CheckTag{});
	 * dependencyVector: NodeRefVec containing dependencies of the checked node.
	 * NodeType: type of the checked node.
	 * CheckTag: a type describing what type dependencies should have, see below.
	 */
	namespace DependencyCheck {
		// Available check tags
		template <typename T> struct AllValueNode {};      // Dynamic sized list of Value<T>
		template <typename... Types> struct ValueNodes {}; // Tuple of Value<T0>,Value<T1>, ...

		// Impl of check for AllValueNode

		template <typename T>
		void check (const NodeRefVec & dependencies, const std::type_info & inNodeType,
		            AllValueNode<T>) {
			for (auto i : index_range (dependencies)) {
				auto & dep = *dependencies[i];
				if (!isValueNode<T> (dep))
					failureDependencyTypeMismatch (inNodeType, i, typeid (Value<T>), dep);
			}
		}

		// Impl of check for ValueNodes

		inline void check (const NodeRefVec &, const std::type_info &, SizeType, ValueNodes<>) {}

		template <typename FirstType, typename... OtherTypes>
		void check (const NodeRefVec & dependencies, const std::type_info & inNodeType, SizeType index,
		            ValueNodes<FirstType, OtherTypes...>) {
			auto & dep = *dependencies[index];
			if (!isValueNode<FirstType> (dep))
				failureDependencyTypeMismatch (inNodeType, index, typeid (Value<FirstType>), dep);

			check (dependencies, inNodeType, index + 1, ValueNodes<OtherTypes...>{});
		}

		template <typename... Types>
		void check (const NodeRefVec & dependencies, const std::type_info & inNodeType,
		            ValueNodes<Types...>) {
			checkDependencyNumber (inNodeType, sizeof...(Types), dependencies.size ());
			check (dependencies, inNodeType, 0, ValueNodes<Types...>{});
		}
	}
}
}

#endif // BPP_NEWPHYL_DATAFLOWUTILS_H
