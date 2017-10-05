//
// File: DataFlowTemplateUtils.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-10-04 00:00:00
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

#ifndef BPP_NEWPHYL_DATAFLOWTEMPLATEUTILS_H
#define BPP_NEWPHYL_DATAFLOWTEMPLATEUTILS_H

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

	// Dependency check

	namespace Impl {
		// Impl of checkDependencies for ReductionOfValue
		template <typename T>
		void checkDependencies (const NodeRefVec & dependencies, const std::type_info & inNodeType,
		                        ReductionOfValue<T>) {
			for (auto i : index_range (dependencies)) {
				auto & dep = *dependencies[i];
				if (!isValueNode<T> (dep))
					failureDependencyTypeMismatch (inNodeType, i, typeid (Value<T>), dep);
			}
		}

		// Impl of checkDependencies for FunctionOfValues
		inline void checkDependenciesRecursive (const NodeRefVec &, const std::type_info &, SizeType,
		                                        FunctionOfValues<>) {}

		template <typename FirstType, typename... OtherTypes>
		void checkDependenciesRecursive (const NodeRefVec & dependencies,
		                                 const std::type_info & inNodeType, SizeType index,
		                                 FunctionOfValues<FirstType, OtherTypes...>) {
			auto & dep = *dependencies[index];
			if (!isValueNode<FirstType> (dep))
				failureDependencyTypeMismatch (inNodeType, index, typeid (Value<FirstType>), dep);

			checkDependenciesRecursive (dependencies, inNodeType, index + 1,
			                            FunctionOfValues<OtherTypes...>{});
		}

		template <typename... Types>
		void checkDependencies (const NodeRefVec & dependencies, const std::type_info & inNodeType,
		                        FunctionOfValues<Types...>) {
			checkDependencyNumber (inNodeType, sizeof...(Types), dependencies.size ());
			checkDependenciesRecursive (dependencies, inNodeType, 0, FunctionOfValues<Types...>{});
		}
	}

	/* Interface for auto generated dependency type checking.
	 * Should be used in node constructors to verify type of dependency nodes.
	 * Thow exceptions on type mismatch, with a detailed message.
	 *
	 * Manual interface:
	 * $ checkDependencies<DependencyTag> (depVector, typeid(CurrentNodeType));
	 * DependencyTag is one of the dependency structure type tags.
	 *
	 * Simplified interface:
	 * $ checkDependencies (*this);
	 * The node class must contain a typedef "Dependencies" pointing to a dependency structure
	 * type tag.
	 * Simply calls the manual interface.
	 */
	template <typename DependencyTag>
	void checkDependencies (const NodeRefVec & dependencies, const std::type_info & inNodeType) {
		Impl::checkDependencies (dependencies, inNodeType, DependencyTag{});
	}
	template <typename NodeType> void checkDependencies (const NodeType & node) {
		checkDependencies<typename NodeType::Dependencies> (node.dependencies (), typeid (NodeType));
	}

	// Wrap compute function

	namespace Impl {
		template<typename Callable, typename ... Types> void callWithValues (const NodeRefVec & dependencies, Callable && callable) {
//TODO WIP
    }
	}

	template <typename DependencyTag, typename Callable>
	void callWithValues (const NodeRefVec & dependencies, Callable && callable);
}
}
#endif // BPP_NEWPHYL_DATAFLOWTEMPLATEUTILS_H
