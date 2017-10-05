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

#include <Bpp/NewPhyl/Cpp14.h>
#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Range.h>
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
		/* Impl of checkDependencies for ReductionOfValue.
		 * Check that each dependency is a Value<T> node.
		 */
		template <typename T>
		void checkDependencies (const NodeRefVec & dependencies, const std::type_info & inNodeType,
		                        ReductionOfValue<T>) {
			for (auto i : index_range (dependencies)) {
				auto & dep = *dependencies[i];
				if (!isValueNode<T> (dep))
					failureDependencyTypeMismatch (inNodeType, i, typeid (Value<T>), dep);
			}
		}

		/* Impl of checkDependencies for FunctionOfValues.
		 * Check that the dependency vector has the required number of dependencies.
		 * Check each dependency is a Value<Types[i]> node (recursively).
		 */
		inline void checkDependencyTypesRecursively (const NodeRefVec &, const std::type_info &,
		                                             FunctionOfValues<>, SizeType) {}

		template <typename FirstType, typename... OtherTypes>
		void checkDependencyTypesRecursively (const NodeRefVec & dependencies,
		                                      const std::type_info & inNodeType,
		                                      FunctionOfValues<FirstType, OtherTypes...>,
		                                      SizeType index) {
			auto & dep = *dependencies[index];
			if (!isValueNode<FirstType> (dep))
				failureDependencyTypeMismatch (inNodeType, index, typeid (Value<FirstType>), dep);

			checkDependencyTypesRecursively (dependencies, inNodeType, FunctionOfValues<OtherTypes...>{},
			                                 index + 1);
		}

		template <typename... Types>
		void checkDependencies (const NodeRefVec & dependencies, const std::type_info & inNodeType,
		                        FunctionOfValues<Types...>) {
			checkDependencyNumber (inNodeType, sizeof...(Types), dependencies.size ());
			checkDependencyTypesRecursively (dependencies, inNodeType, FunctionOfValues<Types...>{}, 0);
		}
	}

	/* Interface for auto generated dependency type checking.
	 * Should be used in node constructors to verify type of dependency nodes.
	 * Thow exceptions on type mismatch, with a detailed message.
	 *
	 * Manual interface:
	 * Usage: checkDependencies<DependencyTag> (depVector, typeid(CurrentNodeType));
	 * DependencyTag is one of the dependency structure type tags.
	 *
	 * Simplified interface:
	 * Usage: checkDependencies (*this);
	 * The node class (for 'this') must contain a typedef "Dependencies" pointing to a dependency
	 * structure type tag.
	 * Simply calls the manual interface.
   *
   * See tests/new_dataflow.cpp
	 */
	template <typename DependencyTag>
	void checkDependencies (const NodeRefVec & dependencies, const std::type_info & inNodeType) {
		// The DependencyTag argument is used to select the correct implementation by overloading.
		Impl::checkDependencies (dependencies, inNodeType, DependencyTag{});
	}
	template <typename NodeType> void checkDependencies (const NodeType & node) {
		checkDependencies<typename NodeType::Dependencies> (node.dependencies (), typeid (NodeType));
	}

	// Call function(s) with values from Value<T> dependencies (wrapper).

	namespace Impl {
		/* Impl of callWithValues for ReductionOfValue.
		 * TODO doc
     * TODO alternative APIs:
     * - init,reduce : init(), foreach { reduce() }
     * - init,first,reduce : if 0 { init } else { first,reduce() }
     * - consume args N by N ?
		 */
		template <typename ReduceFunctionType, typename T>
		void callWithValues (const NodeRefVec & dependencies, ReductionOfValue<T>,
		                     ReduceFunctionType reduce) {
			for (const auto & dep : dependencies)
				reduce (accessValueUnsafe<T> (*dep));
		}

		/* Impl of callWithValues for FunctionOfValues.
		 * Takes a single "function" f(const T0&, const T1&, ...).
		 * TODO doc
		 */
		template <typename FunctionType, typename... Types, SizeType... Indexes>
		void
		callWithValuesWithIndexSequence (const NodeRefVec & dependencies, FunctionOfValues<Types...>,
		                                 Cpp14::IndexSequence<Indexes...>, FunctionType && function) {
			function (accessValueUnsafe<Types> (*dependencies[Indexes])...);
		}

		template <typename FunctionType, typename... Types>
		void callWithValues (const NodeRefVec & dependencies, FunctionOfValues<Types...>,
		                     FunctionType && function) {
			callWithValuesWithIndexSequence (dependencies, FunctionOfValues<Types...>{},
			                                 Cpp14::MakeIndexSequence<sizeof...(Types)>{},
			                                 std::forward<FunctionType> (function));
		}
	}

	/* Interfaces.
	 * TODO doc
   * See tests/new_dataflow.cpp
	 */
	template <typename DependencyTag, typename... Callables>
	void callWithValues (const NodeRefVec & dependencies, Callables &&... callables) {
		Impl::callWithValues (dependencies, DependencyTag{}, std::forward<Callables> (callables)...);
	}
	template <typename NodeType, typename... Callables>
	void callWithValues (const NodeType & node, Callables &&... callables) {
		callWithValues<typename NodeType::Dependencies> (node.dependencies (),
		                                                 std::forward<Callables> (callables)...);
	}
}
}
#endif // BPP_NEWPHYL_DATAFLOWTEMPLATEUTILS_H
