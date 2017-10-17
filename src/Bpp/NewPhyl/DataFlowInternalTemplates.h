//
// File: DataFlowInternalTemplates.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-10-04 00:00:00
// Last modified: 2017-10-17
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

#ifndef BPP_NEWPHYL_DATAFLOWINTERNALTEMPLATES_H
#define BPP_NEWPHYL_DATAFLOWINTERNALTEMPLATES_H

#include <Bpp/NewPhyl/Cpp14.h>
#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/NewPhyl/Signed.h>
#include <algorithm>
#include <cassert>
#include <functional>
#include <type_traits>
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
	[[noreturn]] void failureEmptyDependency (const std::type_info & inNodeType, IndexType depIndex);

	/******************************* Value<T> access *******************************/

	/// Dynamic test that node inherits from Value<T>.
	template <typename T> bool isValueNode (const Node & n) noexcept {
		return dynamic_cast<const Value<T> *> (&n) != nullptr;
	}
	/// Access the value in node assuming this is a Value<T>.
	template <typename T> const T & accessValueUncheckedCast (const Node & n) noexcept {
		assert (isValueNode<T> (n));
		return static_cast<const Value<T> &> (n).accessValue ();
	}

	/// Mutable access to a Value<T> value. Only intended for callWithValues.
	template <typename T> T & accessMutableValue (Value<T> & valueNode) noexcept {
		return valueNode.value_;
	}

	/******************************* Dependency check *******************************/

	namespace Impl {
		/* Impl of checkDependencies for ReductionOfValue.
		 * Check that each dependency is a Value<T> node.
		 */
		template <typename T>
		void checkDependencies (const NodeRefVec & dependencies, const std::type_info & inNodeType,
		                        ReductionOfValue<T>) {
			for (auto i : index_range (dependencies)) {
				auto & dep = dependencies[i];
				if (!dep)
					failureEmptyDependency (inNodeType, i);
				if (!isValueNode<T> (*dep))
					failureDependencyTypeMismatch (inNodeType, i, typeid (Value<T>), *dep);
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
			auto & dep = dependencies[index];
			if (!dep)
				failureEmptyDependency (inNodeType, index);
			if (!isValueNode<FirstType> (*dep))
				failureDependencyTypeMismatch (inNodeType, index, typeid (Value<FirstType>), *dep);

			checkDependencyTypesRecursively (dependencies, inNodeType, FunctionOfValues<OtherTypes...>{},
			                                 index + 1);
		}

		template <typename... Types>
		void checkDependencies (const NodeRefVec & dependencies, const std::type_info & inNodeType,
		                        FunctionOfValues<Types...>) {
			checkDependencyNumber (inNodeType, sizeof...(Types), dependencies.size ());
			checkDependencyTypesRecursively (dependencies, inNodeType, FunctionOfValues<Types...>{}, 0);
		}
	} // namespace Impl

	/** Dependency check interface: standalone.
	 * Usage: call checkDependencies<DependencyStructure> (deps, typeid (NodeTypeForErrorMessage));
	 * Used to check if a dependency vector matches a pattern described by DependencyStructure.
	 */
	template <typename DependencyStructure>
	void checkDependencies (const NodeRefVec & deps, const std::type_info & inNodeType) {
		Impl::checkDependencies (deps, inNodeType, DependencyStructure{});
	}

	/** Dependency check interface: for node constructor.
	 * Usage: call checkDependencies(*this) in node constructor.
	 * Checks that dependencies match what the node excepts (number, types, non-empty).
	 * The node class must contain a typedef "Dependencies" pointing to a dependency structure type
	 * tag.
	 * See tests/new_dataflow.cpp
	 */
	template <typename NodeType> void checkDependencies (const NodeType & node) {
		checkDependencies<typename NodeType::Dependencies> (node.dependencies (), typeid (NodeType));
	}

	/************************ Unpack Value<T> and call function **************************/

	namespace Impl {
		/* Impl of callWithValues for ReductionOfValue.
		 * Takes:
		 * - init(ResultType & value): sets the initial value
		 * - reduce(ResultType & value, const ArgType & arg): called for each argument
		 * TODO doc
		 * TODO alternative APIs reducing with packing
		 */
		template <typename ResultType, typename ArgumentType, typename InitFunc, typename ReduceFunc>
		void callWithValues (ResultType & value, const NodeRefVec & dependencies,
		                     ReductionOfValue<ArgumentType>, InitFunc && init, ReduceFunc reduce) {
			std::forward<InitFunc> (init) (value);
			for (const auto & dep : dependencies)
				reduce (value, accessValueUncheckedCast<ArgumentType> (*dep));
		}

		/* Impl of callWithValues for FunctionOfValues.
		 * Takes a single "function" f(ResultType & value, const T0&, const T1&, ...).
		 * TODO doc
		 */
		template <typename ResultType, typename FunctionType, typename... Types, SizeType... Indexes>
		void callWithValuesWithIndexSequence (ResultType & value, const NodeRefVec & dependencies,
		                                      FunctionOfValues<Types...>,
		                                      Cpp14::IndexSequence<Indexes...>,
		                                      FunctionType && function) {
			std::forward<FunctionType> (function) (
			    value, accessValueUncheckedCast<Types> (*dependencies[Indexes])...);
		}

		template <typename ResultType, typename FunctionType, typename... Types>
		void callWithValues (ResultType & value, const NodeRefVec & dependencies,
		                     FunctionOfValues<Types...>, FunctionType && function) {
			callWithValuesWithIndexSequence (value, dependencies, FunctionOfValues<Types...>{},
			                                 Cpp14::MakeIndexSequence<sizeof...(Types)>{},
			                                 std::forward<FunctionType> (function));
		}
	} // namespace Impl

	/** CallWithValues interface.
	 * TODO doc
	 */
	template <typename NodeType, typename... Callables>
	void callWithValues (NodeType & node, Callables &&... callables) {
		Impl::callWithValues (accessMutableValue (node), node.dependencies (),
		                      typename NodeType::Dependencies{},
		                      std::forward<Callables> (callables)...);
	}

	/******************************* Optimizations *******************************/

	/** Remove dependencies from the list according to a predicate.
	 * Input predicate : const NodeRef & -> bool
	 */
	template <typename Predicate>
	void removeDependenciesIf (NodeRefVec & deps, Predicate && predicate) {
		deps.erase (std::remove_if (deps.begin (), deps.end (), std::forward<Predicate> (predicate)),
		            deps.end ());
	}

	/** Build a predicate testing if a NodeRef is a constant value<T> matching the input predicate.
	 * Input predicate: const T & -> bool
	 * Output predicate: const NodeRef & -> bool
	 * A specialized version takes a T value and compare against it.
	 */
	template <typename T, typename Predicate,
	          typename = typename std::enable_if<!std::is_same<T, Predicate>::value>::type>
	std::function<bool(const NodeRef &)> predicateIsConstantValueMatching (Predicate predicate) {
		return [predicate](const NodeRef & nodeRef) {
			return nodeRef && nodeRef->isConstant () && isValueNode<T> (*nodeRef) &&
			       predicate (accessValueUncheckedCast<T> (*nodeRef));
		};
	}
	template <typename T>
	std::function<bool(const NodeRef &)> predicateIsConstantValueMatching (T t) {
		return predicateIsConstantValueMatching<T> ([t](const T & o) { return t == o; });
	}

} // namespace DF
} // namespace bpp
#endif // BPP_NEWPHYL_DATAFLOWINTERNALTEMPLATES_H
