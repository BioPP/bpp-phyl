//
// File: DataFlowInternal.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-10-27
// Last modified: 2017-10-27
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

#ifndef BPP_NEWPHYL_DATAFLOWINTERNAL_H
#define BPP_NEWPHYL_DATAFLOWINTERNAL_H

/** @file
 * Internal useful primitives for writing DataFlow nodes.
 *
 * Some primitives for unsafe node type cast (fast value access from base pointer).
 * checkDependencies: dependency checking guided by template types (to encode the layout).
 * callWithValues: access dependency values, forward them to a user defined lambda.
 */

#include <Bpp/NewPhyl/Cpp14.h>
#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowInternal.h>
#include <Bpp/NewPhyl/IntegerRange.h>
#include <Bpp/NewPhyl/Signed.h>
#include <algorithm>
#include <cassert>
#include <functional>
#include <type_traits>
#include <typeinfo>

namespace bpp {
namespace DF {
	/******************************* Node cast & access *******************************/

	/// Dynamic test that node inherits (or is) of NodeType.
	template <typename NodeType> bool isNodeType (const Node & n) noexcept {
		return dynamic_cast<const NodeType *> (&n) != nullptr;
	}

	/// Dynamic test that node inherits from Value<T>.
	template <typename T> bool isValueNode (const Node & n) noexcept {
		return isNodeType<Value<T>> (n);
	}

	/// Unsafe cast from Node to Value<T> (faster than dynamic cast).
	template <typename T> const Value<T> & nodeValueCast (const Node & node) noexcept {
		assert (isValueNode<T> (node));
		return static_cast<const Value<T> &> (node);
	}

	/** @brief Access maybe invalid const from raw Node &.
	 *
	 * Typically used to get matrix dimensions during DF graph transformations.
	 * Dimensions are preset in DFNumeric, so they are always valid.
	 */
	template <typename T> const T & accessValueConstCast (const Node & node) {
		return nodeValueCast<T> (node).accessValueConst ();
	}

	/// Access const with assert checking validity.
	template <typename T> const T & accessValidValueConst (const Value<T> & node) {
		assert (node.isValid ());
		return node.accessValueConst ();
	}

	/// Access valid const while casting from a raw Node &.
	template <typename T> const T & accessValidValueConstCast (const Node & node) {
		return accessValidValueConst (nodeValueCast<T> (node));
	}

	/// Access valid const while casting from NodeRef &.
	template <typename T> const T & accessValidValueConstCast (const NodeRef & nodeRef) {
		assert (nodeRef);
		return accessValidValueConstCast<T> (*nodeRef);
	}

	/******************************* Error functions *******************************/

	[[noreturn]] void failureDependencyNumberMismatch (const std::type_info & contextNodeType,
	                                                   SizeType expectedSize, SizeType givenSize);
	[[noreturn]] void failureEmptyDependency (const std::type_info & contextNodeType,
	                                          SizeType depIndex);
	[[noreturn]] void failureDependencyTypeMismatch (const std::type_info & contextNodeType,
	                                                 SizeType depIndex,
	                                                 const std::type_info & expectedType,
	                                                 const std::type_info & givenNodeType);

	/******************************* Dependency check *******************************/

	/// Checks the size of a dependency vector, throws if mismatch.
	void checkDependencyVectorSize (const std::type_info & contextNodeType, const NodeRefVec & deps,
	                                SizeType expectedSize);

	/// Checks that all dependencies are not null, throws if not.
	void checkDependenciesNotNull (const std::type_info & contextNodeType, const NodeRefVec & deps);

	/// Checks that deps[index] is a Value<T> node, throws if not.
	template <typename T>
	void checkNthDependencyIsValue (const std::type_info & contextNodeType, const NodeRefVec & deps,
	                                SizeType index) {
		const auto & dep = *deps[index];
		if (!isValueNode<T> (dep)) {
			failureDependencyTypeMismatch (contextNodeType, index, typeid (Value<T>), typeid (dep));
		}
	}

	/** @brief Check that deps is a ReductionOfValue<T> (selected by type tag).
   *
	 * A reduction is any number of Value<T> nodes.
	 */
	template <typename T>
	void checkDependencyPattern (const std::type_info & contextNodeType, const NodeRefVec & deps,
	                             ReductionOfValue<T>) {
		checkDependenciesNotNull (contextNodeType, deps);
		for (auto i : range (deps.size ()))
			checkNthDependencyIsValue<T> (contextNodeType, deps, i);
	}

	// FunctionOfValues recursion base case
	inline void checkDependencyPatternFunctionImpl (const std::type_info &, const NodeRefVec &,
	                                                SizeType, FunctionOfValues<>) {}

	// FunctionOfValues recursion iteration case
	template <typename FirstType, typename... OtherTypes>
	void checkDependencyPatternFunctionImpl (const std::type_info & contextNodeType,
	                                         const NodeRefVec & deps, SizeType index,
	                                         FunctionOfValues<FirstType, OtherTypes...>) {
		checkNthDependencyIsValue<FirstType> (contextNodeType, deps, index);
		checkDependencyPatternFunctionImpl (contextNodeType, deps, index + 1,
		                                    FunctionOfValues<OtherTypes...>{});
	}

	/** @brief Check that deps is a FunctionOfValues<Types...> (selected by type tag).
   *
	 * A function of values take Value<T> nodes with the exact types specified in the list.
	 * deps[i] must be a Value<Types[i]> node.
	 */
	template <typename... Types>
	void checkDependencyPattern (const std::type_info & contextNodeType, const NodeRefVec & deps,
	                             FunctionOfValues<Types...> tag) {
		checkDependencyVectorSize (contextNodeType, deps, sizeof...(Types));
		checkDependenciesNotNull (contextNodeType, deps);
		checkDependencyPatternFunctionImpl (contextNodeType, deps, 0, tag);
	}

	/** @brief Check that deps is an ArrayOfValues<T> (selected by type tag).
   *
	 * An array of values of size n is a reduction of fixed size.
	 */
	template <typename T>
	void checkDependencyPattern (const std::type_info & contextNodeType, const NodeRefVec & deps,
	                             ArrayOfValues<T> tag) {
		checkDependencyVectorSize (contextNodeType, deps, tag.n);
		checkDependencyPattern (contextNodeType, deps, ReductionOfValue<T>{});
	}

	/** @brief Dependency check interface: out of node class.
   *
	 * Usage: call checkDependencies<NodeType> (deps);
	 * Used to check if a dependency vector matches a pattern described by NodeType::Dependencies.
	 * Checks that dependencies match what the node excepts (number, types, non-empty).
	 * TODO revamp final interface
	 */
	template <typename NodeType> void checkDependencies (const NodeRefVec & deps) {
		checkDependencyPattern (typeid (NodeType), deps, typename NodeType::Dependencies{});
	}

	/** @brief Dependency check interface: for node constructor.
   *
	 * Usage: call checkDependencies(*this) in node constructor.
	 */
	template <typename NodeType> void checkDependencies (const NodeType & node) {
		checkDependencies<NodeType> (node.dependencies ());
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
				reduce (value, accessValidValueConstCast<ArgumentType> (*dep));
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
			    value, accessValidValueConstCast<Types> (*dependencies[Indexes])...);
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
		Impl::callWithValues (node.accessValueMutable (), node.dependencies (),
		                      typename NodeType::Dependencies{},
		                      std::forward<Callables> (callables)...);
	}

	/******************************* Optimizations *******************************/

	/** @brief Remove dependencies from the list according to a predicate.
	 *
	 * Input predicate : const NodeRef & -> bool
	 */
	template <typename Predicate>
	void removeDependenciesIf (NodeRefVec & deps, Predicate && predicate) {
		deps.erase (std::remove_if (deps.begin (), deps.end (), std::forward<Predicate> (predicate)),
		            deps.end ());
	}
} // namespace DF
} // namespace bpp

#endif // BPP_NEWPHYL_DATAFLOWINTERNAL_H
