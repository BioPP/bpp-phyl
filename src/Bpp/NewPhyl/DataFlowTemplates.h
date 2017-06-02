//
// File: DataFlowTemplates.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-09 00:00:00
// Last modified: 2017-05-09
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
#ifndef BPP_NEWPHYL_DATAFLOWTEMPLATES_H
#define BPP_NEWPHYL_DATAFLOWTEMPLATES_H

#include <Bpp/NewPhyl/Cpp14.h>
#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Range.h>
#include <tuple>
#include <typeinfo>

namespace bpp {
namespace DF {
	// Error utils
	void checkDependencyNumber (const std::type_info & computeNodeType, std::size_t expectedSize,
	                            std::size_t givenSize);
	void failureDependencyTypeMismatch (const std::type_info & computeNodeType, std::size_t depIndex,
	                                    const std::type_info & expectedType,
	                                    const Node::Impl & givenNode);

	/** Generic function computation.
	 * Performs a computation with a fixed set of arguments of heterogeneous types.
	 * The computation itself must be encoded as a struct defining multiple elements:
	 * @code{.cpp}
	 * struct MyComputationStruct {
	 *    using ResultType = ...; // return type of the computation
	 *    using ArgumentTypes = std::tuple<Arg1Type, ..., ArgNType>; // list of types of the arguments
	 *    // Computation function
	 *    static void compute (ResultType & result, const Arg1Type & a1, ..., const ArgNType & aN);
	 * };
	 * @endcode
	 *
	 * For each argument type given, a dependency slot is reserved.
	 * Each dependency must be connected to a Value<T> of the matching type.
	 * Dependencies given to the constructor must match.
	 */
	template <typename Op>
	class GenericFunctionComputation : public Value<typename Op::ResultType>::Impl {
	public:
		using ResultType = typename Op::ResultType;
		using ArgumentTypes = typename Op::ArgumentTypes;
		static constexpr auto nbDependencies = std::tuple_size<ArgumentTypes>::value;
		template <std::size_t Index>
		using ArgumentType = typename std::tuple_element<Index, ArgumentTypes>::type;

		template <typename... Args>
		GenericFunctionComputation (NodeVec deps, Args &&... args)
		    : Value<ResultType>::Impl (std::move (deps), std::forward<Args> (args)...) {
			checkDependencies (Cpp14::MakeIndexSequence<nbDependencies>{});
		}

	private:
		/** Runtime checks the dependency for type / number mismatch.
		 */
		template <std::size_t... Is> void checkDependencies (Cpp14::IndexSequence<Is...>) {
			checkDependencyNumber (typeid (GenericFunctionComputation), nbDependencies,
			                       this->dependencies ().size ());
			// check dependency type for each required argument (uses IndexSequence for tuple unpacking)
			static_cast<void> (std::initializer_list<int>{checkDependencyType<Is> ()...});
		}
		template <std::size_t Index> int checkDependencyType () {
			using ArgType = ArgumentType<Index>;
			if (!isValueNode<ArgType> (this->dependencies ()[Index]))
				failureDependencyTypeMismatch (typeid (GenericFunctionComputation), Index,
				                               typeid (typename Value<ArgType>::Impl),
				                               this->dependencies ()[Index].getImpl ());
			return 0; // Dummy int for the initializer_list<int>
		}

		/** Compute implementation.
		 * Apply the compute function on values retrieved from dependent nodes.
		 */
		void compute () override final {
			computeWithUnpackedIndexSequence (Cpp14::MakeIndexSequence<nbDependencies>{});
		}
		template <std::size_t... Is>
		void computeWithUnpackedIndexSequence (Cpp14::IndexSequence<Is...>) {
			// Helper function : unpack and cast each dependency to its type, feed values to compute.
			Op::compute (this->value_, getValueUnsafe<ArgumentType<Is>> (this->dependencies ()[Is])...);
		}
	};

	/** Generic reduction computation.
	 * Performs a reduction by applying a reduction function over an arbitrary set of arguments.
	 * Arguments must all have the same type.
	 *
	 * The reduction operation is defined by implementing a specific struct that defines how to
	 * perform the reduction.
	 * Any reduction defining struct must have the following elements:
	 * @code{.cpp}
	 * struct MyReductionStruct {
	 *   using ResultType = ...; // the type of the result
	 *   using ArgumentType = ...; // the type of the arguments
	 *   static void reset (ResultType &); // a function that resets the storage
	 *   // one step of the reduction
	 *   static void reduce (ResultType & accumulator, const ArgumentType & arg);
	 * };
	 * @endcode
	 */
	template <typename Op>
	class GenericReductionComputation : public Value<typename Op::ResultType>::Impl {
	public:
		using ResultType = typename Op::ResultType;
		using ArgumentType = typename Op::ArgumentType;

		template <typename... Args>
		GenericReductionComputation (NodeVec deps, Args &&... args)
		    : Value<ResultType>::Impl (std::move (deps), std::forward<Args> (args)...) {
			for (auto i : index_range (this->dependencies ()))
				if (!isValueNode<ArgumentType> (this->dependencies ()[i]))
					failureDependencyTypeMismatch (typeid (GenericReductionComputation), i,
					                               typeid (typename Value<ArgumentType>::Impl),
					                               this->dependencies ()[i].getImpl ());
		}

	private:
		void compute () override final {
			Op::reset (this->value_);
			for (const auto & dep : this->dependencies ())
				Op::reduce (this->value_, getValueUnsafe<ArgumentType> (dep));
		}
	};
}
}
#endif // BPP_NEWPHYL_DATAFLOWTEMPLATES_H
