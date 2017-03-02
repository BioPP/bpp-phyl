//
// File: ComputationPatterns.h
// Authors:
//   Francois Gindraud (2017)
// Created on: mercredi 1 mars 2017
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _DF_COMPUTATIONNODES_H_
#define _DF_COMPUTATIONNODES_H_

/**
 * @file ComputationNodes.h
 * Computation nodes helper templates.
 * Provides classes to define common computation patterns:
 * - A fixed arguments heterogeneous computation
 * - A reduction of unknown width
 *
 * TODO fixed reduction ?
 */

#if !(__cplusplus >= 201103L)
#error "Bpp/Phyl/DF/ComputationNodes.h requires C++11 support"
#endif

#include <Bpp/Phyl/DF/Cpp14.h> // IndexSequence
#include <Bpp/Phyl/DF/DataFlow.h>
#include <algorithm> // all_of
#include <cassert>
#include <initializer_list> // foreach on tuple trick
#include <tuple>
#include <type_traits> // result_of
#include <utility>     // move, forward
#include <vector>

namespace bpp
{
  namespace DF
  {
    using bpp::Cpp14::IndexSequence;
    using bpp::Cpp14::IndexSequenceFor;

    /** Heterogeneous computation node.
     * Performs a computation with a fixed set of arguments of possibly different types.
     * This class wraps a function in a node type ready to be used.
     *
     * ArgumentTypes... represent the list of types of arguments used in the computation.
     * Callable must be a functor type (either manual, or a lambda).
     * (functor type is a type with an operator()(Args...) method).
     * The Callable function must have the following prototype:
     * @code{.cpp}
     * void func (ResultType & result, const Arg1Type & arg1, ..., const ArgNType & argN);
     * @endcode
     * With {Arg1Type, ..., ArgNType} being the types in the ArgumentTypes list.
     * The value is modified in place, to avoid reallocations if the type is a complex type like a vector.
     * Callable cannot be a function pointer, please wrap it into a FunctionPointerWrapper class.
     *
     * Each argument generates a dependency (something requesting a value).
     * All dependencies must be connected to producer nodes before we can compute the current node.
     * During computation, arguments values will be retrieved from the producers nodes by following the dependencies.
     * dependency<N> () returns a DependencyConnector object that can be used to connect/disconnect the N-th dependency.
     * 
     * An example of use is available in the test/test_dataflow.cpp file.
     *
     * EBO:
     * In most cases, compute functions will not use any data except their arguments (pure functions).
     * Thus the functors types will have a 0 size.
     * If the functor is an attribute, the compiler must make it 1 byte at least for reasons.
     * This loses space, more than 1 byte due to alignment.
     * Thus we use private inheritance, which lets the compiler make it really 0 size (empty base optimisation).
     * The only downside is that Callable must be a class type.
     * Function pointers are not class types, so we cannot use them.
     * FunctionPointerWrapper allow to use function pointers and have static linking to the function.
     */
    template <typename T, typename Callable, typename... ArgumentTypes>
    class HeterogeneousComputationNode : public ValuedNode<T>, private Callable
    {
    private:
      /// List of dependencies to nodes
      std::tuple<Dependency<ArgumentTypes>...> dependencies_;

    public:
      /// Gives the type of the N-th dependency.
      template <std::size_t Index>
      using DependencyType = typename std::tuple_element<Index, decltype(dependencies_)>::type;

    private:
      /* Helper function: only useful to unpack the index sequence.
       * Unpacking the sequence is needed to apply std::get to the tuple an call getValue on each dependency.
       * For more information, see cppreference std::tuple and std::integer_sequence pages.
       */
      template <std::size_t... Is>
      void computeWithUnpackedIndexSequence(IndexSequence<Is...>)
      {
        // Forward call to functor
        Callable::operator()(ValuedNode<T>::value_, std::get<Is>(dependencies_).getValue()...);
      }
      void compute(void) override final { computeWithUnpackedIndexSequence(IndexSequenceFor<ArgumentTypes...>{}); }

      // Helper function for destructor
      template <std::size_t... Is>
      void destructorWithUnpackedIndexSequence(IndexSequence<Is...>)
      {
        // Awful trick to force calling disconnect on all elements of the tuple.
        // Only a few context allow a parameter pack expansion, and struct initializer list are one of them.
        (void)std::initializer_list<bool>{dependency<Is>().disconnect()...};
      }

    public:
      /** Constructor.
       * A value of the callable must be provided first.
       * This is required for lambdas as they have no default constructor.
       * The value is also initialized using the remaining arguments (but set to invalid to trigger the first computation).
       * This initialization is useful for complex types (like T = std::vector<...>), to set their initial working size.
       */
      template <typename C, typename... Args>
      HeterogeneousComputationNode(C&& c, Args&&... args)
        : ValuedNode<T>(false, std::forward<Args>(args)...)
        , Callable(std::forward<C>(c))
      {
      }

      /** Destructor.
       * Disconnects all dependences.
       */
      ~HeterogeneousComputationNode() { destructorWithUnpackedIndexSequence(IndexSequenceFor<ArgumentTypes...>{}); }

      /** Gets DependencyConnector object linked to the n-th dependency object.
       * This object allows to connect() or disconnect() to any ValuedNode<T>.
       * All dependencies must be connected to a producer before computing this node value.
       * An unconnected dependency during the computation of this node is a runtime error.
       */
      template <std::size_t Index>
      typename DependencyType<Index>::DependencyConnector dependency(void)
      {
        return std::get<Index>(dependencies_).buildConnector(*this);
      }

      /// Get number of dependencies (static function).
      static constexpr std::size_t nbDependencies(void) { return sizeof...(ArgumentTypes); }
    };

    /** Perform a reduction with an arbitrary number of elements.
     * Each element is represented by a dependency.
     * Contrary to HeterogeneousComputationNode, every dependency of this class must be connected.
     * Thus removing an element from the computation is equivalent to removing the dependency completely.
     *
     * The reduction operation is defined by implementing a specific struct that defines how to perform the reduction.
     * Any reduction defining struct must have the following elements:
     * @code{.cpp}
     * struct MyReduction {
     *   using ResultType = ...; // the type of the result
     *   using ArgumentType = ...; // the type of the arguments
     *   static void reset (ResultType &); // a function that resets the storage
     *   static void reduce (ResultType & accumulator, const ArgumentType & arg); // one step of the reduction
     * };
     * @endcode
     */
    template <typename ReductionOp>
    class ReductionComputationNode : public ValuedNode<typename ReductionOp::ResultType>
    {
    public:
      using ResultType = typename ReductionOp::ResultType;
      using ArgumentType = typename ResultType::ArgumentType;

    private:
      std::vector<Dependency<ArgumentType>> dependencies_;
      using ValuedNode<ResultType>::value_;

      void compute(void) override final
      {
        ReductionOp::reset(value_);
        for (auto& dep : dependencies_)
          ReductionOp::reduce(value_, dep.getValue());
      }

    public:
      /// Dependencies manual disconnect.
      ~ReductionComputationNode()
      {
        for (auto& dep : dependencies_)
          dep.buildConnector(*this).disconnect();
      }

      void addDependencyTo(ValuedNode<ArgumentType>& producer)
      {
        dependencies_.emplace_back();
        dependencies_.back().buildConnector(*this).connect(producer);
      }

      // TODO removal and other manipulation tools

      std::size_t nbDependencies(void) const { return dependencies_.size(); }
    };

    /** Creates a functor out of a compile time function pointer.
     * This class is used to give a raw function pointer to the computation classes.
     * Usage:
     * @code{.cpp}
     * void f (int & r, int a, int b) { r = a + b; }
     * using AddIntFunctor = FunctionPointerWrapper<decltype (&f), f>;
     * HeterogeneousComputationNode<int, Functor, int, int> node (AddIntFunctor {});
     * @endcode
     *
     * Most of the time, we know at compile time which function should be used in a node instance.
     * For performance it is better for the node to generate a static call instead of an indirect call.
     * Storing a function pointer generates an indirect call.
     * Thus we use this type which has no storage size, and allow to statically link to the given function.
     */
    template <typename FuncPtrType, FuncPtrType funcPtr>
    struct FunctionPointerWrapper
    {
      template <typename... Args>
      void operator()(Args&&... args)
      {
        static_assert(std::is_same<void, typename std::result_of<FuncPtrType(Args...)>::type>::value,
                      "bpp::DF::FunctionPointerWrapper: function must return void");
        funcPtr(std::forward<Args>(args)...);
      }
    };
  } // end namespace DF
} // end namespace bpp
#endif // _DATAFLOW_H_
