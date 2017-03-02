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
    using bpp::Cpp14::MakeIndexSequence;

    /** Heterogeneous computation node.
     * Performs a computation with a fixed set of arguments of heterogeneous types.
     * The computation itself must be encoded as s struct defining multiple elements:
     * @code{.cpp}
     * struct ComputationStruct {
     *    using ResultType = ...; // return type of the computation
     *    using ArgumentTypes = std::tuple<Arg1Type, ..., ArgNType>; // list of types of the arguments (in order)
     *    // Computation function
     *    static void compute (ResultType & result, const Arg1Type & arg1, ..., const ArgNType & argN);
     * };
     * @endcode
     *
     * For each argument type given, the node class will define a Dependency class instance.
     * This class represents an input to the computation.
     * Each dependency (Dependency class instance) must be connected to a ValuedNode of the matching type.
     * When computing our own value, values of the dependent nodes will be computed and then fed to the compute function.
     *
     * A computation attempt while some dependencies are not connected is a runtime error.
     * dependency<N> () returns a DependencyConnector object that can be used to connect/disconnect the N-th dependency.
     * 
     * An example of use is available in the test/test_dataflow.cpp file.
     */
    template <typename ComputationOp>
    class HeterogeneousComputationNode : public ValuedNode<typename ComputationOp::ResultType>
    {
    private:
      // This struct is used to generate a std::tuple<Dependency<Arg1Type>, ..., std::tuple<Dependency<ArgNType>>
      template <typename T>
      struct MakeDependencyTuple;
      template <typename... ArgumentTypes>
      struct MakeDependencyTuple<std::tuple<ArgumentTypes...>>
      {
        using Type = std::tuple<Dependency<ArgumentTypes>...>;
      };
      /// Generated type of a tuple containing Dependency<Arg> structs for each argument type in order.
      using DependencyTupleType = typename MakeDependencyTuple<typename ComputationOp::ArgumentTypes>::Type;

    public:
      using ResultType = typename ComputationOp::ResultType;

      /// Gives the type of the N-th dependency.
      template <std::size_t Index>
      using DependencyType = typename std::tuple_element<Index, DependencyTupleType>::type;

    private:
      /// Static heterogeneous list of dependencies.
      DependencyTupleType dependencies_;

    public:
      /** Constructor.
       * The stored value is initialized using the given arguments (forwarded to ResultType constructor).
       * This initialization step is useful for complex types like std::vector<...>.
       * Using this, they can be allocated only once to their working size instead of every computation.
       */
      template <typename... Args>
      HeterogeneousComputationNode(Args&&... args)
        : ValuedNode<ResultType>(false, std::forward<Args>(args)...)
      {
      }

      /** Destructor.
       * Disconnects all dependences, as required by the Dependency class.
       * Use a helper function to unpack the tuple (see std::index_sequence).
       */
      ~HeterogeneousComputationNode() { destructorWithUnpackedIndexSequence(MakeIndexSequence<nbDependencies()>{}); }

      /// Get number of dependencies (static function).
      static constexpr std::size_t nbDependencies(void) { return std::tuple_size<DependencyTupleType>::value; }

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

    private:
      /* Helper function: only useful to unpack the index sequence.
       * Unpacking the sequence is needed to apply std::get to the tuple an call getValue on each dependency.
       * For more information, see cppreference std::tuple and std::integer_sequence pages.
       */
      template <std::size_t... Is>
      void computeWithUnpackedIndexSequence(IndexSequence<Is...>)
      {
        ComputationOp::compute(ValuedNode<ResultType>::value_, std::get<Is>(dependencies_).getValue()...);
      }

      /** Compute implementation.
       * Apply the compute function on values retrieved from dependent nodes.
       * Use a helper function to unpack the tuple (see std::index_sequence).
       */
      void compute(void) override final { computeWithUnpackedIndexSequence(MakeIndexSequence<nbDependencies()>{}); }

      // Helper function for destructor
      template <std::size_t... Is>
      void destructorWithUnpackedIndexSequence(IndexSequence<Is...>)
      {
        // Awful trick to force calling disconnect on all elements of the tuple.
        // Only a few context allow a parameter pack expansion, and struct initializer list are one of them.
        (void)std::initializer_list<bool>{dependency<Is>().disconnect()...};
      }
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
      /// Dynamic list of dependencies.
      std::vector<Dependency<ArgumentType>> dependencies_;

    public:
      /** Constructor.
       * The stored value is initialized using the given arguments (forwarded to ResultType constructor).
       * This initialization step is useful for complex types like std::vector<...>.
       * Using this, they can be allocated only once to their working size instead of every computation.
       */
      template <typename... Args>
      ReductionComputationNode(Args&&... args)
        : ValuedNode<ResultType>(false, std::forward<Args>(args)...)
      {
      }

      /** Destructor.
       * Disconnects all dependences, as required by the Dependency class.
       */
      ~ReductionComputationNode()
      {
        for (auto& dep : dependencies_)
          dep.buildConnector(*this).disconnect();
      }

      std::size_t nbDependencies(void) const { return dependencies_.size(); }

      void addDependencyTo(ValuedNode<ArgumentType>& producer)
      {
        dependencies_.emplace_back();
        dependencies_.back().buildConnector(*this).connect(producer);
      }

      // TODO removal and other manipulation tools

    private:
      using ValuedNode<ResultType>::value_;

      /// Compute implementation; applies reduce iteratively for each dependency.
      void compute(void) override final
      {
        ReductionOp::reset(value_);
        for (auto& dep : dependencies_)
          ReductionOp::reduce(value_, dep.getValue());
      }
    };
  } // end namespace DF
} // end namespace bpp
#endif // _DATAFLOW_H_
