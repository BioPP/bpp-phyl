// File: DataFlowBaseClasses.h
// Authors:
//   Francois Gindraud (2017)
// Created: 21/02/2017

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

#ifndef _DF_DATAFLOW_BASE_CLASSES_H_
#define _DF_DATAFLOW_BASE_CLASSES_H_

/**
 * @file DataFlowBaseClasses.h
 * Contains the basic node classes at the root of the dataflow system.
 * These classes are mostly interfaces.
 *
 * Helper classes to define computation nodes easily are found in ComputationNodes.h
 */

#include <algorithm> // remove
#include <cassert>
#include <utility> // move, forward
#include <vector>

namespace bpp
{
  /// namespace for data flow graph (DFG) related classes.
  namespace DF
  {
    /** Data flow graph base class.
     * This class is inherited by all data flow nodes.
     * It represents (and stores data needed by) the graph of invalidations.
     * The graph of invalidations is the reverse of the value computation graph.
     * I.e., an InvalidationNode has pointers to nodes that depend on its value.
     * It also contains a boolean flag for validity.
     *
     * Another feature is a notification of death.
     * When the current node is destroyed, it will notify all dependent nodes so they can remove their links to us.
     *
     * This class is not copyable nor movable, as we will store references to it everywhere.
     * Thus all DFG classes are also non copyable nor movable.
     * TODO might change this, but maybe costly.
     *
     * This class is an internal class of the data flow system.
     * Most of its features are reserved to other DF classes usage (internal api).
     * However due to access constraints, part of the API is public (invalidate and registration).
     * Calling them outside of the DF classes is not recommended though.
     */
    class InvalidableNode
    {
    private:
      std::vector<InvalidableNode*> dependentNodes_; // TODO SSO for vectors ?
      bool valid_{false};

    protected:
      /** Callback for death notification.
       * This function will be called when a node we depend on (dependency node) is destroyed.
       * destroyedNode will contain the pointer to the destroyed dependency node.
       * A computation node should only clear depency pointers to the destroyedNode.
       */
      virtual void dependencyWasDestroyed(const InvalidableNode* destroyedNode) = 0;

      /// Set the valid flag.
      void makeValid(void) { valid_ = true; }

      /// Protected constructor.
      InvalidableNode(bool initialState)
        : valid_(initialState)
      {
      }

    public:
      // Non copyable nor movable
      InvalidableNode(const InvalidableNode&) = delete;
      InvalidableNode(InvalidableNode&&) = delete;
      InvalidableNode& operator=(const InvalidableNode&) = delete;
      InvalidableNode& operator=(InvalidableNode&&) = delete;

      /** On destruction, notify dependent nodes so they can clear dependencies to this node.
       * This prevents undetected pointers to destroyed objects.
       * Notification is done by calling the virtul dependencyWasDestroyed function.
       *
       * Dependent nodes must not deregister themselves from us (dependentNodes_ list).
       * First, it is unneeded as we are destroying the node (the list will be destroyed anyway).
       * Second, it can modify the list while we are iterating on it... bad.
       */
      virtual ~InvalidableNode()
      {
        for (auto dependentNode : dependentNodes_)
          dependentNode->dependencyWasDestroyed(this);
      }

      /** Test if the node has a valid value already computed.
       * Other nodes should not have to query the flag directly.
       * They should just use getValue from ValuedNode<T> and let the class manage the rest.
       * However, it does not hurt to let this public, and it is useful for testing (internal info).
       */
      bool isValid(void) const { return valid_; }

      /** Recursively invalidate nodes.
       * Does nothing if already invalid, as its dependent nodes will already be invalid too.
       */
      void invalidate(void)
      {
        if (isValid())
        {
          valid_ = false;
          for (auto n : dependentNodes_)
            n->invalidate();
        }
      }

      /** Adds node to the list of dependent nodes.
       * No duplicate checks performed.
       */
      void registerDependent(InvalidableNode& node) { dependentNodes_.emplace_back(&node); }

      /** Removes node from list of dependent nodes.
       * Only one instance will be removed if duplicates.
       */
      void unregisterDependent(InvalidableNode& node)
      {
        dependentNodes_.erase(std::remove(dependentNodes_.begin(), dependentNodes_.end(), &node),
                              dependentNodes_.end());
      }
    };

    /** Data flow node storing a value (interface).
     * This class is a node storing a value of type T.
     * Constructor is still protected as we are not supposed to instance this class manually.
     *
     * This class provides an interface used by all node subclasses with a value.
     * It allows a computation that requires an input value of type T to not care about how the value is provided.
     * I.e., there is no difference in interface between getValue on a parameter and a computation.
     *
     * The only indirection we need is on how to compute the value (hence the virtual compute).
     * No other function is required to be virtual in order for the data flow system to work.
     */
    template <typename T>
    class ValuedNode : public InvalidableNode
    {
    protected:
      /// Value is protected to let subclasses modify it.
      T value_;

      /** Initialize validity state, and value by forwarding arguments to constructor.
       * Allow in-place initilization of T without overhead.
       */
      template <typename... Args>
      ValuedNode(bool initialState, Args&&... args)
        : InvalidableNode(initialState)
        , value_(std::forward<Args>(args)...)
      {
      }
      /** Default constructor, with an invalid default value.
       */
      ValuedNode()
        : InvalidableNode(false)
        , value_()
      {
      }

      /** Virtual function to recompute value.
       * Computation nodes will subclass it with the computation process.
       * ParameterNode should subclass it to a failure (parameters should always be valid).
       * This function must only compute and set the new value (valid flag is managed by getValue).
       */
      virtual void compute(void) = 0;

    public:
      /** Access value, recompute if needed.
       *
       * This function is not virtual, so only 1 indirection (isValid) in the best case.
       * In the worst case, two indirection with the virtual function.
       * It would have been more natural to make this function virtual, but it would cost 2 indirections for all cases.
       */
      const T& getValue(void)
      {
        if (!isValid())
        {
          compute();
          makeValid();
        }
        return value_;
      }
    };

    /** Represents a dependency to a ValuedNode<T>.
     *
     * A Dependency is supposed to have a fixed owner node.
     * The owner node is a computation node that uses the value of the dependency.
     * The Dependency object then abstracts the (un)registering of the owner node to its producers.
     *
     * A Dependency can have two states:
     * - unset: no dependency, pointer set to nullptr.
     * - set: a dependency to the producer node pointer to by producer_; owner node is registered to producer_.
     *
     * For a reduced size, the Dependency does not store a pointer to its owner node.
     * The owner node must be passed to some API functions when needed (when it involves registering).
     * While possible, a Dependency must NOT be transfered between owner nodes (registering would not be updated).
     *
     * TODO transfer that updates registering instead of unreg+reg ?
     *
     *
     *
     * Never return a mutable ref to user code... would allow breaking the invariant.
     */
    template <typename T>
    class Dependency
    {
    public:
      using ProducerNodeType = ValuedNode<T>;

    private:
      ValuedNode<T>* producer_{};

    public:
      /// Create unset dependency.
      Dependency() = default;
      /// Create set dependency.
      Dependency(ValuedNode<T>& producer, InvalidableNode& ownerNode) { set(producer, ownerNode); }
      /** Destruction: check that dependency was unset.
       * We cannot perform unregistering without the owner node, so all dependencies must be unset at destruction.
       * The assert is used to check that in debug builds.
       */
      ~Dependency() { assert(!isSet()); }

      // Non copyable
      Dependency(const Dependency&) = delete;
      Dependency& operator=(const Dependency&) = delete;

      // Movable (between same owner node !)
      Dependency(Dependency&& other)
        : Dependency(other.producer_)
      {
        other.producer_ = nullptr;
      }
      Dependency& operator=(Dependency&& other)
      {
        assert(!isSet()); // Only movable assignable if currently unset.
        producer_ = other.producer_;
        other.producer_ = nullptr;
      }

      bool isSet(void) const { return producer_ != nullptr; }
      ValuedNode<T>* producer(void) const { return producer_; }

      const T& getValue(void)
      {
        assert(isSet());
        return producer_->getValue();
      }

      void set(ValuedNode<T>& producer, InvalidableNode& ownerNode)
      {
        assert(!isSet());
        producer_ = &producer;
        producer.registerDependent(ownerNode);
      }
      void unset(InvalidableNode& ownerNode)
      {
        assert(isSet());
        producer_->unregisterDependent(ownerNode);
        clear();
        ownerNode.invalidate(); // Node needs to be recomputed
      }
      void clear(void) { producer_ = nullptr; }
      void clearIfMatching(const InvalidableNode* producerToClear)
      {
        if (producer_ == producerToClear)
          clear();
      }
    };

    /** Node storing a parameter.
     * Leaf of the data flow dag.
     * A parameter is valid as long as it has been initialized.
     * So the compute function is a runtime error (accessing uninitialized value).
     */
    template <typename T>
    class ParameterNode : public ValuedNode<T>
    {
    private:
      /// Runtime error: should not have to recompute.
      void compute(void) override final { assert(false); }

      /// Runtime error: should never have any dependency to other nodes.
      void dependencyWasDestroyed(const InvalidableNode*) override final { assert(false); }

    public:
      /// Default constructor will be invalid state.
      ParameterNode() = default;
      /// Construct the T in place by forwarding arguments, valid state.
      template <typename... Args>
      ParameterNode(Args&&... args)
        : ValuedNode<T>(true, std::forward<Args>(args)...)
      {
      }

      /** Change parameter value.
       * For now, interface that replaces the value.
       * TODO: add in place change ?
       */
      template <typename U>
      void setValue(U&& value)
      {
        InvalidableNode::invalidate();
        ValuedNode<T>::value_ = std::forward<U>(value);
        InvalidableNode::makeValid();
      }
    };
  } // end namespace DF
} // end namespace bpp
#endif // _DF_DATAFLOW_BASE_CLASSES_H_
