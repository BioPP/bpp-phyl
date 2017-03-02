//
// File: DataFlow.h
// Authors:
//   Francois Gindraud (2017)
// Created on: mardi 21 février 2017
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

#ifndef _DF_DATAFLOW_H_
#define _DF_DATAFLOW_H_

/**
 * @file DataFlow.h
 * Contains the basic node classes at the root of the dataflow system.
 * These classes are mostly interfaces.
 *
 * Helper classes to define computation nodes easily are found in ComputationNodes.h
 */

#if !(__cplusplus >= 201103L)
#error "Bpp/Phyl/DF/DataFlow.h requires C++11 support"
#endif

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
     * This class is not copyable nor movable, as we will store references to it everywhere.
     * Thus all DFG classes are also non copyable nor movable.
     * TODO might change this, but maybe costly.
     * TODO virtual func to notify dependentNodes_ on death ?
     *
     * This class contains implementation details, and should not be used directly.
     * Hence the constructor is protected.
     * However some of the API has to be public (registering)...
     *
     * This base class makes invalidation indendent on types, and 100% non virtual.
     * The API of this class should not be called by anything else than other DF classes.
     * However it cannot be private or protected as many DF classes rely on them.
     */
    class InvalidableNode
    {
    private:
      std::vector<InvalidableNode*> dependentNodes_; // TODO SSO for vectors ?
      bool valid_{false};

    protected:
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
      ~InvalidableNode() = default;

      /* Other nodes should not have to query the flag directly.
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
      virtual ~ValuedNode() = default;

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

    /** This class represent a consumer of a type T value generated in a ValuedNode<T>.
     * It should be used by computation nodes to represent a dependence to a value node.
     * It basically wraps a pointer to a ValuedNode<T> with a nice interface.
     * In descriptions, owner node refers to the computation node that stores and uses this dependency.
     *
     * For memory size concerns, it only stores a pointer to the producer node, and not to the owner node.
     * When (dis)connecting, it needs to register/unregister the owner node.
     * Thus the current class cannot perform (dis)connection operations by itself.
     * Instead, this class is used for storage (small) in computation nodes.
     * To perform a (dis)connection, the owner class should create a temporary DependencyConnector from this class.
     * A DependencyConnector just bundles both Dependency and owner node pointer.
     */
    template <typename T>
    class Dependency
    {
    private:
      ValuedNode<T>* producer_{nullptr};
      friend class DependencyConnector; // Allow the small connector subclass to modify the producer_.

    public:
      Dependency() = default;
      /** On destruction check that we have been disconnected properly.
       * A dependency cannot disconnect itself, as it does not store which node owns it.
       * Owner nodes must disconnect all dependencies in their destructors.
       * An assert checks that a Dependency has been disconnected before destruction in debug mode.
       */
      ~Dependency() { assert(!isConnected()); }

      bool isConnected(void) const { return producer_ != nullptr; }

      /// Get the value from the producer node linked to this dependency (check it exists).
      const T& getValue(void)
      {
        assert(isConnected());
        return producer_->getValue();
      }

      /** Little temporary class to allow connection operations.
       * Instances of this class should be built by owner nodes and returned to allow users to manage connections.
       * This class stores a reference to the owner node and to the dependency, so we can connect/disconnect.
       * It can be copied but this should not be stored...
       */
      class DependencyConnector
      {
      private:
        Dependency& dep_;
        InvalidableNode& ownerNode_;

      public:
        DependencyConnector(Dependency& dep, InvalidableNode& ownerNode)
          : dep_(dep)
          , ownerNode_(ownerNode)
        {
        }

        /// Check if we are connected.
        bool isConnected (void) const { return dep_.isConnected (); }

        /// Connects, clears any previous connection.
        void connect(ValuedNode<T>& producer)
        {
          disconnect();
          dep_.producer_ = &producer;
          dep_.producer_->registerDependent(ownerNode_);
        }

        /// Pointer overload of connect
        void connect(ValuedNode<T>* producer)
        {
          assert(producer != nullptr);
          connect(*producer);
        }

        /// Disconnect, returns true if it was connected, and invalidates its owner node.
        bool disconnect(void)
        {
          if (isConnected())
          {
            dep_.producer_->unregisterDependent(ownerNode_);
            dep_.producer_ = nullptr;
            ownerNode_.invalidate(); // Node value is now out of date
            return true;
          }
          return false;
        }
      };

      /// Builds a DependencyConnector by packing the current dependency and the owner node.
      DependencyConnector buildConnector(InvalidableNode& ownerNode) { return {*this, ownerNode}; }
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
      /// Generates a runtime error.
      void compute(void) override final { assert(false); }

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
#endif // _DF_DATAFLOW_H_
