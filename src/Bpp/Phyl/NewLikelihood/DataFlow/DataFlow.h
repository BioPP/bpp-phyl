//
// File: DataFlow.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-18
// Last modified: 2017-04-18
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

#ifndef BPP_DATAFLOW_H
#define BPP_DATAFLOW_H

#include <cassert>
#include <cstddef>
#include <memory>
#include <string>
#include <typeinfo>
#include <unordered_set>
#include <utility>
#include <vector>
#include <iostream>
#include <sstream>

#include <unordered_map> // For recreateWithSubstitution

#include <Bpp/Exceptions.h>

/** @file Defines the basic types of data flow nodes.
 */

namespace bpp {
/// Debug: return a readable name for a C++ type descriptor (from typeid operator).
  std::string prettyTypeName (const std::type_info & type_info);

  /// Combine hashable value to a hash, from Boost library
  template <typename T> void combineHash (std::size_t & seed, const T & t) {
    std::hash<T> hasher;
    seed ^= hasher (t) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  class Node_DF;
  template <typename T> class Value;
  class Context;


  /////  Dot output
  /// Flags to control dot output of dataflow graphs.
  enum class DotOptions {
    None = 0,
    DetailedNodeInfo = 1 << 0,
    FollowUpwardLinks = 1 << 1,
    ShowDependencyIndex = 1 << 2};
    
  DotOptions operator| (DotOptions a, DotOptions b);
  bool operator& (DotOptions a, DotOptions b);

  /// Write dataflow graph starting at nodes to output stream.
  void writeGraphToDot (std::ostream & os, const std::vector<const Node_DF *> & nodes, DotOptions opt = DotOptions::None);

  /// Write dataflow graph starting at nodes to file at filename, overwriting it.
  void writeGraphToDot (const std::string & filename, const std::vector<const Node_DF *> & nodes,
                        DotOptions opt  = DotOptions::None);

  /// Output with given precision
  template <typename T>
  std::string to_string_with_precision(const T a_value, const int n = 6)
  {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
  }
  
  ///////////////////////////////////
  /// 
  /// Node instances are always manipulated as shared pointers: provide a short alias.
    
  using NodeRef = std::shared_ptr<Node_DF>;

  /// Alias for a dependency vector (of NodeRef).
  using NodeRefVec = std::vector<NodeRef>;

  /// Shared pointer alias for Value<T>.
  template <typename T> using ValueRef = std::shared_ptr<Value<T>>;

  /// Numerical properties for DF Node.
  enum class NumericalProperty {
    Constant,         // Is a constant value
    ConstantZero,     // Is a zero value for its type
    ConstantOne,      // Is a one value for its type
    ConstantIdentity, // Is identity (for matrices mostly ; double(1.) is considered identity)
  };

  /// @name Error functions (generate a message and throw exceptions).
  ///@{
  [[noreturn]] void failureComputeWasCalled (const std::type_info & nodeType);
  [[noreturn]] void failureNodeConversion (const std::type_info & handleType, const Node_DF & node);
  [[noreturn]] void failureDependencyNumberMismatch (const std::type_info & contextNodeType,
                                                     std::size_t expectedSize, std::size_t givenSize);
  [[noreturn]] void failureEmptyDependency (const std::type_info & contextNodeType, std::size_t depIndex);
  [[noreturn]] void failureDependencyTypeMismatch (const std::type_info & contextNodeType,
                                                   std::size_t depIndex,
                                                   const std::type_info & expectedType,
                                                   const std::type_info & givenNodeType);
  ///@}

  /** @brief Base dataflow Node class.
   *
   * All data flow nodes inherit from this class.
   * Instances of this class must be used through std::shared_ptr<T>.
   *
   * A base node represents something that can be computed from its dependencies.
   * The computation is defined by the compute() method in derived classes.
   *
   * The base class stores a boolean indicating if the current value is valid.
   * This boolean is used to cache intermediate results in the graph, and avoid recomputings.
   * When data flow graph leaves change, they invalidate (transitively) dependent nodes.
   * This means setting the valid flag to false, to trigger a recomputation at next access.
   *
   * Each node has dependencies: other nodes providing data it needs to compute.
   * These dependencies are stored as NodeRef (shared_ptr<Node>).
   * Thus, multiple nodes can depend on the same value, to shared intermediate results.
   * The list of dependencies is FIXED at creation, and cannot be changed.
   * Changing the meaning of a node (dependencies, config) is done by creating a different node.
   *
   * Nodes also store pointers to dependent nodes.
   * Dependent nodes are nodes which have a direct dependency to the current node.
   * This is used to propagate invalidations to all transitive dependent nodes of a leaf.
   * Dependent nodes are stored as raw pointers (no ownership).
   *
   * Base nodes thus store the data flow graph topology, in terms of pointer to base Node.
   * This is useful for code reuse: no templates in topology walks code.
   * It is useful for debugging: can print the whole graph structure.
   * Accessing the derived classes (for compute) can be done efficiently with static_cast.
   *
   * Two main invariants are maintained at all times by this class.
   * 1: If a node is valid, all its transitive dependencies are valid.
   * 2: If a node is invalid, all transitively dependent nodes are invalid.
   *
   * In addition, two properties are implemented for performance:
   * 1: During computation of a node dependencies, do not recompute valid nodes.
   * 2: Only invalidate transitively dependent node, not others (still valid).
   *
   * Specific features are present in the base class as virtual functions.
   * This include derivation (numerical values), debug, etc.
   * These features have no-op or failure defaults which can be overriden in derived classes.
   */
    
  class Node_DF : public std::enable_shared_from_this<Node_DF> {
  public:
    Node_DF () = default;
    Node_DF (const Node_DF &) = delete;
    Node_DF (Node_DF &&) = delete;
    Node_DF & operator= (const Node_DF &) = delete;
    Node_DF & operator= (Node_DF &&) = delete;
    virtual ~Node_DF (); // Deregisters from dependencies

    // Sets dependencies + register
    Node_DF (const NodeRefVec & dependenciesArg);
    Node_DF (NodeRefVec && dependenciesArg);

    // Accessors
    bool isValid () const noexcept { return isValid_; }
    const NodeRefVec & dependencies () const noexcept { return dependencyNodes_; }
    const std::vector<Node_DF *> & dependentNodes () const noexcept { return dependentNodes_; }
    std::size_t nbDependencies () const noexcept { return dependencyNodes_.size (); }
    const NodeRef & dependency (std::size_t i) const noexcept { return dependencyNodes_[i]; }

    /// Node pretty name (default = type name).
    virtual std::string description () const;

    virtual std::string color () const
    {
      return "white";
    }

    virtual std::string shape () const
    {
      return "box";
    }

    /// Node debug info (default = ""): user defined detailed info for DF graph debug.
    virtual std::string debugInfo () const;

    /** @brief Test if the node has the given numerical property.
     *
     * This is an optional indication only, used for optimisations.
     * If unsure, leave it to always false (the default implementation).
     * This should be non recursive, to ensure a constant time check.
     */
    virtual bool hasNumericalProperty (NumericalProperty prop) const;

    /** @brief Compare node-specific configuration to another.
     *
     * Required by Context for merging.
     * It must compare everything in the node configuration except its type and dependencies.
     * See DataFlowNumeric for examples.
     * The default returns false, so nodes are considered different and not merged.
     */
    virtual bool compareAdditionalArguments (const Node_DF & other) const;
    /** @brief  Return the hash of node-specific configuration.
     *
     * Compute a hash from additional arguments of the node.
     * The hashed values must the same as the ones compared by compareAdditionalArguments.
     * The default returns 0, which is ok if there are no additional arguments.
     */
    virtual std::size_t hashAdditionalArguments () const;

    /** @brief Returns a node computing d(this_node_expression)/d(node_expression).
     *
     * The expression represented by 'node' is considered as a variable.
     * Event if 'node' is a constant value node, d(node)/d(node) == 1.
     * The derivative of a matrix is the matrix of the derivatives.
     * Derivation is undefined by default, and this function will throw an exception.
     * Implementations will usually recursively derive sub-expressions and combine them.
     */
    virtual NodeRef derive (Context & c, const Node_DF & node);

    /// Recreate the node with different dependencies.
    virtual NodeRef recreate (Context & c, NodeRefVec && deps);

    /** @brief Compute this node value, recomputing dependencies (transitively) as needed.
     *
     * Not thread safe !
     */
    void computeRecursively ();

  protected:
    /** @brief Computation implementation.
     *
     * This functions is defined by derived classes.
     * It should compute the new node value from dependency node values.
     * When called, dependency node are guaranteed to have valid values.
     *
     * This function is private to prevent use for invalid dependencies.
     * Higher level functions like computeRecursively call it while ensuring dependency validity.
     *
     * Compute has access to dependencies as a NodeRefVec (base Node_DF classes only).
     * The recommended usage is to check dependency types at Node construction.
     * Then use static_cast to access derived classes efficiently from the NodeRefVec.
     * See DataFlowNumeric.h for examples.
     */
    virtual void compute () = 0;

      
    /** @brief Invalidate (transitively) dependent nodes from this one.
     *
     * Not thread safe !
     */
    void invalidateRecursively () noexcept;

    void makeInvalid () noexcept { isValid_ = false; }
    void makeValid () noexcept { isValid_ = true; }

  protected:

    // void setDependencies_(NodeRefVec && dependenciesArg)
    // {
    //   if (dependencyNodes_.size()!=0)
    //     throw bpp::Exception("Node::setDependencies_ possible only for Nodes without dependencies.");
        
    //   dependencyNodes_=std::move(dependenciesArg);

    //   for (auto & n : dependencyNodes_)
    //       n->registerNode (this);

    //   invalidateRecursively ();
    //   makeInvalid ();
    // }

    void resetDependencies_(NodeRefVec && dependenciesArg)
    {
      for (auto & n : dependencyNodes_)
        if (n)
          n->unregisterNode (this);
        
      dependencyNodes_=dependenciesArg;

      for (auto & n : dependencyNodes_)
        if (n)
          n->registerNode (this);
        
      invalidateRecursively ();
      makeInvalid ();

    }

  private:
    void registerNode (Node_DF * n);
    void unregisterNode (const Node_DF * n);

    NodeRefVec dependencyNodes_{};         // Nodes that we depend on.
    std::vector<Node_DF *> dependentNodes_{}; // Nodes that depend on us.
    bool isValid_{false};
  };

  /// Convert a node ref with runtime type check.
  template <typename T, typename U> std::shared_ptr<T> convertRef (const std::shared_ptr<U> & from) {
    auto p = std::dynamic_pointer_cast<T> (from);
    if (!p)
      failureNodeConversion (typeid (T), *from);
    return p;
  }

  /// Check if searchedDependency if a transitive dependency of node.
  bool isTransitivelyDependentOn (const Node_DF & searchedDependency, const Node_DF & node);

  /// Recreate node by transitively replacing dependencies according to substitutions mapping.
  NodeRef recreateWithSubstitution (Context & c, const NodeRef & node,
                                    const std::unordered_map<const Node_DF *, NodeRef> & substitutions);


    
  /** @brief Abstract Node storing a value of type T.
   *
   * Represents a DataFlow node containing a T value, but still abstract (no compute()).
   * This intermediate class is used everywhere to indicate a typed dependency to a T value.
   *
   * For performance, the T value is stored directly in the node.
   * Access to the value are made without any virtual call.
   * Types can not be forward-declared.
   *
   * accessValueConst() returns the current raw value, which may be invalid.
   * getTargetValue() returns the value, recomputing it recursively if invalid.
   *
   * The Value<T> constructor forwards dependencies to the base Node, and other arguments to the
   * T value constructor.
   */
  template <typename T> class Value : public Node_DF {
  public:
    // Init deps
    template <typename... Args>
    Value (const NodeRefVec & deps, Args &&... args) : Node_DF (deps), value_ (std::forward<Args> (args)...) {}

    template <typename... Args>
    Value (NodeRefVec && deps, Args &&... args)
      : Node_DF (std::move (deps)), value_ (std::forward<Args> (args)...) {}

    /** @brief Access value, recompute if needed.
     *
     * Recompute the value if it is not up to date.
     * Then access it as const.
     * Recomputation is single threaded and not thread safe.
     */
    const T & getTargetValue () {
      this->computeRecursively ();
      return accessValueConst ();
    }

    /** @brief Raw value access (const).
     *
     * Value is not guaranteed to be valid (no recomputation).
     */
    const T & accessValueConst () const noexcept { return value_; }

    /// Derive and cast result as Value<T> (most nodes derive to the
    /// same value type).
    ValueRef<T> deriveAsValue (Context & c, const Node_DF & node) {
      return convertRef<Value<T>> (this->derive (c, node));
    }

    /** @brief General case for modification of the T object.
     *
     * Takes a callable object (lamda, function pointer) that performs the modification.
     * It must take a single T& as argument, which will refer to the T object to modify.
     * The callable is called exactly once.
     *
     *@param modifier callable to modify this object
     *@param makeValid boolean if this object is valid after
     *modification (which means no call to this->compute()
     */
      
    template <typename Callable> void modify (Callable && modifier, bool makeValid) {
      this->invalidateRecursively ();
      std::forward<Callable> (modifier) (this->accessValueMutable ());
      if (makeValid)
        this->makeValid ();
    }

  protected:
    /// Raw value access (mutable). Should only be used by
    /// subclasses to implement compute().
    T & accessValueMutable () noexcept { return value_; }

  private:
    T value_;
  };

  /// Helper: access value of Node as a Value<T> with unchecked
  /// cast.
  template <typename T> const T & accessValueConstCast (const Node_DF & node) {
    //assert (dynamic_cast<const Value<T> *> (&node) != nullptr);      // Check type in debug mode
    return static_cast<const Value<T> &> (node).accessValueConst (); // Fast cast access
  }

  /// @name Basic dependency check primitives
  ///@{
  /// Checks the size of a dependency vector, throws if mismatch.
  void checkDependencyVectorSize (const std::type_info & contextNodeType, const NodeRefVec & deps,
                                  std::size_t expectedSize);

  /// Checks the minimum size of a dependency vector, throws if mismatch.
  void checkDependencyVectorMinSize (const std::type_info & contextNodeType, const NodeRefVec & deps,
                                     std::size_t expectedMinSize);

  /// Checks that all dependencies are not null, throws if not.
  void checkDependenciesNotNull (const std::type_info & contextNodeType, const NodeRefVec & deps);

  void checkNthDependencyNotNull (const std::type_info & contextNodeType, const NodeRefVec & deps, std::size_t index);

  /// Checks that deps[index] is a T node, throws if not.
  template <typename T>
  void checkNthDependencyIs (const std::type_info & contextNodeType, const NodeRefVec & deps,
                             std::size_t index) {
    const auto & dep = *deps[index];
    if (dynamic_cast<const T *> (&dep) == nullptr)
      failureDependencyTypeMismatch (contextNodeType, index, typeid (T), typeid (dep));
  }

  /// Checks that deps[index] is a T1 or a T2 node, throws if not.
  template <typename T1, typename T2>
  void checkNthDependencyIs (const std::type_info & contextNodeType, const NodeRefVec & deps,
                             std::size_t index){
    const auto & dep = *deps[index];
    if (dynamic_cast<const T1 *> (&dep) == nullptr &&
        (dynamic_cast<const T2 *> (&dep) == nullptr))

      throw Exception (prettyTypeName (contextNodeType) + ": expected class derived from " +
                       prettyTypeName (typeid(T1))
                       + " or " + prettyTypeName (typeid(T2))
                       + " as " + std::to_string (index)
                       + "-th dependency, got " + prettyTypeName (typeid (dep)));
  }

  /// Check that deps[start, end[ contains T, throws if not.
  template <typename T>
  void checkDependencyRangeIs (const std::type_info & contextNodeType, const NodeRefVec & deps,
                               std::size_t start, std::size_t end) {
    for (std::size_t i = start; i < end; ++i) {
      checkNthDependencyIs<T> (contextNodeType, deps, i);
    }
  }
    
  /// Checks that deps[index] is a Value<T> node, throws if not.
  template <typename T>
  void checkNthDependencyIsValue (const std::type_info & contextNodeType, const NodeRefVec & deps,
                                  std::size_t index) {
    checkNthDependencyIs<Value<T>> (contextNodeType, deps, index);
  }

  /// Check that deps[start, end[ contains Value<T> nodes, throws if not.
  template <typename T>
  void checkDependencyRangeIsValue (const std::type_info & contextNodeType, const NodeRefVec & deps,
                                    std::size_t start, std::size_t end) {
    for (std::size_t i = start; i < end; ++i) {
      checkNthDependencyIsValue<T> (contextNodeType, deps, i);
    }
  }
  ///@}

  /** @brief Context for dataflow node construction.
   *
   * A context argument is passed to every function constructing dataflow nodes.
   * This class can thus be used to provide construction-time features.
   *
   * The only feature for now is merging of nodes representing the same value.
   * Each Context instance stores the set of nodes created with it.
   * When passed to a node creation function (create, derive, recreate),
   * the set is updated and used to prevent duplicate nodes.
   * Internally, creation functions will call Context::cached(newlyCreatedNode),
   * which will return an old node and discard the new if it already exists,
   * or update the set and return the new node.
   * Note that the Context owns shared_ptr to nodes,
   * so they will not be destroyed unless the context is destroyed.
   * Multiple context can be used for independent computations.
   *
   * Nodes are merged if they represent the same value.
   * As the value is not computed yet, two nodes are merged if they have:
   * - the same derived class type (same computation code),
   * - the same dependencies (same input value),
   * - same additional arguments (constants, etc).
   * Which is equivalent to the value if all additional arguments are compared.
   * (compare|hash)AdditionalArguments implement polymorphic comparison of additional arguments.
   */
  class Context {
  public:
    Context () = default;

    /** For a newly created node, return its equivalent from the cache.
     * If not already present in the cache, add it and return newNode.
     * If already present in the cache, return the stored one.
     *
     * The returned node is always of the same derived class than newNode.
     * It represents the same value.
     */
    NodeRef cached (NodeRef && newNode);

    NodeRef cached (NodeRef& newNode);

    size_t size() const
    {
      return nodeCache_.size();
    }
  private:
    /* NodeRef is hashable and comparable as a pointer.
     * CachedNodeRef is hashable and comparable, by comparing the node configuration:
     * - Derived class type,
     * - Dependencies,
     * - Additional values (for constants, Model, etc).
     * Thus a set of CachedNodeRef will merge nodes by represented value.
     */
    struct CachedNodeRef {
      NodeRef ref;

      CachedNodeRef (NodeRef && r) : ref (std::move (r)) {}

      CachedNodeRef (NodeRef & r) : ref (r) {}

      bool operator== (const CachedNodeRef & other) const;
    };
    struct CachedNodeRefHash {
      std::size_t operator() (const CachedNodeRef & ref) const;
    };
      
    std::unordered_set<CachedNodeRef, CachedNodeRefHash> nodeCache_;
  };

  /// Helper: Same as Context::cached but with a shared_ptr<T> node.
  template <typename T> std::shared_ptr<T> cachedAs (Context & c, std::shared_ptr<T> && newNode) {
    // We can use the faster static_cast due to Context::cached(): node types are conserved.
    return std::static_pointer_cast<T> (c.cached (std::move (newNode)));
  }

  template <typename T> std::shared_ptr<T> cachedAs (Context & c, std::shared_ptr<T>& newNode) {
    // We can use the faster static_cast due to Context::cached(): node types are conserved.
    return std::static_pointer_cast<T> (c.cached (newNode));
  }
} // namespace bpp

#endif // BPP_NEWPHYL_DATAFLOW_H
