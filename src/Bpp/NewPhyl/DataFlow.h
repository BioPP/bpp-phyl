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

#ifndef BPP_NEWPHYL_DATAFLOW_H
#define BPP_NEWPHYL_DATAFLOW_H

#include <Bpp/NewPhyl/Dimension.h>
#include <Bpp/NewPhyl/Signed.h>
#include <Bpp/NewPhyl/Vector.h>
#include <map>
#include <memory>
#include <string>   // description
#include <typeinfo> // convertRef
#include <utility>

/** @file Defines the basic types of data flow nodes.
 */
namespace bpp {
namespace DF {
	class Node;
	template <typename T> class Value;

	/// Node instances are always manipulated as shared pointers, use an alias.
	using NodeRef = std::shared_ptr<Node>;
	/// Shared pointer alias for Value<T>.
	template <typename T> using ValueRef = std::shared_ptr<Value<T>>;

	/// Alias for a dependency vector (of NodeRef).
	using NodeRefVec = Vector<NodeRef>;

	/** Base data flow Node class.
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
	 * Specific features are present in the base class as virtual functions.
	 * This include derivation (numerical values), information (debug, name, isConstant...).
	 * These features have no-op defaults which can be overriden if meaningful in derived classes.
	 */
	class Node {
	public:
		Node () = default;
		Node (const Node &) = delete;
		Node (Node &&) = delete;
		Node & operator= (const Node &) = delete;
		Node & operator= (Node &&) = delete;
		virtual ~Node (); // Deregisters from dependencies

		// Sets dependencies + register
		Node (const NodeRefVec & dependenciesArg);
		Node (NodeRefVec && dependenciesArg);

		// Accessors
		bool isValid () const noexcept { return isValid_; }
		const NodeRefVec & dependencies () const noexcept { return dependencyNodes_; }
		const Vector<Node *> & dependentNodes () const noexcept { return dependentNodes_; }
		SizeType nbDependencies () const noexcept { return dependencyNodes_.size (); }
		const NodeRef & dependency (SizeType i) const noexcept { return dependencyNodes_[i]; }

		/// Node string description (default = type name): shorter name for in debug.
		virtual std::string description () const;

		/// Node debug info (default = ""): user defined detailed info for DF graph debug.
		virtual std::string debugInfo () const;

		/** Indicates if the node represents a constant value.
     * This is an optional indication only, used for optimisations.
     * If unsure, leave it to false (the default implementation).
     * This should be non recursive (only test the local node, not the sub tree).
     */
		virtual bool isConstant () const;

    /** Is the node a zero constant of its type (indication, like isConstant).
     * For vector / matrices, this indicates a constant filled with zeroes.
     */
    virtual bool isConstantZero () const;

    /** Is the node a zero constant of its type (indication, like isConstant).
     * For vector / matrices, this indicates a constant filled with ones.
     */
    virtual bool isConstantOne () const;

		/// Derive with respect to node (default = not implemented), recursive.
		virtual NodeRef derive (const Node & node);
		virtual bool isDerivable (const Node & node) const;
		// TODO derivation doc somewhere

		// Check if node is part of transitive dependencies
		bool isTransitivelyDependentOn (const Node & node) const;

		// Rebuild the node with different dependencies
		virtual NodeRef rebuild (NodeRefVec && deps) const;

	protected:
		/** Computation implementation.
		 * This functions is defined by derived classes.
		 * It should compute the new node value from dependency node values.
		 * When called, dependency node are guaranteed to have valid values.
		 *
		 * This function is private to prevent use for invalid dependencies.
		 * Higher level functions like computeRecursively() call it while ensuring dependency validity.
     *
     * Compute has access to dependencies as a NodeRefVec (base Node classes only).
     * The recommended usage is to check dependency types at Node construction.
     * Then use static_cast to access derived classes efficiently from the NodeRefVec.
     * Several helper functions in DataFlowInternal.h simplify this.
     * See DataFlowNumeric.h for examples.
		 */
		virtual void compute () = 0;

    /** Invalidate (transitively) dependent nodes from this one.
     * Not thread safe !
     */
		void invalidateRecursively () noexcept;

    /** Compute this node value, recomputing dependencies (transitively) as needed.
     * Not thread safe !
     */
		void computeRecursively ();

		void makeInvalid () noexcept { isValid_ = false; }
		void makeValid () noexcept { isValid_ = true; }

	private:
		void registerNode (Node * n);
		void unregisterNode (const Node * n);

	private:
		NodeRefVec dependencyNodes_{};    // Nodes that we depend on.
		Vector<Node *> dependentNodes_{}; // Nodes that depend on us.
		bool isValid_{false};
	};

	/* Free functions.
	 */
	NodeRef rebuildWithSubstitution (const NodeRef & node,
	                                 const std::map<const Node *, NodeRef> & substitutions);

	/** Node construction class.
	 * Node are by convention manipulated as shared_ptr<NodeType>.
	 * Node construction should always be done through the Builder<NodeType>::make function.
	 * This class defaults to simply using std::make_shared.
	 * However it can be specialised for some node types.
	 * Specialisations may perform merging, numeric optimisation / simplification.
	 *
	 * The shorter makeNode<NodeType> (...) function is provided (it calls Builder<NodeType>::make).
	 */
	template <typename NodeType> struct Builder {
		template <typename... Args> static std::shared_ptr<NodeType> make (Args &&... args) {
			return std::make_shared<NodeType> (std::forward<Args> (args)...);
		}
	};

	/// Defers to the Builder<NodeType>::make.
	template <typename NodeType, typename... Args>
	auto makeNode (Args &&... args)
	    -> decltype (Builder<NodeType>::make (std::forward<Args> (args)...)) {
		return Builder<NodeType>::make (std::forward<Args> (args)...);
	}
	/// Overload that accepts a NodeRef initializer list (commonly used).
	template <typename NodeType, typename... Args>
	auto makeNode (std::initializer_list<NodeRef> deps, Args &&... args)
	    -> decltype (Builder<NodeType>::make (deps, std::forward<Args> (args)...)) {
		return Builder<NodeType>::make (deps, std::forward<Args> (args)...);
	}

	// Used in Value<T> to select constructors
	struct NoDependencyTag {
		constexpr NoDependencyTag () = default;
	};
	constexpr NoDependencyTag noDependency{};

	/* Valued node. TODO doc
	 * Represents a DataFlow node containing a T value, but still abstract (no compute()).
	 */
	template <typename T> class Value : public Node {
	public:
		// Init deps
		template <typename... Args>
		Value (const NodeRefVec & deps, Args &&... args)
		    : Node (deps), value_ (std::forward<Args> (args)...) {}
		template <typename... Args>
		Value (NodeRefVec && deps, Args &&... args)
		    : Node (std::move (deps)), value_ (std::forward<Args> (args)...) {}

		// Without dependencies
		template <typename... Args>
		Value (NoDependencyTag, Args &&... args) : Node (), value_ (std::forward<Args> (args)...) {}

		/** Access value, recompute if needed.
		 * Recompute the value if it is not up to date.
		 * Then access it as const.
		 * Recomputation is single threaded.
		 */
		const T & getValue () {
			this->computeRecursively ();
			return accessValueConst ();
		}

		/** Raw value access (const).
		 * Value is not guaranteed to be valid (no recomputation).
		 */
		const T & accessValueConst () const noexcept { return value_; }

		/** Raw value access (mutable).
		 * Access the value as a mutable reference (may not be valid, no recomputation).
		 * Modifying a computed value through this functions does NOT invalidate dependent values.
		 * This breaks data flow invariants, so only use it in internal data flow node code.
		 */
		T & accessValueMutable () noexcept { return value_; }

		/** Access the node target dimension.
		 * Target dimension is a node property, useful for sized data types (vector, matrices).
		 * It is used as documentation of the node.
		 * Target dimension can be compared to actual dimension (after computation) as a safety check.
		 * It is also used for numerical graph simplification (creating simplified node of same size).
		 * This property should be set by using setTargetDimension.
		 */
		const Dimension<T> & getTargetDimension () const noexcept { return targetDimension_; }

		/// Set the target dimension of the node.
		void setTargetDimension (const Dimension<T> & dim) { targetDimension_ = dim; }

		// Defined as default to enable specialisation
		std::string debugInfo () const override { return Node::debugInfo (); }

	protected:
		Dimension<T> targetDimension_{};
		T value_;
	};

	// Debug info override for double (in DataFlowNumeric.cpp)
	template <> std::string Value<double>::debugInfo () const;
	// overrides for VectorDouble and MatrixDouble are declared in LinearAlgebra.h

	/* Dependency structure description.
	 * These type tags are used to specify compute node dependency types.
	 * This can serve as documentation about what arguments node expect.
	 * Helper functions in DataFlowInternal.h act depending on these type tags.
	 */
	template <typename T> struct ReductionOfValue {};           // Dynamic sized list of Value<T>
	template <typename... Types> struct FunctionOfValues {};    // Tuple of Value<T0>, Value<T1>, ...
	template <typename T> struct ArrayOfValues { SizeType n; }; // n Value<T>

	// Error function
	[[noreturn]] void failureNodeConversion (const std::type_info & handleType, const Node & node);

	// Convert handles with check
	template <typename T, typename U>
	std::shared_ptr<T> convertRef (const std::shared_ptr<U> & from) {
		auto p = std::dynamic_pointer_cast<T> (from);
		if (!p)
			failureNodeConversion (typeid (T), *from);
		return p;
	}
} // namespace DF
} // namespace bpp

#endif // BPP_NEWPHYL_DATAFLOW_H
