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

#include <Bpp/NewPhyl/Signed.h>
#include <Bpp/NewPhyl/Vector.h>
#include <memory>
#include <string>   // description
#include <typeinfo> // convertRef
#include <utility>

namespace bpp {
namespace DF {
	/* In this file:
	 * - Basic node class definition.
	 */

	// Declarations
	class Node;
	template <typename T> class Value;
	class InternalAccessor;
	using NodeRef = std::shared_ptr<Node>;
	using NodeRefVec = Vector<NodeRef>;
	template <typename T> using ValueRef = std::shared_ptr<Value<T>>;

	/* Base Node class.
	 * Abstract : compute() needs to be defined to the actual computation.
	 * It is supposed to be used with std::shared_ptr for node class ownership.
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

		// Node string description (default = type name): shorter name for in debug.
		virtual std::string description () const;

		// Node debug info (default = ""): user defined detailed info for DF graph debug.
		virtual std::string debugInfo () const;

		// Is the node returning a constant value ? (default = false)
		virtual bool isConstant () const;

		// Derive with respect to node (default = error)
		virtual NodeRef derive (const Node & node);
		virtual bool isDerivable (const Node & node);

	protected:
		// Computation implementation
		virtual void compute () = 0;

		void invalidateRecursively () noexcept;
		void computeRecursively ();

		void makeInvalid () noexcept { isValid_ = false; }
		void makeValid () noexcept { isValid_ = true; }

		// FIXME delete later
		// Only use is in model, but auto creation of Parameters will move to create.
		void appendDependency (NodeRef node);

	private:
		void registerNode (Node * n);
		void unregisterNode (const Node * n);

	protected:
		// TODO small opt vector ?
		NodeRefVec dependencyNodes_{};    // Nodes that we depend on.
		Vector<Node *> dependentNodes_{}; // Nodes that depend on us.

	private:
		bool isValid_{false};

		friend class InternalAccessor;
	};

	/* Node are by convention manipulated as shared_ptr<NodeType>.
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

	// Defers to the Builder<NodeType>::make.
	template <typename NodeType, typename... Args>
	auto makeNode (Args &&... args)
	    -> decltype (Builder<NodeType>::make (std::forward<Args> (args)...)) {
		return Builder<NodeType>::make (std::forward<Args> (args)...);
	}
	// Overload that accepts a NodeRef initializer list (commonly used)
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

	/* Valued node.
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

		// Get updated value
		const T & getValue () {
			this->computeRecursively ();
			return value_;
		}
		// TODO support for // computation... type tags ?

		// Defined as default to enable specialisation
		std::string debugInfo () const override { return Node::debugInfo (); }

	protected:
		T value_;

	private:
		friend class InternalAccessor;
	};

	// Debug info override for double (in DataFlowNumeric.cpp)
	template <> std::string Value<double>::debugInfo () const;
	// overrides for VectorDouble and MatrixDouble are declared in LinearAlgebra.h

	/* Dependency structure description.
	 * These type tags are used to specify compute node dependency types.
	 * This can serve as documentation about what arguments node expect.
	 * Helper functions in DataFlowInternalTemplates.h act depending on these type tags.
	 */
	template <typename T> struct ReductionOfValue {};        // Dynamic sized list of Value<T>
	template <typename... Types> struct FunctionOfValues {}; // Tuple of Value<T0>, Value<T1>, ...

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
