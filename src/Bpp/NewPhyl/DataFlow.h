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

#include <Bpp/NewPhyl/Debug.h> // description
#include <Bpp/NewPhyl/Vector.h>
#include <cassert>
#include <memory>
#include <string> // description
#include <typeinfo>
#include <utility>

namespace bpp {
namespace DF {
	// Fwd declaration
	class Node;
  template<typename T> class Value;
  struct NumericProperties; // In DataFlowNumeric.h

	// Error functions
	[[noreturn]] void failureComputeWasCalled (const std::type_info & paramType);
	[[noreturn]] void failureNodeConversion (const std::type_info & handleType, const Node & node);

	// Convenient typedefs : Node is supposed to be used as shared_ptr instances.
	using NodeRef = std::shared_ptr<Node>;
	using NodeRefVec = Vector<NodeRef>;
	template <typename T> using ValueRef = std::shared_ptr<Value<T>>;

	/* Base Node class.
	 * Abstract : compute() needs to be defined to the actual computation.
	 * TODO determine what to remove from the class API (computeRecursively ?)
	 * → depends if Node is high level (safe) or low level
	 * → choose invariants carefully : tend to think these are low level classes
	 */
	class Node {
	public:
		Node () = default;
		Node (const Node &) = delete;
		Node (Node &&) = delete;
		Node & operator= (const Node &) = delete;
		Node & operator= (Node &&) = delete;
		Node (NodeRefVec dependencies); // Sets dependencies + register
		virtual ~Node ();               // Deregisters from dependencies

		// Accessors
		const Vector<Node *> & dependentNodes () const noexcept { return dependentNodes_; }
		const NodeRefVec & dependencies () const noexcept { return dependencyNodes_; }
		bool isValid () const noexcept { return isValid_; }

		// Computation
		void invalidate () noexcept;
		virtual void compute () = 0;
		void computeRecursively ();

		// Derivation stuff
		virtual NodeRef derive (const Node & variable); // Defaults to error
    virtual void numericProperties (NumericProperties &) const;
    bool isConstant () const;

		// Debug information (smaller graph)
		virtual std::string description () const { return "Node"; }

	protected:
		void makeValid () noexcept { isValid_ = true; }

		// FIXME replace ? Only for complex init
		void appendDependency (NodeRef node);

	private:
		void registerNode (Node * n);
		void unregisterNode (const Node * n);

	protected:
		// TODO small opt vector ?
		Vector<Node *> dependentNodes_{}; // Nodes that depend on us.
		NodeRefVec dependencyNodes_{};    // Nodes that we depend on.

	private:
		bool isValid_{false};
	};

	// Create node with make shared.
	template <typename T, typename... Args> std::shared_ptr<T> createNode (Args &&... args) {
		return std::make_shared<T> (std::forward<Args> (args)...);
	}
	template <typename T, typename... Args>
	std::shared_ptr<T> createNode (std::initializer_list<NodeRef> ilist, Args &&... args) {
		return std::make_shared<T> (std::move (ilist), std::forward<Args> (args)...);
	}

	/* Valued node.
	 */
	struct NoDependencyTag {
		constexpr NoDependencyTag () = default;
	};
	constexpr NoDependencyTag noDependency{};

	template <typename T> class Value : public Node {
	public:
		// Init deps
		template <typename... Args>
		Value (NodeRefVec deps, Args &&... args)
		    : Node (std::move (deps)), value_ (std::forward<Args> (args)...) {}

		// Without dependencies
		template <typename... Args>
		Value (NoDependencyTag, Args &&... args) : Node (), value_ (std::forward<Args> (args)...) {}

		// Access the stored value (no recomputation !)
		const T & value () const noexcept {
			assert (this->isValid ());
			return value_;
		}

		std::string description () const override { return "Value<" + prettyTypeName<T> () + ">"; }

	protected:
		T value_;
	};

	// Convert handles with check
	template <typename T, typename U>
	std::shared_ptr<T> convertRef (const std::shared_ptr<U> & from) {
		auto p = std::dynamic_pointer_cast<T> (from);
		if (!p)
			failureNodeConversion (typeid (T), *from);
		return p;
	}

	/* Node value access.
	 * TODO improve
	 */
	template <typename T> bool isValueNode (const Node & n) noexcept {
		return dynamic_cast<const Value<T> *> (&n) != nullptr;
	}
	template <typename T> const T & accessValueUnsafe (const Node & n) noexcept {
		assert (isValueNode<T> (n));
		return static_cast<const Value<T> &> (n).value ();
	}

	// Get Value<T> value with recompute TODO improve
	template <typename T> const T & getUpToDateValue (Value<T> & node) {
		node.computeRecursively ();
		return node.value ();
	}
	template <typename T>
	auto getUpToDateValue (const std::shared_ptr<T> & node) -> decltype (getUpToDateValue (*node)) {
		return getUpToDateValue (*node);
	}
}
}

#endif // BPP_NEWPHYL_DATAFLOW_H
