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

// FIXME improve (test with eigen::vec to decide impl)
template <typename T> struct NumericInfo { using Derivable = std::false_type; };
template <> struct NumericInfo<double> {
	using Derivable = std::true_type;
	static double zero () noexcept { return 0.0; }
	static double one () noexcept { return 1.0; }
};

namespace DF {
	// Fwd declaration
	class Node;

	// Error functions
	[[noreturn]] void failureComputeWasCalled (const std::type_info & paramType);
	[[noreturn]] void failureNodeConversion (const std::type_info & handleType, const Node & node);
	[[noreturn]] void failureDerivationNotSupportedForType (const std::type_info & type);

	// Convenient typedefs : Node is supposed to be used as shared_ptr instances.
	using NodeRef = std::shared_ptr<Node>;
	using NodeRefVec = Vector<NodeRef>;

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
		virtual bool isConstant () const { return false; }

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

	/* Valued node.
	 */
	template <typename T> class Value : public Node {
	public:
		// Init deps
		template <typename... Args>
		Value (NodeRefVec deps, Args &&... args)
		    : Node (std::move (deps)), value_ (std::forward<Args> (args)...) {}

		// No deps TODO add type tag to guarantee template choice ?
		template <typename... Args>
		Value (Args &&... args) : Node (), value_ (std::forward<Args> (args)...) {}

		// Access the stored value (no recomputation !)
		const T & value () const noexcept {
			assert (this->isValid ());
			return value_;
		}

		std::string description () const override { return "Value<" + prettyTypeName<T> () + ">"; }

	protected:
		T value_;
	};

	/* Constant value.
	 */
	template <typename T> class Constant : public Value<T> {
	public:
		template <typename... Args>
		Constant (Args &&... args) : Value<T> (std::forward<Args> (args)...) {
			this->makeValid ();
		}

		void compute () override final { failureComputeWasCalled (typeid (Constant<T>)); }

		// Deriving a constant returns 0
		NodeRef derive (const Node &) override final {
			return std::make_shared<Constant<T>> (NumericInfo<T>::zero ());
		}
		bool isConstant () const override final { return true; }

		std::string description () const override final {
			return "Constant<" + prettyTypeName<T> () + ">(" + debug_to_string (this->getValue ()) + ")";
		}
	};

	/* Param node.
	 */
	template <typename T> class Parameter : public Value<T> {
	public:
		template <typename... Args>
		Parameter (Args &&... args) : Value<T> (std::forward<Args> (args)...) {
			this->makeValid ();
		}

		void setValue (T t) noexcept {
			this->invalidate ();
			this->value_ = std::move (t);
			this->makeValid ();
		}

		NodeRef derive (const Node & variable) override final;

		std::string description () const final { return "Parameter<" + prettyTypeName<T> () + ">"; }

	private:
		void compute () override final { failureComputeWasCalled (typeid (Parameter<T>)); }
	};

	// Derivation only make sense for some types (like not for Sequence*)
	// Use template trick to only generate an error for unsupported types
	// FIXME improve too
	template <typename T>
	NodeRef deriveParameterHelper (const Parameter<T> & param, const Node & variable,
	                               std::true_type) {
		auto && v = (&variable == &param) ? NumericInfo<T>::one () : NumericInfo<T>::zero ();
		return std::make_shared<Constant<T>> (v);
	}
	template <typename T>
	NodeRef deriveParameterHelper (const Parameter<T> &, const Node &, std::false_type) {
		failureDerivationNotSupportedForType (typeid (T));
	}
	template <typename T> NodeRef Parameter<T>::derive (const Node & variable) {
		return deriveParameterHelper<T> (*this, variable, typename NumericInfo<T>::Derivable{});
	}

	// Convert handles
	template <typename T, typename U> std::shared_ptr<T> convert (const std::shared_ptr<U> & from) {
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
}
}

#endif // BPP_NEWPHYL_DATAFLOW_H
