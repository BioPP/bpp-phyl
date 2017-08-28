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
  /* TODO reorganize
   * Separate Handle types from node types.
   */

	// Fwd declaration
	class Node;
	template <typename T> class Value;
	template <typename T> class Parameter;

	/* Base Node.
	 * Always contains a node (as do all handle classes.
	 */
	class Node {
	public:
		class Impl;

		template <typename T> Node (const Value<T> & v) noexcept;
		template <typename T> Node (const Parameter<T> & p) noexcept;
		const std::shared_ptr<Impl> & getShared () const noexcept { return pImpl_; }

		template <typename T, typename... Args> static Node create (Args &&... args) {
			return Node (std::make_shared<T> (std::forward<Args> (args)...));
		}
		template <typename T, typename... Args>
		static Node create (std::initializer_list<Node> deps, Args &&... args) {
			return Node (std::make_shared<T> (std::move (deps), std::forward<Args> (args)...));
		}

		Impl & getImpl () const noexcept { return *pImpl_; }

		bool operator== (const Node & other) const noexcept { return pImpl_ == other.pImpl_; }
		bool operator== (const Node::Impl * other) const noexcept { return pImpl_.get () == other; }
		std::size_t hashCode () const noexcept { return std::hash<std::shared_ptr<Impl>>{}(pImpl_); }

	private:
		explicit Node (std::shared_ptr<Impl> p) noexcept : pImpl_ (std::move (p)) {}
		std::shared_ptr<Impl> pImpl_;
	};

	using NodeVec = Vector<Node>;

	class Node::Impl {
	public:
		Impl () = default;
		Impl (const Impl &) = delete;
		Impl (Impl &&) = delete;
		Impl & operator= (const Impl &) = delete;
		Impl & operator= (Impl &&) = delete;

		// Constructor that sets the dependencies
		Impl (NodeVec dependencies) : dependencyNodes_ (std::move (dependencies)) {
			for (auto & n : dependencyNodes_)
				n.getImpl ().registerNode (this);
		}

		virtual ~Impl () {
			for (auto & n : dependencyNodes_)
				n.getImpl ().unregisterNode (this);
		}

		bool isValid () const noexcept { return isValid_; }
		void invalidate () noexcept;

		virtual void compute () = 0;
		void computeRecursively ();

		// Derivation stuff
		virtual Node derive (const Node & variable); // Defaults to error
		virtual bool isConstant () const { return false; }

		// TODO Remove ?
		template <typename F> void foreachDependentNode (F f) const {
			for (auto & p : dependentNodes_)
				f (p);
		}
		template <typename F> void foreachDependencyNode (F f) const {
			for (auto & p : dependencyNodes_)
				f (&p.getImpl ());
		}

		const Vector<Impl *> & dependentNodes () const noexcept { return dependentNodes_; }
		const NodeVec & dependencies () const noexcept { return dependencyNodes_; }

		// Debug information (smaller graph)
		virtual std::string description () const { return "Node"; }

	protected:
		void makeValid () noexcept { isValid_ = true; }

		// FIXME replace ? This should only be used for complex initialisation !
		void appendDependency (Node node) noexcept {
			node.getImpl ().registerNode (this);
			dependencyNodes_.emplace_back (std::move (node));
			invalidate ();
		}

	private:
		void registerNode (Impl * n);
		void unregisterNode (const Impl * n);

	protected:
		// TODO small opt vector ?
		Vector<Impl *> dependentNodes_{}; // Nodes that depend on us.
		NodeVec dependencyNodes_{};       // Nodes that we depend on.

	private:
		bool isValid_{false};
	};

	// Error functions
	[[noreturn]] void failureComputeWasCalled (const std::type_info & paramType);
	[[noreturn]] void failureNodeHandleConversion (const std::type_info & handleType, const Node::Impl & node);
	[[noreturn]] void failureDerivationNotSupportedForType (const std::type_info & type);

	/* Valued node.
	 */
	template <typename T> class Value {
	public:
		class Impl;

		explicit Value (const Node & n);
		Value (const Parameter<T> & p) noexcept;
		const std::shared_ptr<Impl> & getShared () const noexcept { return pImpl_; }

		template <typename U, typename... Args> static Value create (Args &&... args) {
			return Value (std::make_shared<U> (std::forward<Args> (args)...));
		}
		template <typename U, typename... Args>
		static Value create (std::initializer_list<Node> deps, Args &&... args) {
			return Value (std::make_shared<U> (std::move (deps), std::forward<Args> (args)...));
		}

		Impl & getImpl () const noexcept { return *pImpl_; }
		const T & getValue () {
			// The value class is the interface, perform the recomputation.
			pImpl_->computeRecursively ();
			return pImpl_->getValue ();
		}

	private:
		explicit Value (std::shared_ptr<Impl> p) noexcept : pImpl_ (std::move (p)) {}
		std::shared_ptr<Impl> pImpl_;
	};

	template <typename T> class Value<T>::Impl : public Node::Impl {
	public:
		// Init deps
		template <typename... Args>
		Impl (NodeVec deps, Args &&... args)
		    : Node::Impl (std::move (deps)), value_ (std::forward<Args> (args)...) {}

		// No deps TODO add type tag to guarantee template choice ?
		template <typename... Args>
		Impl (Args &&... args) : Node::Impl (), value_ (std::forward<Args> (args)...) {}

		const T & getValue () const noexcept {
			// Implementation, do not recompute (should be done using Node::Impl functions before).
			assert (this->isValid ());
			return value_;
		}

		std::string description () const override { return "Value<" + prettyTypeName<T> () + ">"; }

	protected:
		T value_;
	};

	/* Constant value.
	 */
	template <typename T> class Constant : public Value<T>::Impl {
	public:
		template <typename... Args>
		Constant (Args &&... args) : Value<T>::Impl (std::forward<Args> (args)...) {
			this->makeValid ();
		}

		bool isConstant () const override final { return true; }

		// Create a constant of default value (usually zeroes for must numeric types).
		Node derive (const Node &) override final {
			return Node::create<Constant<T>> (NumericInfo<T>::zero ());
		}

		std::string description () const override final {
			return "Constant<" + prettyTypeName<T> () + ">(" + debug_to_string (this->getValue ()) + ")";
		}

	private:
		void compute () override final { failureComputeWasCalled (typeid (Constant<T>)); }
	};

	/* Param node.
	 */
	template <typename T> class Parameter {
	public:
		class Impl;

		explicit Parameter (const Node & n);
		explicit Parameter (const Value<T> & v);
		const std::shared_ptr<Impl> & getShared () const noexcept { return pImpl_; }

		template <typename... Args> static Parameter create (Args &&... args) {
			return Parameter (std::make_shared<Impl> (std::forward<Args> (args)...));
		}

		Impl & getImpl () const noexcept { return *pImpl_; }
		const T & getValue () const noexcept { return pImpl_->getValue (); }
		void setValue (T t) noexcept { pImpl_->setValue (std::move (t)); }

	private:
		explicit Parameter (std::shared_ptr<Impl> p) noexcept : pImpl_ (std::move (p)) {}
		std::shared_ptr<Impl> pImpl_;
	};

	template <typename T> class Parameter<T>::Impl : public Value<T>::Impl {
	public:
		template <typename... Args>
		Impl (Args &&... args) : Value<T>::Impl (std::forward<Args> (args)...) {
			this->makeValid ();
		}

		void setValue (T t) noexcept {
			this->invalidate ();
			this->value_ = std::move (t);
			this->makeValid ();
		}

		Node derive (const Node & variable) override final;

		std::string description () const final { return "Parameter<" + prettyTypeName<T> () + ">"; }

	private:
		void compute () override final { failureComputeWasCalled (typeid (Parameter<T>)); }
	};

	// Derivation only make sense for some types (like not for Sequence*)
	// Use template trick to only generate an error for unsupported types
  // FIXME improve too
	template <typename T>
	Node deriveParameterHelper (const typename Parameter<T>::Impl * param, const Node & variable,
	                            std::true_type) {
		auto v = (variable == param) ? NumericInfo<T>::one () : NumericInfo<T>::zero ();
		return Node::create<Constant<T>> (v);
	}
	template <typename T>
	Node deriveParameterHelper (const typename Parameter<T>::Impl *, const Node &, std::false_type) {
		failureDerivationNotSupportedForType (typeid (T));
	}
	template <typename T> Node Parameter<T>::Impl::derive (const Node & variable) {
		return deriveParameterHelper<T> (this, variable, typename NumericInfo<T>::Derivable{});
	}

	/* Conversion constructors
	 */
	template <typename T> Node::Node (const Value<T> & v) noexcept : pImpl_ (v.getShared ()) {}
	template <typename T> Node::Node (const Parameter<T> & p) noexcept : pImpl_ (p.getShared ()) {}

	template <typename T>
	Value<T>::Value (const Node & n) : pImpl_ (std::dynamic_pointer_cast<Impl> (n.getShared ())) {
		if (!pImpl_)
			failureNodeHandleConversion (typeid (Value<T>), n.getImpl ());
	}
	template <typename T>
	Value<T>::Value (const Parameter<T> & p) noexcept : pImpl_ (p.getShared ()) {}

	template <typename T>
	Parameter<T>::Parameter (const Node & n)
	    : pImpl_ (std::dynamic_pointer_cast<Impl> (n.getShared ())) {
		if (!pImpl_)
			failureNodeHandleConversion (typeid (Parameter<T>), n.getImpl ());
	}
	template <typename T>
	Parameter<T>::Parameter (const Value<T> & v)
	    : pImpl_ (std::dynamic_pointer_cast<Impl> (v.getShared ())) {
		if (!pImpl_)
			failureNodeHandleConversion (typeid (Parameter<T>), v.getImpl ());
	}

	/* Node value access.
	 */
	template <typename T> bool isValueNode (const Node & n) noexcept {
		return dynamic_cast<const typename Value<T>::Impl *> (&n.getImpl ()) != nullptr;
	}
	template <typename T> const T & getValueUnsafe (const Node & n) noexcept {
		assert (isValueNode<T> (n));
		return static_cast<const typename Value<T>::Impl &> (n.getImpl ()).getValue ();
	}
}
}

namespace std {
template <> struct hash<bpp::DF::Node> {
	using argument_type = bpp::DF::Node;
	using result_type = std::size_t;
	result_type operator() (const argument_type & node) const { return node.hashCode (); }
};
}

#endif // BPP_NEWPHYL_DATAFLOW_H
